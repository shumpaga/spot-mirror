// -*- coding: utf-8 -*-
// Copyright (C) 2015, 2016, 2017, 2018 Laboratoire de Recherche et
// Developpement de l'Epita
//
// This file is part of Spot, a model checking library.
//
// Spot is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Spot is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <atomic>
#include <thread>
#include <vector>

#include <spot/bricks/brick-hashset>
#include <spot/kripke/kripke.hh>
#include <spot/misc/common.hh>
#include <spot/misc/fixpool.hh>
#include <spot/misc/timer.hh>
#include <spot/twacube/twacube.hh>

namespace spot
{
  /// \brief This object is returned by the algorithm below
  struct SPOT_API cndfs_stats
    {
     unsigned states;            ///< \brief Number of states visited
     unsigned transitions;       ///< \brief Number of transitions visited
     unsigned instack_dfs;       ///< \brief Maximum DFS stack
     bool is_empty;              ///< \brief Is the model empty
     unsigned walltime;          ///< \brief Walltime for this thread in ms
    };

  template<typename State, typename SuccIterator,
           typename StateHash, typename StateEqual>
  class swarmed_cndfs
  {
    struct local_colors
    {
      bool cyan;
      bool is_in_Rp;
    };

    /// \brief The colors of a state
    struct cndfs_colors
    {
      std::atomic<bool> blue;
      std::atomic<bool> red;
      local_colors l[1];
    };

    struct product_state
    {
      State st_kripke;
      unsigned st_prop;
      cndfs_colors* colors;
    };

    /// \brief The hasher for the previous state.
    struct state_hasher
    {
      state_hasher(const product_state&)
      { }

      state_hasher() = default;

      brick::hash::hash128_t
      hash(const product_state& lhs) const
      {
        StateHash hash;
        // Not modulo 31 according to brick::hashset specifications.
        unsigned u = hash(lhs.st_kripke) % (1<<30);
        u = wang32_hash(lhs.st_prop) ^ u;
        u = u % (1<<30);
        return {u, u};
      }

      bool equal(const product_state& lhs,
                 const product_state& rhs) const
      {
        StateEqual equal;
        return (lhs.st_prop == rhs.st_prop)
          && equal(lhs.st_kripke, rhs.st_kripke);
      }
    };

    struct todo__element
    {
      product_state st;
      SuccIterator* it_kripke;
      std::shared_ptr<trans_index> it_prop;
      bool from_accepting;
    };

  public:

    ///< \brief Shortcut to ease shared map manipulation
    using shared_map = brick::hashset::FastConcurrent <product_state,
                                                       state_hasher>;

    swarmed_cndfs(kripkecube<State, SuccIterator>& sys, twacube_ptr twa,
                  shared_map& map, unsigned tid, std::atomic<bool>& stop):
      sys_(sys), twa_(twa), tid_(tid), map_(map),
      nb_th_(std::thread::hardware_concurrency()),
      p_colors_(sizeof(cndfs_colors) +
                sizeof(local_colors)*(std::thread::hardware_concurrency() - 1)),
      stop_(stop)
    {
      SPOT_ASSERT(is_a_kripkecube(sys));
    }

    virtual ~swarmed_cndfs()
    {
      while (!todo_blue_.empty())
        {
          sys_.recycle(todo_blue_.back().it_kripke, tid_);
          todo_blue_.pop_back();
        }
      while (!todo_red_.empty())
        {
          sys_.recycle(todo_red_.back().it_kripke, tid_);
          todo_red_.pop_back();
        }
    }

    void setup()
    {
      tm_.start("DFS thread " + std::to_string(tid_));
    }

    std::pair<bool,  product_state>
    push_blue(product_state s, bool from_accepting)
    {
      cndfs_colors* c = (cndfs_colors*) p_colors_.allocate();
      c->red = false;
      c->blue = false;
      for (unsigned i = 0; i < nb_th_; ++i)
        {
          c->l[i].cyan = false;
          c->l[i].is_in_Rp = false;
        }

      s.colors = c;

      // Try to insert the new state in the shared map.
      auto it = map_.insert(s);
      bool b = it.isnew();

      // Insertion failed, delete element
      // FIXME Should we add a local cache to avoid useless allocations?
      if (!b)
        {
          p_colors_.deallocate(c);
          bool blue = it->colors->blue.load();
          bool cyan = it->colors->l[tid_].cyan;
          if (blue || cyan)
            return {false, *it};
        }

      // Mark state as visited.
      it->colors->l[tid_].cyan = true;
      ++states_;
      todo_blue_.push_back({*it,
                            sys_.succ(it->st_kripke, tid_),
                            twa_->succ(it->st_prop),
                            from_accepting});
      return {true, *it};
    }

    bool pop_blue()
    {
      // Track maximum dfs size
      dfs_ = todo_blue_.size()  > dfs_ ? todo_blue_.size() : dfs_;

      todo_blue_.back().st.colors->l[tid_].cyan = false;
      sys_.recycle(todo_blue_.back().it_kripke, tid_);
      todo_blue_.pop_back();
      return true;
    }

    std::pair<bool, product_state>
    push_red(product_state s, bool ignore_cyan)
    {
      // Try to insert the new state in the shared map.
      auto it = map_.insert(s);
      bool b = it.isnew();

      SPOT_ASSERT(!b); // should never be new in a red DFS
      bool red = it->colors->red.load();
      bool cyan = it->colors->l[tid_].cyan;
      bool in_Rp = it->colors->l[tid_].is_in_Rp;
      if (red || (cyan && !ignore_cyan) || in_Rp)
        return {false, *it}; // couldn't insert

      // Mark state as visited.
      it->colors->l[tid_].is_in_Rp = true;
      Rp_.push_back(*it);
      ++states_;
      todo_red_.push_back({*it,
                           sys_.succ(it->st_kripke, tid_),
                           twa_->succ(it->st_prop),
                           false});
      return {true, *it};
    }

    bool pop_red()
    {
      // Track maximum dfs size
      dfs_ = todo_blue_.size() + todo_red_.size() > dfs_ ?
        todo_blue_.size() + todo_red_.size() : dfs_;


      sys_.recycle(todo_red_.back().it_kripke, tid_);
      todo_red_.pop_back();
      return true;
    }

    void finalize()
    {
      tm_.stop("DFS thread " + std::to_string(tid_));
    }

    unsigned states()
    {
      return states_;
    }

    unsigned transitions()
    {
      return transitions_;
    }

    void run()
    {
      setup();
      blue_dfs();
      finalize();
    }

    void blue_dfs()
    {
      product_state initial = {sys_.initial(tid_),
                               twa_->get_initial(),
                               nullptr};
      if (!push_blue(initial, false).first)
        return;

      // Property automaton has only one state
      if (todo_blue_.back().it_prop->done())
        return;

      forward_iterators(todo_blue_, true);

      while (!todo_blue_.empty() && !stop_.load(std::memory_order_relaxed))
        {
          auto current = todo_blue_.back();

          if (!current.it_kripke->done())
            {
              ++transitions_;
              product_state s = {
                                 current.it_kripke->state(),
                                 twa_->trans_storage(current.it_prop, tid_).dst,
                                 nullptr
              };

              bool acc = (bool) twa_->trans_storage(current.it_prop, tid_).acc_;
              forward_iterators(todo_blue_, false);

              auto tmp = push_blue(s, acc);
              if (tmp.first)
                forward_iterators(todo_blue_, true);
              else if (acc)
                {
                  // The state cyan and we can reach it throught an
                  // accepting transition, a accepting cycle has been
                  // found without launching a red dfs
                  if (tmp.second.colors->l[tid_].cyan)
                    {
                      cycle_start_ = s;
                      is_empty_ = false;
                      stop_.store(true);
                      return;
                    }

                  SPOT_ASSERT(tmp.second.colors->blue);

                  red_dfs(s);
                  if (!is_empty_)
                    return;
                  post_red_dfs();
                }
            }
          else
            {
              current.st.colors->blue.store(true);

              // backtracked an accepting transition; launch red DFS
              if (current.from_accepting)
                {
                  red_dfs(todo_blue_.back().st);
                  if (!is_empty_)
                    return;
                  post_red_dfs();
                }

              pop_blue();
            }
        }
    }

    void post_red_dfs()
    {
      for (product_state& s : Rp_acc_)
        {
          while (s.colors->red.load() && !stop_.load())
            {
              // await
            }
        }
      for (product_state& s : Rp_)
        {
          s.colors->red.store(true);
          s.colors->l[tid_].is_in_Rp = false; // empty Rp
        }

      Rp_.clear();
      Rp_acc_.clear();
    }

    void red_dfs(product_state initial)
    {
      auto init_push = push_red(initial, true);
      SPOT_ASSERT(init_push.second.colors->blue);

      if (!init_push.first)
        return;

      forward_iterators(todo_red_, true);

      while (!todo_red_.empty() && !stop_.load(std::memory_order_relaxed))
        {
          auto current = todo_red_.back();

          if (!current.it_kripke->done())
            {
              ++transitions_;
              product_state s = {
                                 current.it_kripke->state(),
                                 twa_->trans_storage(current.it_prop, tid_).dst,
                                 nullptr
              };
              bool acc = (bool) twa_->trans_storage(current.it_prop, tid_).acc_;
              forward_iterators(todo_red_, false);

              auto res = push_red(s, false);
              if (res.first) // could push properly
                {
                  forward_iterators(todo_red_, true);

                  SPOT_ASSERT(res.second.colors->blue);

                  // The transition is accepting, we want to keep
                  // track of this state
                  if (acc)
                    {
                      // Do not insert twice a state
                      bool found = false;
                      for (auto& st: Rp_acc_)
                        {
                          if (st.colors == res.second.colors)
                            {
                              found = true;
                              break;
                            }
                        }
                      if (!found)
                        Rp_acc_.push_back(Rp_.back());
                    }
                }
              else
                {
                  if (res.second.colors->l[tid_].cyan)
                    {
                      // color pointers are unique to each element,
                      // comparing them is equivalent (but faster) to comparing
                      // st_kripke and st_prop individually.
                      if (init_push.second.colors == res.second.colors && !acc)
                        continue;
                      is_empty_ = false;
                      cycle_start_ = s;
                      stop_.store(true);
                      return;
                    }
                  else if (acc && res.second.colors->l[tid_].is_in_Rp)
                    {
                      auto it = map_.insert(s);
                      Rp_acc_.push_back(*it);
                    }
                }
            }
          else
            {
              pop_red();
            }
        }
    }

    std::string trace()
    {
      SPOT_ASSERT(!is_empty());
      StateEqual equal;
      auto state_equal = [equal](product_state a, product_state b)
                         {
                           return a.st_prop == b.st_prop
                             && equal(a.st_kripke, b.st_kripke);
                         };

      std::string res = "Prefix:\n";

      auto it = todo_blue_.begin();
      while (it != todo_blue_.end())
        {
          if (state_equal(it->st, cycle_start_))
            break;
          res += "  " + std::to_string(it->st.st_prop)
            + "*" + sys_.to_string(it->st.st_kripke) + "\n";
          ++it;
        }

      res += "Cycle:\n";
      while (it != todo_blue_.end())
        {
          res += "  " + std::to_string(it->st.st_prop)
            + "*" + sys_.to_string(it->st.st_kripke) + "\n";
          ++it;
        }

      if (!todo_red_.empty())
        {
          it = todo_red_.begin() + 1; // skip first element, also in blue
          while (it != todo_red_.end())
            {
              res += "  " + std::to_string(it->st.st_prop)
                + "*" + sys_.to_string(it->st.st_kripke) + "\n";
              ++it;
            }
        }
      res += "  " + std::to_string(cycle_start_.st_prop)
        + "*" + sys_.to_string(cycle_start_.st_kripke) + "\n";

      return res;
    }

    bool is_empty()
    {
      return is_empty_;
    }

    unsigned walltime()
    {
      return tm_.timer("DFS thread " + std::to_string(tid_)).walltime();
    }

    cndfs_stats stats()
    {
      return {states(), transitions(), dfs_, is_empty(), walltime()};
    }

  protected:
    void forward_iterators(std::vector<todo__element>& todo, bool just_pushed)
    {
      SPOT_ASSERT(!todo.empty());

      auto top = todo.back();

      SPOT_ASSERT(!(top.it_prop->done() &&
                    top.it_kripke->done()));

      // Sometimes kripke state may have no successors.
      if (top.it_kripke->done())
        return;

      // The state has just been push and the 2 iterators intersect.
      // There is no need to move iterators forward.
      SPOT_ASSERT(!(top.it_prop->done()));
      if (just_pushed && twa_->get_cubeset()
          .intersect(twa_->trans_data(top.it_prop, tid_).cube_,
                     top.it_kripke->condition()))
        return;

      // Otherwise we have to compute the next valid successor (if it exits).
      // This requires two loops. The most inner one is for the twacube since
      // its costless
      if (top.it_prop->done())
        top.it_prop->reset();
      else
        top.it_prop->next();

      while (!top.it_kripke->done())
        {
          while (!top.it_prop->done())
            {
              if (twa_->get_cubeset()
                  .intersect(twa_->trans_data(top.it_prop, tid_).cube_,
                             top.it_kripke->condition()))
                return;
              top.it_prop->next();
            }
          top.it_prop->reset();
          top.it_kripke->next();
        }
    }

  private:
    kripkecube<State, SuccIterator>& sys_;
    twacube_ptr twa_;
    std::vector<todo__element> todo_blue_;
    std::vector<todo__element> todo_red_;
    unsigned transitions_ = 0;         ///< \brief Number of transitions
    unsigned tid_;                     ///< \brief Thread's current ID
    shared_map map_;                       ///< \brief Map shared by threads
    spot::timer_map tm_;                   ///< \brief Time execution
    unsigned states_ = 0;                  ///< \brief Number of states
    unsigned dfs_ = 0;                     ///< \brief Maximum DFS stack size
    /// \brief Maximum number of threads that can be handled by this algorithm
    unsigned nb_th_ = 0;
    fixed_size_pool<pool_type::Unsafe> p_colors_;
    bool is_empty_ = true;               ///< \brief Accepting cycle detected?
    std::atomic<bool>& stop_;              ///< \brief Stop-the-world boolean
    std::vector<product_state> Rp_;
    std::vector<product_state> Rp_acc_;
    product_state cycle_start_;
  };
}
