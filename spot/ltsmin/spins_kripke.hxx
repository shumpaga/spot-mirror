// -*- coding: utf-8 -*-
// Copyright (C) 2017, 2018 Laboratoire de Recherche et Développement de
// l'Epita (LRDE)
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

#include <algorithm>
#include <bricks/brick-hash>
#include <bricks/brick-hashset>
#include <spot/ltsmin/spins_interface.hh>
#include <spot/misc/fixpool.hh>
#include <spot/misc/mspool.hh>
#include <spot/misc/intvcomp.hh>
#include <spot/misc/intvcmp2.hh>
#include <spot/twacube/cube.hh>

namespace spot
{
  cspins_state_manager::cspins_state_manager(unsigned int state_size,
                                             int compress)
    : p_((state_size+2)*sizeof(int)), compress_(compress),
      /* reserve one integer for the hash value and size */
      state_size_(state_size),
      fn_compress_(compress == 0 ? nullptr
                   : compress == 1 ? int_array_array_compress
                   : int_array_array_compress2),
      fn_decompress_(compress == 0 ? nullptr
                     : compress == 1 ? int_array_array_decompress
                     : int_array_array_decompress2)
  { }

  int* cspins_state_manager::unbox_state(cspins_state s) const
  {
    return s+2;
  }

  // cmp is the area we can use to compute the compressed rep of dst.
  cspins_state
  cspins_state_manager::alloc_setup(int *dst, int* cmp, size_t cmpsize)
  {
    cspins_state out = nullptr;
    size_t size = state_size_;
    int* ref = dst;
    if (compress_)
      {
        size_t csize = cmpsize;
        fn_compress_(dst, state_size_, cmp, csize);
        ref = cmp;
        size = csize;
        out = (cspins_state) msp_.allocate((size+2)*sizeof(int));
      }
    else
      out = (cspins_state) p_.allocate();
    int hash_value = 0;
    memcpy(unbox_state(out), ref, size * sizeof(int));
    for (unsigned int i = 0; i < state_size_; ++i)
      hash_value = wang32_hash(hash_value ^ dst[i]);
    out[0] = hash_value;
    out[1] = size;
    return out;
  }

  void cspins_state_manager::decompress(cspins_state s, int* uncompressed,
                                        unsigned size) const
  {
    fn_decompress_(s+2, s[1], uncompressed, size);
  }

  void cspins_state_manager::dealloc(cspins_state s)
  {
    if (compress_)
      msp_.deallocate(s, (s[1]+2)*sizeof(int));
    else
      p_.deallocate(s);
  }

  unsigned int cspins_state_manager::size() const
  {
    return state_size_;
  }


  // This class provides an iterator over the successors of a state.
  // All successors are computed once when an iterator is recycled or
  // created.
  cspins_iterator::cspins_iterator(cspins_state s,
                                   const spot::spins_interface* d,
                                   cspins_state_manager& manager,
                                   inner_callback_parameters& inner,
                                   cube cond,
                                   bool compress,
                                   bool selfloopize,
                                   cubeset& cubeset,
                                   int dead_idx, unsigned tid,
                                   bool use_por,
                                   porinfos* porinfos,
                                   std::vector<bool>* reduced)
  :  current_(0), cond_(cond), tid_(tid), use_por_(use_por)
  {
    successors_.reserve(10);
    inner.manager = &manager;
    inner.succ = &successors_;
    inner.compress = compress;
    inner.selfloopize = selfloopize;
    int* ref = s;

    // Local variable since we only need it to compute
    // reduced set
    std::vector<int> transitions_id;
    inner.transitions_id = &transitions_id;

    if (compress)
      // already filled by compute_condition
      ref = inner.uncompressed_;

    int n = d->get_successors
      (nullptr, manager.unbox_state(ref),
       [](void* arg, transition_info_t* t, int *dst){
        inner_callback_parameters* inner =
        static_cast<inner_callback_parameters*>(arg);
        cspins_state s =
        inner->manager->alloc_setup(dst, inner->compressed_,
                                    inner->manager->size() * 2);
        inner->succ->push_back(s);
        inner->transitions_id->push_back(t->group);
       },
       &inner);
    if (!n && selfloopize)
      {
        successors_.push_back(s);
        if (dead_idx != -1)
          cubeset.set_true_var(cond, dead_idx);
      }

    if (use_por_) // FIXME should be constexpr!
      setup_por_iterator(s, porinfos, transitions_id, reduced);
  }

  void cspins_iterator::recycle(cspins_state s,
                                const spot::spins_interface* d,
                                cspins_state_manager& manager,
                                inner_callback_parameters& inner,
                                cube cond,
                                bool compress,
                                bool selfloopize,
                                cubeset& cubeset, int dead_idx,
                                unsigned tid,
                                bool use_por,
                                porinfos* porinfos,
                                std::vector<bool>* reduced)
  {
    tid_ = tid;
    cond_ = cond;
    current_ = 0;
    // Constant time since int* is is_trivially_destructible
    successors_.clear();
    inner.manager = &manager;
    inner.succ = &successors_;
    inner.compress = compress;
    inner.selfloopize = selfloopize;

    // Local variable since we only need it to compute
    // reduced set
    std::vector<int> transitions_id;
    inner.transitions_id = &transitions_id;
    // Reset POR variables
    fireall_ = false;
    first_pass_ = true;
    use_por_ = use_por;

    int* ref  = s;

    if (compress)
      // Already filled by compute_condition
      ref = inner.uncompressed_;

    int n = d->get_successors
      (nullptr, manager.unbox_state(ref),
       [](void* arg, transition_info_t* t, int *dst){
        inner_callback_parameters* inner =
        static_cast<inner_callback_parameters*>(arg);
        cspins_state s =
        inner->manager->alloc_setup(dst, inner->compressed_,
                                    inner->manager->size() * 2);
        inner->succ->push_back(s);
        inner->transitions_id->push_back(t->group);
       },
       &inner);
    if (!n && selfloopize)
      {
        successors_.push_back(s);
        if (dead_idx != -1)
          cubeset.set_true_var(cond, dead_idx);
      }

    if (use_por_) // FIXME should be constexpr!
      setup_por_iterator(s, porinfos, transitions_id, reduced);
  }

  cspins_iterator::~cspins_iterator()
  {
    // Do not release successors states, the manager
    // will do it on time.
  }

  void cspins_iterator::next()
  {
    if (!use_por_)
      ++current_;
    else
      next_por();
  }

  bool cspins_iterator::done() const
  {
    return current_ >= successors_.size();
  }

  cspins_state cspins_iterator::state() const
  {
    SPOT_ASSERT(successors_.size());
    if (SPOT_UNLIKELY(!tid_))
      return successors_[current_];
    return  successors_[compute_index()];
  }

  cube cspins_iterator::condition() const
  {
    return cond_;
  }

  void cspins_iterator::fireall()
  {
    // Fireall has been done after all successors have
    // been visited
    if (done() && use_por_ && !fireall_)
      {
        first_pass_ = false;
        fireall_ = true;
        current_ = 0;           // reset iterator

        // Check if current_ is referencing an (enabled\reduced) transition
        if ((tid_ == 0 && !reduced_[current_]) ||
            (tid_ != 0 && !reduced_[compute_index()]))
          return;
        next_por();
      }
    else
      fireall_ = true;

  }

  bool cspins_iterator::naturally_expanded() const
  {
    return !use_por_;
  }

  const std::vector<bool>& cspins_iterator::reduced() const
  {
    return reduced_;
  }

  void cspins_iterator::setup_por_iterator(cspins_state s, porinfos* porinfos,
                                           std::vector<int>& transitions_id,
                                           std::vector<bool>* red)
  {
    SPOT_ASSERT(use_por_);

    // Compute the reduced set for s
    if (red == nullptr)
      reduced_ = porinfos->compute_reduced_set(transitions_id, s);
    else
      reduced_ = *red;

    // Check wether the state is naturally enabled.
    if (std::find(std::begin(reduced_), std::end(reduced_), false)
        == std::end(reduced_))
      {
        use_por_ = false;
        return;
      }

    SPOT_ASSERT(successors_.size());

    // Now we have to update current to the first valid position.
    // distinguish tid = 0 from other.
    if ((tid_ == 0 && reduced_[current_]) ||
        (tid_ != 0 && reduced_[compute_index()]))
      return;
    next_por();
  }

  unsigned cspins_iterator::compute_index() const
  {
    unsigned long long c = current_ + 1;
    unsigned long long p = primes[tid_];
    unsigned long long s = successors_.size();
    return (unsigned)  ((c*p) % s);
  }


  void cspins_iterator::next_por()
  {
    SPOT_ASSERT(use_por_);

    // We are walking the iterator during the first phase.
    if (first_pass_)
      {
        // tid_ == 0, means the "default" order
        if (tid_ == 0)
          {
            do
              {
                ++current_;
              }
            while (!done() && !reduced_[current_]);
          }
        else
          {
            do
              {
                ++current_;
                SPOT_ASSERT(successors_.size());
              }
            while (!done() && !reduced_[compute_index()]);
          }

        // We have reach some enabled transitions in the reduced set,
        // so we can stop.
        if (!done())
          return;

        // Otherwise all the transitions in the reduced set have been explored
        first_pass_ = false;

        if (!fireall_)
          return;

        // Reset the iterator over transitions
        current_ = 0;

        SPOT_ASSERT(successors_.size());

        // Check if current_ is referencing an (enabled\reduced) transition
        if ((tid_ == 0 && !reduced_[current_]) ||
            (tid_ != 0 && !reduced_[compute_index()]))
          return;
      }

    // All the transitions in the reduced set have been explored BUT
    // we asked for expansion, i.e. now we ignore transitions all in
    // reduced and only consider those that are not in reduced
    if (!first_pass_ && fireall_)
      {
        // tid_ == 0, means the "default" order
        if (tid_ == 0)
          {
            do
              {
                ++current_;
              }
            while (!done() && reduced_[current_]);
          }
        else
          {
            do
              {
                ++current_;
                SPOT_ASSERT(successors_.size());
              }
            while (!done() && reduced_[compute_index()]);
          }
      }
  }

  kripkecube<cspins_state, cspins_iterator>
  ::kripkecube (spins_interface_ptr sip,
                bool compress,
                const atomic_prop_set* visible_aps,
                bool selfloopize, std::string dead_prop,
                unsigned int nb_threads, bool use_por,
                const porinfos_options& por_opt)
    : sip_(sip), d_(sip.get()),
      compress_(compress), cubeset_(visible_aps->size()),
      selfloopize_(selfloopize), aps_(visible_aps),
      nb_threads_(nb_threads), use_por_(use_por)
  {
    manager_ = static_cast<cspins_state_manager*>
      (::operator new(sizeof(cspins_state_manager) * nb_threads));
    inner_ = new inner_callback_parameters[nb_threads_];
    for (unsigned i = 0; i < nb_threads_; ++i)
      {
        recycle_.push_back(std::vector<cspins_iterator*>());
        recycle_.back().reserve(2000000);
        new (&manager_[i])
          cspins_state_manager(d_->get_state_size(), compress);
        inner_[i].compressed_ = new int[d_->get_state_size() * 2];
        inner_[i].uncompressed_ = new int[d_->get_state_size()+30];
      }
    dead_idx_ = -1;
    match_aps(visible_aps, dead_prop);

    porinfos_options por_opt_with_ap_set = porinfos_options(por_opt.method,
                                             por_opt.ltl_condition, &pset_);

    porinfos_ = new porinfos(d_, por_opt_with_ap_set);
  }

  kripkecube<cspins_state, cspins_iterator>::~kripkecube()
  {
    for (auto& i: recycle_)
      {
        for (auto& j: i)
          {
            cubeset_.release(j->condition());
            delete j;
          }
      }

    for (unsigned i = 0; i < nb_threads_; ++i)
      {
        manager_[i].~cspins_state_manager();
        delete[] inner_[i].compressed_;
        delete[] inner_[i].uncompressed_;
      }
    ::operator delete(manager_);
    delete[] inner_;
  }

  cspins_state kripkecube<cspins_state, cspins_iterator>::initial(unsigned tid)
  {
    d_->get_initial_state(inner_[tid].uncompressed_);
    return manager_[tid].alloc_setup(inner_[tid].uncompressed_,
                                     inner_[tid].compressed_,
                                     manager_[tid].size() * 2);
  }

  std::string
  kripkecube<cspins_state, cspins_iterator>::to_string(const cspins_state s,
                                                       unsigned tid) const
  {
    std::string res = "";
    cspins_state out = manager_[tid].unbox_state(s);
    cspins_state ref = out;
    if (compress_)
      {
        manager_[tid].decompress(s, inner_[tid].uncompressed_,
                                 manager_[tid].size());
        ref = inner_[tid].uncompressed_;
      }
    for (int i = 0; i < d_->get_state_size(); ++i)
      res += (d_->get_state_variable_name(i)) +
        ("=" + std::to_string(ref[i])) + ",";
    res.pop_back();
    return res;
  }

  cspins_iterator*
  kripkecube<cspins_state, cspins_iterator>::succ(const cspins_state s,
                                                  unsigned tid,
                                                  std::vector<bool>* reduced)
  {
    if (SPOT_LIKELY(!recycle_[tid].empty()))
      {
        auto tmp  = recycle_[tid].back();
        recycle_[tid].pop_back();
        compute_condition(tmp->condition(), s, tid);
        tmp->recycle(s, d_, manager_[tid], inner_[tid],
                     tmp->condition(), compress_, selfloopize_,
                     cubeset_, dead_idx_, tid, use_por_, porinfos_,
                     reduced);
        return tmp;
      }
    cube cond = cubeset_.alloc();
    compute_condition(cond, s, tid);
    return new cspins_iterator(s, d_, manager_[tid], inner_[tid],
                               cond, compress_, selfloopize_,
                               cubeset_, dead_idx_, tid, use_por_, porinfos_,
                               reduced);
  }

  void
  kripkecube<cspins_state, cspins_iterator>::recycle(cspins_iterator* it,
                                                     unsigned tid)
  {
    recycle_[tid].push_back(it);
  }

  const std::vector<std::string>
  kripkecube<cspins_state, cspins_iterator>::get_ap_str()
  {
    return aps_str_;
  }


  const atomic_prop_set*
  kripkecube<cspins_state, cspins_iterator>::get_ap()
  {
    return aps_;
  }

  unsigned kripkecube<cspins_state, cspins_iterator>::get_threads()
  {
    return nb_threads_;
  }

  void
  kripkecube<cspins_state, cspins_iterator>::compute_condition
  (cube c, cspins_state s, unsigned tid)
  {
    int i = -1;
    int *vars = manager_[tid].unbox_state(s);

    // State is compressed, uncompress it
    if (compress_)
      {
        manager_[tid].decompress(s, inner_[tid].uncompressed_+2,
                                 manager_[tid].size());
        vars = inner_[tid].uncompressed_ + 2;
      }

    for (auto& ap: pset_)
      {
        ++i;
        bool cond = false;
        switch (ap.op)
          {
          case kripke_cube::OP_EQ_VAR:
            cond = (ap.lval.val == vars[ap.rval.var]);
            break;
          case kripke_cube::OP_NE_VAR:
            cond = (ap.lval.val != vars[ap.rval.var]);
            break;
          case kripke_cube::OP_LT_VAR:
            cond = (ap.lval.val < vars[ap.rval.var]);
            break;
          case kripke_cube::OP_GT_VAR:
            cond = (ap.lval.val > vars[ap.rval.var]);
            break;
          case kripke_cube::OP_LE_VAR:
            cond = (ap.lval.val <= vars[ap.rval.var]);
            break;
          case kripke_cube::OP_GE_VAR:
            cond = (ap.lval.val >= vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_EQ:
            cond = (vars[ap.lval.var] == ap.rval.val);
            break;
          case kripke_cube::VAR_OP_NE:
            cond = (vars[ap.lval.var] != ap.rval.val);
            break;
          case kripke_cube::VAR_OP_LT:
            cond = (vars[ap.lval.var] < ap.rval.val);
            break;
          case kripke_cube::VAR_OP_GT:
            cond = (vars[ap.lval.var] > ap.rval.val);
            break;
          case kripke_cube::VAR_OP_LE:
            cond = (vars[ap.lval.var] <= ap.rval.val);
            break;
          case kripke_cube::VAR_OP_GE:
            cond = (vars[ap.lval.var] >= ap.rval.val);
            break;
          case kripke_cube::VAR_OP_EQ_VAR:
            cond = (vars[ap.lval.var] == vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_NE_VAR:
            cond = (vars[ap.lval.var] != vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_LT_VAR:
            cond = (vars[ap.lval.var] < vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_GT_VAR:
            cond = (vars[ap.lval.var] > vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_LE_VAR:
            cond = (vars[ap.lval.var] <= vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_OP_GE_VAR:
            cond = (vars[ap.lval.var] >= vars[ap.rval.var]);
            break;
          case kripke_cube::VAR_DEAD:
            break;
          default:
            SPOT_ASSERT(false);
          }

        if (cond)
          cubeset_.set_true_var(c, i);
        else
          cubeset_.set_false_var(c, i);
      }
  }

  // FIXME I think we only need visbles aps, i.e. if the system has
  // following variables, i.e. P_0.var1 and P_0.var2 but the property
  // automaton only mention P_0.var2, we do not need to capture (in
  // the resulting cube) any atomic proposition for P_0.var1
  void
  kripkecube<cspins_state,
             cspins_iterator>::match_aps(const atomic_prop_set* aps,
                                         std::string dead_prop)
  {
    // Keep trace of errors
    int errors = 0;
    std::ostringstream err;

    // First we capture state name of each processes.
    int type_count = d_->get_type_count();
    typedef std::map<std::string, int> enum_map_t;
    std::vector<enum_map_t> enum_map(type_count);
    std::unordered_map<std::string, int> matcher;
    for (int i = 0; i < type_count; ++i)
      {
        matcher[d_->get_type_name(i)] = i;
        int enum_count = d_->get_type_value_count(i);
        for (int j = 0; j < enum_count; ++j)
          enum_map[i].emplace(d_->get_type_value_name(i, j), j);
      }

    // Then we extract the basic atomics propositions from the Kripke
    std::vector<std::string> k_aps;
    int state_size = d_->get_state_size();
    for (int i = 0; i < state_size; ++i)
      k_aps.push_back(d_->get_state_variable_name(i));

    int i  = -1;
    for (auto ap_f: *aps)
      {
        ++i;

        bool positive = true;
        std::string ap;
        if (ap_f.is(op::Not))
          {
            positive = false;
            ap = ap_f[0].ap_name();
          }
        else
          ap = ap_f.ap_name();

        aps_str_.push_back(ap);

        // Grab dead property
        if (ap.compare(dead_prop) == 0)
          {
            dead_idx_ = i;
            pset_.push_back({i, kripke_cube::VAR_DEAD, 0, false}); // FIXME false
            continue;
          }

        // Get ap name and remove all extra whitespace
        ap.erase(std::remove_if(ap.begin(), ap.end(),
                                [](char x){
                                  return std::isspace(x);
                                }),
                 ap.end());

        // Look if it is a well known atomic proposition
        auto it = std::find(k_aps.begin(), k_aps.end(), ap);
        if (it != k_aps.end())
          {
            // The aps is directly an AP of the system, we will just
            // have to detect if the variable is 0 or not.
            pset_.push_back({(int)std::distance(k_aps.begin(), it),
                  kripke_cube::VAR_OP_NE, 0, positive});
            continue;
          }

        // The ap is not known. We distinguish many cases:
        //     - It is a State name, i.e P_0.S or P_0 == S
        //     - It refers a specific variable value, i.e. P_0.var == 2,
        //       P_0.var < 2, P_0.var != 2, ...
        //     - It's an unknown variable
        // Note that we do not support P_0.state1 == 12 since we do not
        // know how to interpret such atomic proposition.

        // We split the formula according to operators
        std::size_t found_op_first = ap.find_first_of("=<>!");
        std::size_t found_op_last = ap.find_last_of("=<>!");
        std::string left;
        std::string right;
        std::string ap_error;
        std::string op;

        if (found_op_first == 0 || found_op_last == ap.size()-1)
          {
            err << "Invalid operator use in " << ap << '\n';
            ++errors;
            continue;
          }

        if (std::string::npos == found_op_first)
          {
            left = ap;
            right = "";
            op = "";
          }
        else
          {
            left = ap.substr(0, found_op_first);
            right = ap.substr(found_op_last+1, ap.size()-found_op_last);
            op  = ap.substr(found_op_first, found_op_last+1-found_op_first);
          }

        // Variables to store the left part of the atomic proposition
        bool left_is_digit = false;
        int lval;

        // Variables to store the right part of the atomic proposition
        bool right_is_digit = false;
        int rval;

        // And finally the operator
        kripke_cube::relop oper;


        // Now, left and (possibly) right may refer atomic
        // propositions or specific state inside of a process.
        // First check if it is a known atomic proposition
        it = std::find(k_aps.begin(), k_aps.end(), left);
        if (it != k_aps.end())
          {
            // The ap is directly an AP of the system, we will just
            // have to detect if the variable is 0 or not.
            lval = std::distance(k_aps.begin(), it);
          }
        else
          {
            // Detect if it is a process state
            std::size_t found_dot = left.find_first_of('.');
            if (std::string::npos != found_dot)
              {
                std::string proc_name = left.substr(0, found_dot);
                std::string st_name = left.substr(found_dot+1,
                                                  left.size()-found_dot);

                auto ni = matcher.find(proc_name);
                if (ni == matcher.end())
                  {
                    ap_error = left;
                    goto error_ap_unknown;
                  }
                int type_num = ni->second;
                enum_map_t::const_iterator ei =
                  enum_map[type_num].find(st_name);
                if (ei == enum_map[type_num].end())
                  {
                    ap_error = left;
                    goto error_ap_unknown;
                  }

                if (right.compare("") != 0)
                  {
                    // We are in the case P.state1 == something.. We don't
                    // know how to interpret this.
                    ap_error = op + right;
                    err << "\nOperation " << op  << " in \"" << ap_error
                        << "\" is not available for process's state"
                        << " (i.e. " <<  left << ")\n";
                    ++errors;
                    continue;
                  }

                pset_.push_back(
                  { (int) std::distance(k_aps.begin(), std::find(k_aps.begin(),
                  k_aps.end(), proc_name)), kripke_cube::VAR_OP_EQ, ei->second,
                  positive });
                continue;
              }
            else
              {
                // Finally, it's a number...
                left_is_digit = true;
                for (auto c: left)
                  if (!isdigit(c))
                    left_is_digit = false;

                if (left_is_digit)
                  lval = std::strtol (left.c_str(), nullptr, 10);
                else
                  {
                    // ... or something like: State1 == P_0
                    // so it doesn't contains '.'
                    if (std::string::npos != right.find_first_of('.'))
                      {
                        err << "\nOperation \"" << right
                            << "\" does not refer a process"
                            << " (i.e. " << left << " is not valid)\n";
                        ++errors;
                        continue;
                      }

                    // or something like: P_0 == State1
                    auto ni = matcher.find(right);
                    if (ni == matcher.end())
                      {
                        ap_error = ap;
                        goto error_ap_unknown;
                      }
                    int type_num = ni->second;
                    enum_map_t::const_iterator ei =
                      enum_map[type_num].find(left);
                    if (ei == enum_map[type_num].end())
                      {
                        ap_error = left;
                        goto error_ap_unknown;
                      }
                    pset_.push_back(
                      { (int) std::distance(k_aps.begin(),
                      std::find(k_aps.begin(), k_aps.end(), right)),
                      kripke_cube::VAR_OP_EQ, ei->second, positive });
                    continue;
                  }
              }
          }

        // Here Left is known. Just detect cases where left is digit there is
        // no right part.
        if (left_is_digit && right.empty())
          {
            ap_error = ap;
            goto error_ap_unknown;
          }

        SPOT_ASSERT(!right.empty() && !op.empty());

        // Parse right part of the atomic proposition
        // Check if it is a known atomic proposition
        it = std::find(k_aps.begin(), k_aps.end(), right);
        if (it != k_aps.end())
          {
            // The aps is directly an AP of the system, we will just
            // have to detect if the variable is 0 or not.
            rval = std::distance(k_aps.begin(), it);
          }
        else
          {
            // We are is the right part, so  if it is a process state
            // we do not know how to interpret (xxx == P.state1). Abort
            std::size_t found_dot = right.find_first_of('.');
            if (std::string::npos != found_dot)
              {
                ap_error = left + op;
                err << "\nOperation " << op  << " in \"" << ap_error
                    << "\" is not available for process's state"
                    << " (i.e. " <<  right << ")\n";
                ++errors;
                continue;
              }
            else
              {
                // Finally, it's a number
                right_is_digit = true;
                for (auto c: right)
                  if (!isdigit(c))
                    right_is_digit = false;

                if (right_is_digit)
                  rval = std::strtol (right.c_str(), nullptr, 10);
                else
                  {
                    if (std::string::npos != left.find_first_of('.'))
                      {
                        err << "\nProposition \"" << ap
                            << "\" cannot be interpreted"
                            << " (i.e. " << op + right << " is not valid)\n";
                        ++errors;
                        continue;
                      }

                    // or something like: P_0 == State1
                    auto ni = matcher.find(left);
                    if (ni == matcher.end())
                      {

                        ap_error = left;
                        goto error_ap_unknown;
                      }
                    int type_num = ni->second;
                    enum_map_t::const_iterator ei =
                      enum_map[type_num].find(right);
                    if (ei == enum_map[type_num].end())
                      {
                        ap_error = right;
                        goto error_ap_unknown;
                      }
                    pset_.push_back({
                        (int) std::distance(k_aps.begin(),
                                            std::find(k_aps.begin(),
                                                      k_aps.end(), left)),
                          kripke_cube::VAR_OP_EQ,  ei->second, positive});
                    continue;
                  }
              }
          }

        if (left_is_digit && right_is_digit)
          {
            err << "\nOperation \"" << op
                << "\" between two numbers not available"
                << " (i.e. " << right << " and, "
                << left << ")\n";
            ++errors;
            continue;
          }

        // Left and Right are know, just analyse the operator.
        if (op.compare("==") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_EQ_VAR :
            (left_is_digit? kripke_cube::OP_EQ_VAR : kripke_cube::VAR_OP_EQ);
        else if (op.compare("!=") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_NE_VAR :
            (left_is_digit? kripke_cube::OP_NE_VAR : kripke_cube::VAR_OP_NE);
        else if (op.compare("<") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_LT_VAR :
            (left_is_digit? kripke_cube::OP_LT_VAR : kripke_cube::VAR_OP_LT);
        else if (op.compare(">") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_GT_VAR :
            (left_is_digit? kripke_cube::OP_GT_VAR : kripke_cube::VAR_OP_GT);
        else if (op.compare("<=") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_LE_VAR :
            (left_is_digit? kripke_cube::OP_LE_VAR : kripke_cube::VAR_OP_LE);
        else if (op.compare(">=") == 0)
          oper = !left_is_digit && !right_is_digit? kripke_cube::VAR_OP_GE_VAR :
            (left_is_digit? kripke_cube::OP_GE_VAR : kripke_cube::VAR_OP_GE);
        else
          {
            err << "\nOperation \"" << op
                << "\" is unknown\n";
            ++errors;
            continue;
          }

        pset_.push_back({lval, oper, rval, positive});
        continue;

      error_ap_unknown:
        err << "\nProposition \"" << ap_error << "\" does not exist\n";
        ++errors;
        continue;
      }

    if (errors)
      throw std::runtime_error(err.str());
  }
}
