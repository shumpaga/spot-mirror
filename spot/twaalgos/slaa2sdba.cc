// -*- coding: utf-8 -*-
// Copyright (C) 2018 Laboratoire de Recherche et Développement de
// l'Epita (LRDE).
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

#include "config.h"
#include <spot/twaalgos/slaa2sdba.hh>
#include <spot/twa/twagraph.hh>
#include <spot/misc/minato.hh>
#include <spot/misc/bitvect.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/totgba.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/isweakscc.hh>
#include <sstream>
#include <tuple>

namespace spot
{
  enum slaa_to_sdba_state_type { State_Leave, State_MayStay, State_MustStay };

  static std::vector<slaa_to_sdba_state_type>
  slaa_to_sdba_state_types(const_twa_graph_ptr aut)
  {
    if (!(aut->acc().is_co_buchi() || aut->acc().is_t()))
      throw std::runtime_error
        ("slaa_to_sdba only works for co-Büchi acceptance");

    unsigned n = aut->num_states();

    std::vector<slaa_to_sdba_state_type> res;
    res.reserve(n);

    for (unsigned s = 0; s < n; ++s)
      {
        bool must_stay = true;
        bool may_stay = false;
        for (auto& e : aut->out(s))
          {
            bool have_self_loop = false;
            for (unsigned d: aut->univ_dests(e))
              if (d == s)
                {
                  have_self_loop = true;
                  if (!e.acc)
                    may_stay = true;
                }
            if (!have_self_loop)
              must_stay = false;
          }
        res.push_back(must_stay ? State_MustStay :
                      may_stay ? State_MayStay :
                      State_Leave);
      }

    return res;
  }

  twa_graph_ptr
  slaa_to_sdba_highlight_state_types(twa_graph_ptr aut)
  {
    auto v = slaa_to_sdba_state_types(aut);

    auto hs =
      aut->get_named_prop<std::map<unsigned, unsigned>>("highlight-states");
    if (!hs)
      {
        hs = new std::map<unsigned, unsigned>;
        aut->set_named_prop("highlight-states", hs);
      }

    unsigned n = v.size();
    for (unsigned i = 0; i < n; ++i)
      {
        slaa_to_sdba_state_type t = v[i];
        (*hs)[i] = (t == State_Leave) ? 4 : (t == State_MayStay) ? 2 : 5;
      }

    return aut;
  }


  namespace
  {
    typedef std::tuple<bdd, bdd, bdd> triplet_t;

    struct triplet_lt final
    {
      bool operator()(const triplet_t& left, const triplet_t& right) const
      {
        bdd l_left;
        bdd l_right;
        bdd l_component;
        std::tie(l_left, l_right, l_component) = left;
        bdd r_left;
        bdd r_right;
        bdd r_component;
        std::tie(r_left, r_right, r_component) = right;

        if (l_left.id() < r_left.id())
          return true;
        if (l_left.id() > r_left.id())
          return false;
        if (l_right.id() < r_right.id())
          return true;
        if (l_right.id() > r_right.id())
          return false;
        return l_component.id() < r_component.id();
      }
    };

    typedef std::pair<unsigned, bdd> pair_state_component_t;
    struct pair_state_component_lt final
    {
      bool operator()(const pair_state_component_t& left,
                      const pair_state_component_t& right) const
      {
        if (left.first < right.first)
          return true;
        if (left.first > right.first)
          return false;
        return left.second.id () < right.second.id();
      }
    };

    static bdd positive_bdd(bdd in)
    {
      assert(in != bddfalse);
      bdd res = bddtrue;
      while (in != bddtrue)
        {
          bdd low = bdd_low(in);
          if (low == bddfalse)
            {
              res &= bdd_ithvar(bdd_var(in));
              in = bdd_high(in);
            }
          else
            {
              in = low;
            }
        }
      return res;
    }

    static bool is_positive_conjunct(bdd in)
    {
      bdd res = bddtrue;
      while (in != bddtrue)
        {
          if (bdd_low(in) == bddfalse)
            in = bdd_high(in);
          else
            return false;
        }
      return true;
    }

    class slaa_to_sdba_runner final
    {
      const_twa_graph_ptr aut_;
      // We use n BDD variables to represent the n states.
      // The state i is BDD variable state_base_+i.
      unsigned state_base_;
      std::vector<unsigned> dots_;
      std::map<int, acc_cond::mark_t> var_to_mark_;
      std::vector<bdd> mark_to_state_as_bdd_;
      bdd all_states_and_dots_;
      std::map<triplet_t, unsigned, triplet_lt> triplet_to_unsigned_;
      std::vector<triplet_t> unsigned_to_triplet_;
      std::vector<slaa_to_sdba_state_type> types_;
      std::vector<bool> true_state_;
      bitvect_array* reachable_;
      bitvect* det_from_now_;
      bdd stay_states;
      bdd neg_must_stay_states;
    public:
      slaa_to_sdba_runner(const_twa_graph_ptr aut)
        : aut_(aut)
      {
        types_ = slaa_to_sdba_state_types(aut);

        // Scan for true states.  Since we work with co-Büchi SLAA,
        // any state with a non-rejecting true self-loop is a true
        // state.
        unsigned ns = aut->num_states();
        true_state_.resize(ns, false);
        for (unsigned s = 0; s < ns; ++s)
          for (auto& e: aut_->out(s))
            if (e.dst == s && e.cond == bddtrue && !e.acc)
              {
                true_state_[s] = true;
                break;
              }

        reachable_ = make_bitvect_array(ns, ns);
        det_from_now_ = make_bitvect(ns);
        scc_info si(aut, scc_info_options::TRACK_SUCCS
                    | scc_info_options::TRACK_STATES);
        unsigned nscc = si.scc_count();
        for (unsigned n = 0; n < nscc; ++n)
          {
            unsigned s = si.one_state_of(n);
            bitvect& svec = reachable_->at(s);
            svec.set(s);
            bool isdet = true;
            for (unsigned d: si.succ(n))
              {
                unsigned dd = si.one_state_of(d);
                svec |= reachable_->at(dd);
                isdet &= det_from_now_->get(dd);
              }
            if (isdet && is_weak_scc(si, n))
              {
                bdd available = bddtrue;
                for (auto&t : aut->out(s))
                  if (!bdd_implies(t.cond, available))
                    {
                      isdet = false;
                      break;
                    }
                  else
                    {
                      available -= t.cond;
                    }
                if (isdet)
                  det_from_now_->set(s);
              }
          }
      }

      ~slaa_to_sdba_runner()
      {
        aut_->get_dict()->unregister_all_my_variables(this);
        delete reachable_;
        delete det_from_now_;
      }

      bdd state_bdd(unsigned num)
      {
        return bdd_ithvar(state_base_ + num);
      }

      bdd dot_bdd(unsigned num)
      {
        return bdd_ithvar(dots_[num]);
      }

      unsigned state_num(bdd b)
      {
        unsigned num = bdd_var(b) - state_base_;
        assert(num < aut_->num_states());
        return num;
      }


      std::map<pair_state_component_t, bdd,
               pair_state_component_lt> state_comp_as_bdd_cache_;
      bdd state_comp_as_bdd(unsigned num, bdd component)
      {
        auto p = state_comp_as_bdd_cache_.emplace
          (pair_state_component_t{num, component}, bddfalse);
        if (p.second)
          {
            bdd res = bddfalse;
            if (true_state_[num])
              res = bddtrue;
            else
              {
                bool stay_here = bdd_implies(component, state_bdd(num));

                for (auto& e: aut_->out(num))
                  {
                    bdd dest = bddtrue;
                    bool seen_self = false;
                    for (unsigned d: aut_->univ_dests(e.dst))
                      {
                        if (d == num)
                          seen_self = true;
                        if (!true_state_[d])
                          dest &= state_bdd(d);
                        if (d == e.src && e.acc && component == bddtrue)
                          dest &= dot_bdd(d);
                      }
                    if (!stay_here || (seen_self && !e.acc))
                      res |= e.cond & dest;
                  }
                if (component != bddtrue)
                  // We cannot enter must stay_states that do not belong
                  // to this component.
                  res = bdd_restrict(res,
                                     bdd_exist(neg_must_stay_states,
                                               component));
              }
            p.first->second = res;
          }
        return p.first->second;
      }

      std::string print_state_formula(bdd f)
      {
        if (f == bddtrue)
          return "⊤";
        if (f == bddfalse)
          return "⊥";
        std::ostringstream out;
        minato_isop isop(f);
        bdd cube;
        bool outer_first = true;
        unsigned ns = aut_->num_states();
        while ((cube = isop.next()) != bddfalse)
          {
            if (!outer_first)
              out << '+';
            bool inner_first = true;
            while (cube != bddtrue)
              {
                assert(bdd_low(cube) == bddfalse);
                unsigned v = bdd_var(cube);
                if ((state_base_ <= v) && (v < state_base_ + ns))
                  {
                    if (!inner_first)
                      out << "·";
                    out << state_num(cube);
                  }
                cube = bdd_high(cube);
                inner_first = false;
              }
            outer_first = false;
          }
        return out.str();
      }

      std::string print_pair(const triplet_t& p)
      {
        std::string res =  print_state_formula(std::get<0>(p));
        bdd comp = std::get<2>(p);
        if (comp != bddtrue)
          res += ("\n" + print_state_formula(std::get<1>(p))
                  + "\n" + print_state_formula(std::get<2>(p)));
        return res;
      }


      twa_graph_ptr run(bool three_plus, bool force_jump, bool fb_opt,
                        bool no_dealternation)
      {
        auto d = aut_->get_dict();

        // We need one BDD variable per possible output state.  If
        unsigned ns = aut_->num_states();
        state_base_ = d->register_anonymous_variables(ns, this);
        dots_.resize(ns);
        all_states_and_dots_ = bddtrue;
        bdd all_dots = bddtrue;
        bdd leave_states = bddtrue;
        stay_states = bddtrue;
        bdd may_stay_states = bddtrue;
        bdd must_stay_states = bddtrue;
        neg_must_stay_states = bddtrue;
        bitvect* stay_bitvect = make_bitvect(ns);
        bitvect* tmp_bitvect = make_bitvect(ns);
        unsigned mark_pos = 0;
        acc_cond::mark_t all_marks = {};
        std::vector<bool> has_rejecting_loop;
        has_rejecting_loop.reserve(ns);

        for (unsigned s = 0; s < ns; ++s)
          {
            {
              bool rl = false;
              for (auto& e: aut_->out(s))
                if (!aut_->acc().accepting(e.acc))
                  {
                    rl = true;
                    break;
                  }
              has_rejecting_loop.push_back(rl);
            }
            bdd sb = state_bdd(s);
            all_states_and_dots_ &= sb;
            if (types_[s] == State_Leave)
              {
                leave_states &= sb;
              }
            else
              {
                stay_states &= sb;
                stay_bitvect->set(s);
                if (types_[s] == State_MayStay)
                  {
                    may_stay_states &= sb;
                  }
                else
                  {
                    must_stay_states &= sb;
                    neg_must_stay_states &= !sb;
                  }
              }

            // If this state has a rejecting loop, we will need some
            // mark.
            for (auto& e: aut_->out(s))
              {
                if (!e.acc)
                  continue;
                for (unsigned dst: aut_->univ_dests(e))
                  if (dst == s)
                    {
                      int v = d->register_anonymous_variables(1, this);
                      dots_[s] = v;
                      all_marks.set(mark_pos);
                      var_to_mark_.emplace(v, acc_cond::mark_t({mark_pos++}));
                      all_dots &= dot_bdd(s);
                      mark_to_state_as_bdd_.push_back(state_bdd(s));
                      goto needmark_done;
                    }
              }
          needmark_done:;
          }
        all_states_and_dots_ &= all_dots;
        // Make sure we use at least one acceptance set, otherwise we
        // cannot distinguish the non-accepting part from the
        // accepting one.
        if (mark_pos == 0)
          all_marks.set(mark_pos++);

        twa_graph_ptr res = make_twa_graph(d);
        res->copy_ap_of(aut_);
        res->prop_copy(aut_, {false, false, false, false, false, true});
        res->set_generalized_buchi(mark_pos);

        // Compute a BDD representing the successors of one
        // state-formula.  The successor representation is a
        // BDD formula containing letters and state numbers.
        auto delta = [&](bdd left, bdd component)
          {
            bdd res = bddfalse;
            minato_isop isop(left);
            bdd cube;
            while ((cube = isop.next()) != bddfalse)
              {
                bdd inner_res = bddtrue;
                while (cube != bddtrue)
                  {
                    assert(bdd_low(cube) == bddfalse);
                    unsigned v = bdd_var(cube);
                    if ((state_base_ <= v) && (v < state_base_ + ns))
                      {
                        unsigned num = state_num(cube);
                        inner_res &= state_comp_as_bdd(num, component);
                      }
                    cube = bdd_high(cube);
                  }
                res |= inner_res;
              }
            return res;
          };

        std::vector<unsigned> todo;

        auto new_state = [&](triplet_t s)
          {
            auto p = triplet_to_unsigned_.emplace(s, 0U);
            if (p.second)
              {
                unsigned_to_triplet_.push_back(s);
                unsigned q = res->new_state();
                assert(unsigned_to_triplet_.size() == q + 1);
                todo.push_back(q);
                p.first->second = q;
              }
            return p.first->second;
          };

        auto determinizable = [&](bdd config)
          {
            if (no_dealternation)
              return false;
            bool alldet = true;
            while (config != bddtrue)
              {
                unsigned num = state_num(config);
                if (!det_from_now_->get(num))
                  {
                    alldet = false;
                    break;
                  }
                config = bdd_high(config);
              }
            return alldet;
          };

        auto may_not_force_jump = [&](bdd config)
          {
            tmp_bitvect->clear_all();
            // Get all the states reachable from config.
            bdd oldconfig = config;
            while (config != bddtrue)
              {
                unsigned num = state_num(config);
                // Current config should not have any orange.
                if (types_[num] == State_MayStay)
                  return true;
                if (types_[num] == State_MustStay && has_rejecting_loop[num])
                  return true;
                *tmp_bitvect |= reachable_->at(num);
                config = bdd_high(config);
              }
            // We are only after new red states, so erase the old
            // ones.
            while (oldconfig != bddtrue)
              {
                unsigned num = state_num(oldconfig);
                tmp_bitvect->clear(num);
                oldconfig = bdd_high(oldconfig);
              }
            // Did we introduce some new red or orange state?
            for (unsigned n = 0; n < ns; ++n)
              if (types_[n] != State_Leave && tmp_bitvect->get(n)
                  && !true_state_[n])
                return true;
            return false;
          };

        auto stay_states_can_reach_all = [&](bdd staystates, bdd dst)
          {
            tmp_bitvect->clear_all();
            while (staystates != bddtrue)
              {
                unsigned num = state_num(staystates);
                staystates = bdd_high(staystates);
                *tmp_bitvect |= reachable_->at(num);
              }
            while (dst != bddtrue)
              {
                unsigned num = state_num(dst);
                dst = bdd_high(dst);
                if (!tmp_bitvect->get(num))
                  return false;
              }
            return true;
          };

        unsigned init_state = aut_->get_init_state_number();
        bdd init = bddtrue;
        for (unsigned d: aut_->univ_dests(init_state))
          init &= state_bdd(d);
        bdd init_right = bddfalse;
        bdd init_comp = bddtrue;

        // Do we want to enter the breakpoint construction initially?
        // We only do that if
        // (1) the standard dealternation cannot be applied to
        // get a deterministic automaton.
        // (2) we have some red states,
        // (3) we do not have any orange state, cannot
        // reach new orange or red states, and the current
        // red states have no rejecting self-loops.
        // (4) those red states can reach all states in the
        // initial configuration (assuming optimization three_plus
        // is turned on).
        if (!determinizable(init)) // (1)
          {
            bdd local_stay_states = bdd_exist(init, leave_states);
            bdd local_must_stay_states =
              bdd_exist(local_stay_states, may_stay_states);
            if (local_must_stay_states != bddtrue) // 2
              {
                if (!may_not_force_jump(init) // (3)
                    && (!three_plus
                        || stay_states_can_reach_all(local_stay_states,
                                                     init))) // (4)
                  {
                    init_comp = init_right = local_stay_states;
                    init = bdd_restrict(init, init_right);
                  }
              }
          }

        new_state(triplet_t{init, init_right, init_comp});

        while (!todo.empty())
          {
            unsigned src_state = todo.back();
            triplet_t src = unsigned_to_triplet_[src_state];
            todo.pop_back();
            bdd src_left, src_right, src_comp;
            std::tie(src_left, src_right, src_comp) = src;

            if (src_comp == bddtrue)
              {
                if (src_left == bddtrue)
                  {
                    res->new_edge(src_state, src_state, bddtrue, all_marks);
                    continue;
                  }
                 // If all the states of src_left can only reach
                 // deterministic (and weak) states, we want to do the
                 // standard dealternation construction.
                else if (determinizable(src_left))
                  {
                    bdd dests = delta(src_left, bddtrue);
                    bdd aps = bdd_exist(bdd_support(dests),
                                        all_states_and_dots_);
                    bdd all = bddtrue;
                    while (all != bddfalse)
                      {
                        bdd letter = bdd_satoneset(all, aps, bddfalse);
                        all -= letter;
                        bdd trans = bdd_restrict(dests, letter);
                        if (trans == bddfalse)
                          continue;
                        bdd dest = bdd_exist(trans, all_dots);
                        bdd dots = bdd_existcomp(trans, all_dots);
                        acc_cond::mark_t m = all_marks;
                        while (dots != bddtrue)
                          {
                            int v = bdd_var(dots);
                            m -= var_to_mark_[v];
                            dots = bdd_high(dots);
                          }
                        if (!mark_to_state_as_bdd_.empty())
                          for (unsigned m1: m.sets())
                            {
                              bdd m1state = mark_to_state_as_bdd_[m1];
                              if (bdd_implies(dest, m1state)
                                  && !bdd_implies(src_left, m1state))
                                m.clear(m1);
                            }
                        res->new_edge(src_state,
                                      new_state(triplet_t{dest,
                                            bddfalse, bddtrue}),
                                      letter, m);

                        }
                    continue;
                  }
              }

            bdd succs_left = bdd_restrict(delta(src_left, src_comp), src_comp);
            bdd succs_right = delta(src_right, src_comp);
            // decompose successors
            bdd aps = bdd_exist(bdd_support(succs_left), all_states_and_dots_) &
              bdd_exist(bdd_support(succs_right), all_states_and_dots_);
            bdd all = bddtrue;
            while (all != bddfalse)
              {
                bdd letter = bdd_satoneset(all, aps, bddfalse);
                all -= letter;
                bdd left = bdd_restrict(succs_left, letter);
                assert(bdd_exist(left, aut_->ap_vars()) == left);
                if (left == bddfalse)
                  continue;
                if (src_comp != bddtrue)
                  {
                    bdd right = bdd_restrict(succs_right, letter);
                    assert(bdd_exist(right, aut_->ap_vars()) == right);
                    if (right == bddfalse && src_comp != bddtrue)
                      continue;
                    acc_cond::mark_t m = {};
                    unsigned dst;
                    if (left == bddtrue)
                      {
                        bdd dst_left = bdd_restrict(right, src_comp);
                        bdd dst_right = src_comp;
                        if (is_positive_conjunct(dst_left))
                          dst_right = bdd_restrict(dst_right, dst_left);
                        dst = new_state(triplet_t{dst_left, dst_right,
                                                  src_comp});
                        m = all_marks;
                        if (fb_opt && src_left == bddtrue && src_state != dst)
                          m = {};
                      }
                    else
                      {
                        if (is_positive_conjunct(left))
                          right = bdd_restrict(right, left);
                        dst = new_state(triplet_t
                                        {left, right, src_comp});
                      }
                    res->new_edge(src_state, dst, letter, m);
                  }
                else
                  {
                    minato_isop isop(left);
                    bdd cube;
                    while ((cube = isop.next()) != bddfalse)
                      {
                        bdd dst = bdd_exist(cube, all_dots);

                        bool did_jump = false;
                        if (bdd_restrict(dst, stay_states) != dst
                            && !determinizable(dst))
                          {
                            bdd local_stay_states =
                              bdd_exist(dst, leave_states);
                            do
                              {
                                if (local_stay_states == bddfalse)
                                  break;
                                bdd local_may_stay_states =
                                  bdd_exist(local_stay_states,
                                            must_stay_states);
                                bdd local_must_stay_states =
                                  bdd_exist(local_stay_states, may_stay_states);

                                bdd all = bddtrue;
                                while (all != bddfalse)
                                  {
                                    bdd c = bdd_satoneset(all,
                                                          local_may_stay_states,
                                                          bddfalse);
                                    all -= c;
                                    c = positive_bdd(c);
                                    c &= local_must_stay_states;
                                    if (c == bddtrue)
                                      continue;

                                    if (three_plus
                                        && !stay_states_can_reach_all(c, dst))
                                      continue;

                                    bdd succ = bdd_restrict(dst, c);
                                    res->new_edge(src_state,
                                                  new_state(triplet_t{succ,
                                                                      c, c}),
                                                  letter);

                                    did_jump = true;
                                  }
                              }
                            while (false);
                          }
                        // If there is new red or orange reachable
                        // from dst, and if there is no orange in dst,
                        // then there is no point in staying in the
                        // non-det part.  We call this the force_jump
                        // optimization.
                        if (!force_jump || !did_jump
                            || may_not_force_jump(dst))
                          res->new_edge(src_state,
                                        new_state(triplet_t{dst,
                                              bddfalse, bddtrue}),
                                        letter);
                      }
                  }
              }
            }

        delete stay_bitvect;
        delete tmp_bitvect;
        res->merge_edges();

        auto* names = new std::vector<std::string>;
        res->set_named_prop("state-names", names);
        unsigned output_size = unsigned_to_triplet_.size();
        names->reserve(output_size);
        for (unsigned s = 0; s < output_size; ++s)
          names->emplace_back(print_pair(unsigned_to_triplet_[s]));

        return res;
      }


    };
  }

  twa_graph_ptr
  slaa_to_sdba(const_twa_graph_ptr aut, bool force_build,
               bool three_plus, bool force_jump, bool fb_opt,
               bool no_dealternation)
  {
    if (!force_build && is_deterministic(aut))
      return to_generalized_buchi(aut);

    slaa_to_sdba_runner runner(aut);
    return runner.run(three_plus, force_jump, fb_opt, no_dealternation);
  }

}
