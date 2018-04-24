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
#include <sstream>
#include <tuple>

namespace spot
{
  static std::vector<slaa_to_sdba_state_type>
  slaa_to_sdba_state_types(const_twa_graph_ptr aut)
  {
    if (!aut->acc().is_co_buchi())
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

  SPOT_API twa_graph_ptr
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
      bool operator()(const pair_state_component_t& left, const pair_state_component_t& right) const
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

    class slaa_to_sdba_runner final
    {
      const_twa_graph_ptr aut_;
      // We use n BDD variables to represent the n states.
      // The state i is BDD variable state_base_+i.
      unsigned state_base_;
      unsigned dot_base_;
      bdd all_states_and_dots_;
      bool cutdet_;
      std::map<triplet_t, unsigned, triplet_lt> triplet_to_unsigned_;
      std::vector<triplet_t> unsigned_to_triplet_;
      std::vector<slaa_to_sdba_state_type> types_;
      std::vector<bool> true_state_;
      bitvect_array* reachable_;
      bdd stay_states;
    public:
      slaa_to_sdba_runner(const_twa_graph_ptr aut, bool cutdet)
        : aut_(aut), cutdet_(cutdet)
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
        scc_info si(aut, scc_info_options::TRACK_SUCCS);
        unsigned nscc = si.scc_count();
        for (unsigned n = 0; n < nscc; ++n)
          {
            unsigned s = si.one_state_of(n);
            bitvect& svec = reachable_->at(s);
            svec.set(s);
            for (unsigned d: si.succ(n))
              svec |= reachable_->at(d);
          }
      }

      ~slaa_to_sdba_runner()
      {
        aut_->get_dict()->unregister_all_my_variables(this);
        delete reachable_;
      }

      bdd state_bdd(unsigned num)
      {
        return bdd_ithvar(state_base_ + num);
      }

      bdd dot_bdd(unsigned num)
      {
        return bdd_ithvar(dot_base_ + num);
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
              out << "+";
            bool inner_first = true;
            while (cube != bddtrue)
              {
                assert(bdd_low(cube) == bddfalse);
                unsigned v = bdd_var(cube);
                if (!inner_first)
                  out << "·";
                if ((dot_base_ <= v) && (v < dot_base_ + ns))
                  out << '(' << (v - dot_base_) << ')';
                else
                  out << state_num(cube);
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


      twa_graph_ptr run()
      {
        auto d = aut_->get_dict();

        twa_graph_ptr res = make_twa_graph(d);
        res->copy_ap_of(aut_);
        res->prop_copy(aut_, {false, false, false, false, false, true});
        res->set_buchi();

        // We need one BDD variable per possible output state.  If
        unsigned ns = aut_->num_states();
        state_base_ = d->register_anonymous_variables(ns, this);
        dot_base_ = d->register_anonymous_variables(ns, this);
        all_states_and_dots_ = bddtrue;
        bdd all_dots = bddtrue;
        bdd leave_states = bddtrue;
        stay_states = bddtrue;
        bdd may_stay_states = bddtrue;
        bdd must_stay_states = bddtrue;
        bitvect* stay_bitvect = make_bitvect(ns);
        bitvect* tmp_bitvect = make_bitvect(ns);

        for (unsigned s = 0; s < ns; ++s)
          {
            bdd sb = state_bdd(s);
            bdd db = dot_bdd(s);
            all_states_and_dots_ &= sb & db;
            all_dots &= db;
            if (types_[s] == State_Leave)
              {
                leave_states &= sb;
              }
            else
              {
                stay_states &= sb;
                stay_bitvect->set(s);
                if (types_[s] == State_MayStay)
                  may_stay_states &= sb;
                else
                  must_stay_states &= sb;
              }
          }
        if (cutdet_)
          {
            leave_states &= all_dots;
            must_stay_states &= all_dots;
            may_stay_states &= all_dots;
          }

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
                    if ((dot_base_ <= v) && (v < dot_base_ + ns))
                      {
                        cube = bdd_high(cube);
                        continue;
                      }
                    unsigned num = state_num(cube);
                    if (num < ns)
                      inner_res &= state_comp_as_bdd(num, component);
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

        bdd init = state_bdd(aut_->get_init_state_number());
        new_state({init, bddfalse, bddtrue});

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
                    res->new_edge(src_state, src_state, bddtrue, {0});
                    continue;
                  }
                else if (bdd_restrict(src_left, stay_states) != src_left)
                  {
                    minato_isop isop(src_left);
                    bdd cube;
                    while ((cube = isop.next()) != bddfalse)
                      {
                        bdd local_stay_states = bdd_exist(cube, leave_states);
                        if (local_stay_states == bddfalse)
                          continue;
                        bdd local_may_stay_states =
                          bdd_exist(local_stay_states, must_stay_states);
                        bdd local_must_stay_states =
                          bdd_exist(local_stay_states, may_stay_states);

                        bdd all = bddtrue;
                        while (all != bddfalse)
                          {
                            bdd c = bdd_satoneset(all, local_may_stay_states,
                                                  bddfalse);
                            all -= c;
                            c = positive_bdd(c);
                            c &= local_must_stay_states;
                            if (c == bddtrue)
                              continue;
                            bdd succs = delta(cube, c);
                            bdd aps =
                              bdd_exist(bdd_support(succs), all_states_and_dots_);
                            succs = bdd_restrict(succs, c);
                            bdd all_aps = bddtrue;
                            while (all_aps != bddfalse)
                              {
                                bdd letter = bdd_satoneset(all_aps, aps,
                                                           bddfalse);
                                all_aps -= letter;
                                bdd dest = bdd_restrict(succs, letter);
                                dest = bdd_exist(dest, all_dots);
                                if (dest != bddfalse)
                                  res->new_edge(src_state,
                                                new_state({dest, c, c}),
                                                letter);
                              }
                          }
                      }
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
                if (left == bddfalse)
                  continue;
                if (cutdet_ || src_comp != bddtrue)
                  {
                    bdd right = bdd_restrict(succs_right, letter);
                    if (right == bddfalse && src_comp != bddtrue)
                      continue;
                    acc_cond::mark_t m = {};
                    unsigned dst;
                    if (left == bddtrue && src_comp != bddtrue)
                      {
                        dst = new_state({bdd_restrict(right, src_comp), src_comp, src_comp});
                        m = {0};
                      }
                    else
                      {
                        dst = new_state({left, right, src_comp});
                      }
                    res->new_edge(src_state, dst, letter, m);
                  }
                else
                  {
                    minato_isop isop(left);
                    bdd cube;
                    while ((cube = isop.next()) != bddfalse)
                      res->new_edge(src_state,
                                    new_state({bdd_exist(cube, all_dots), bddfalse, bddtrue}),
                                    letter);
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



  SPOT_API twa_graph_ptr
  slaa_to_sdba(const_twa_graph_ptr aut, bool cutdet)
  {
    slaa_to_sdba_runner runner(aut, cutdet);
    return runner.run();
  }

}
