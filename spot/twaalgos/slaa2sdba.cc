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
#include <sstream>

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
    typedef std::pair<bdd, bdd> pair_t;

    struct pair_lt
    {
      bool operator()(const pair_t& left, const pair_t& right) const
      {
        if (left.first.id() < right.first.id())
          return true;
        if (left.first.id() > right.first.id())
          return false;
        return left.second.id() < right.second.id();
      }
    };

    class slaa_to_sdba_runner final
    {
      const_twa_graph_ptr aut_;
      // We use n BDD variables to represent the n states.
      // The state i is BDD variable state_base_+i.
      unsigned state_base_;
      bdd all_states_;
      bool cutdet_;
      std::map<pair_t, unsigned, pair_lt> pair_to_unsigned;
      std::vector<pair_t> unsigned_to_pair;

    public:
      slaa_to_sdba_runner(const_twa_graph_ptr aut, bool cutdet)
        : aut_(aut), cutdet_(cutdet)
      {
      }

      ~slaa_to_sdba_runner()
      {
        aut_->get_dict()->unregister_all_my_variables(this);
      }

      bdd state_bdd(unsigned num)
      {
        return bdd_ithvar(state_base_ + num);
      }

      unsigned state_num(bdd b)
      {
        unsigned num = bdd_var(b) - state_base_;
        assert(num < aut_->num_states());
        return num;
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
        while ((cube = isop.next()) != bddfalse)
          {
            if (!outer_first)
              out << "+";
            bool inner_first = true;
            while (cube != bddtrue)
              {
                if (!inner_first)
                  out << "·";
                out << state_num(cube);
                assert(bdd_low(cube) == bddfalse); // expect only positive vars
                cube = bdd_high(cube);
                inner_first = false;
              }
            outer_first = false;
          }
        return out.str();
      }

      std::string print_pair(const pair_t& p)
      {
        std::string res =  print_state_formula(p.first);
        if (p.second != bddfalse)
          res += "\n" + print_state_formula(p.second);
        return res;
      }


      twa_graph_ptr run()
      {
        auto d = aut_->get_dict();

        twa_graph_ptr res = make_twa_graph(d);
        res->copy_ap_of(aut_);
        res->prop_copy(aut_, {false, false, false, false, true, true});
        res->set_buchi();

        // We need one BDD variable per possible output state.  If
        unsigned ns = aut_->num_states();
        state_base_ = d->register_anonymous_variables(ns, this);
        bdd all_states = bddtrue;
        for (unsigned s = 0; s < ns; ++s)
          all_states &= state_bdd(s);
        all_states_ = all_states;


        // Scan for true states.  Since we work with co-Büchi SLAA,
        // any state with a non-rejecting true self-loop is a true
        // state.
        std::vector<bool> true_state(ns, false);
        for (unsigned s = 0; s < ns; ++s)
          for (auto& e: aut_->out(s))
            if (e.dst == s && e.cond == bddtrue && !e.acc)
              {
                true_state[s] = true;
                break;
              }

        // Compute the BDD representation of each state.
        std::vector<bdd> state_as_bdd;
        state_as_bdd.reserve(ns);
        for (unsigned s = 0; s < ns; ++s)
          {
            bdd res = bddfalse;
            if (true_state[s])
              res = bddtrue;
            else
              for (auto& e: aut_->out(s))
                {
                  bdd dest = bddtrue;
                  for (unsigned d: aut_->univ_dests(e.dst))
                    if (!true_state[d])
                      dest &= state_bdd(d);
                  res |= e.cond & dest;
                }
            state_as_bdd.push_back(res);
          }

        // Compute a BDD representing the successors of one
        // state-formula.  The successor representation is a
        // BDD formula containing letters and state numbers.
        auto delta = [&](bdd left)
          {
            bdd res = bddfalse;
            minato_isop isop(left);
            bdd cube;
            while ((cube = isop.next()) != bddfalse)
              {
                bdd inner_res = bddtrue;
                while (cube != bddtrue)
                  {
                    inner_res &= state_as_bdd[state_num(cube)];
                    // expect only positive vars
                    assert(bdd_low(cube) == bddfalse);
                    cube = bdd_high(cube);
                  }
                res |= inner_res;
              }
            return res;
          };

        std::vector<unsigned> todo;

        auto new_state = [&](pair_t s)
          {
            auto p = pair_to_unsigned.emplace(s, 0U);
            if (p.second)
              {
                unsigned_to_pair.push_back(s);
                unsigned q = res->new_state();
                assert(unsigned_to_pair.size() == q + 1);
                todo.push_back(q);
                p.first->second = q;
              }
            return p.first->second;
          };

        new_state({state_bdd(aut_->get_init_state_number()), bddfalse});

        while (!todo.empty())
          {
            unsigned src_state = todo.back();
            pair_t src = unsigned_to_pair[src_state];
            todo.pop_back();
            bdd succs = delta(src.first);
            // decompose successors
            bdd aps = bdd_exist(bdd_support(succs), all_states_);
            bdd all = bddtrue;
            while (all != bddfalse)
              {
                bdd letter = bdd_satoneset(all, aps, bddfalse);
                all -= letter;
                bdd left = bdd_restrict(succs, letter);
                if (cutdet_)
                  {
                    if (left != bddfalse)
                      res->new_edge(src_state, new_state({left, bddfalse}), letter);
                  }
                else
                  {
                    minato_isop isop(left);
                    bdd cube;
                    while ((cube = isop.next()) != bddfalse)
                      res->new_edge(src_state,
                                    new_state({cube, bddfalse}), letter);
                  }
              }
          }

        res->merge_edges();

        auto* names = new std::vector<std::string>;
        res->set_named_prop("state-names", names);
        unsigned output_size = unsigned_to_pair.size();
        names->reserve(output_size);
        for (unsigned s = 0; s < output_size; ++s)
          names->emplace_back(print_pair(unsigned_to_pair[s]));

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
