// -*- coding: utf-8 -*-
// Copyright (C) 2010, 2011, 2013, 2014, 2015 Laboratoire de Recherche
// et Développement de l'Epita (LRDE)
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

#include "strength.hh"
#include "misc/hash.hh"
#include "isweakscc.hh"
#include "mask.hh"
#include <deque>

namespace spot
{
  namespace
  {
    template <bool terminal, bool set = false>
    bool is_type_automaton(const twa_graph_ptr& aut, scc_info* si)
    {
      // Create an scc_info if the user did not give one to us.
      bool need_si = !si;
      if (need_si)
	si = new scc_info(aut);
      si->determine_unknown_acceptance();

      bool is_weak = true;
      bool is_term = true;
      unsigned n = si->scc_count();
      for (unsigned i = 0; i < n; ++i)
	{
	  if (si->is_trivial(i))
	    continue;
	  bool first = true;
	  acc_cond::mark_t m = 0U;
	  for (auto src: si->states_of(i))
	    for (auto& t: aut->out(src))
	      if (si->scc_of(t.dst) == i)
		{
		  if (first)
		    {
		      first = false;
		      m = t.acc;
		    }
		  else if (m != t.acc)
		    {
		      is_weak = false;
		      goto exit;
		    }
		}
	  if (terminal && si->is_accepting_scc(i) && !is_complete_scc(*si, i))
	    {
	      is_term = false;
	      if (!set)
		break;
	    }
	}
    exit:
      if (need_si)
	delete si;
      if (set)
	{
	  if (terminal)
	    aut->prop_terminal(is_term && is_weak);
	  aut->prop_weak(is_weak);
	}
      return is_weak && is_term;
    }
  }

  bool
  is_terminal_automaton(const const_twa_graph_ptr& aut, scc_info* si)
  {
    return (aut->prop_terminal() ||
	    is_type_automaton<true>(std::const_pointer_cast<twa_graph>(aut),
				    si));
  }

  bool
  is_weak_automaton(const const_twa_graph_ptr& aut, scc_info* si)
  {
    return (aut->prop_weak() ||
	    is_type_automaton<false>(std::const_pointer_cast<twa_graph>(aut),
				     si));
  }

  void check_strength(const twa_graph_ptr& aut, scc_info* si)
  {
    is_type_automaton<true, true>(aut, si);
  }

  bool is_safety_mwdba(const const_twa_graph_ptr& aut)
  {
    if (!(aut->acc().is_buchi() || aut->acc().is_tt()))
      throw std::runtime_error
	("is_safety_mwdba() should be called on a Buchi automaton");

    for (auto& t: aut->edges())
      if (!aut->acc().accepting(t.acc))
	return false;
    return true;
  }


  twa_graph_ptr
  decompose_strength(const const_twa_graph_ptr& aut, const char* keep_opt)
  {
    if (keep_opt == nullptr || *keep_opt == 0)
      throw std::runtime_error
	(std::string("option for decompose_strength() should not be empty"));

    enum strength {
      Ignore = 0,
      Terminal = 1,
      WeakStrict = 2,
      Weak = Terminal | WeakStrict,
      Strong = 4,
      Needed = 8,		// Needed SCCs are those that lead to
				// the SCCs we want to keep.
    };
    unsigned char keep = Ignore;
    while (auto c = *keep_opt++)
      switch (c)
	{
	case 's':
	  keep |= Strong;
	  break;
	case 't':
	  keep |= Terminal;
	  break;
	case 'w':
	  keep |= WeakStrict;
	  break;
	default:
	  throw std::runtime_error
	    (std::string("unknown option for decompose_strength(): ") + c);
	}

    acc_cond::mark_t wacc = 0U;	// Acceptance for weak SCCs
    acc_cond::mark_t uacc = 0U;	// Acceptance for "needed" SCCs, that
				// we only want to traverse.

    // If the acceptance condition is always satisfiable, we will
    // consider the automaton has weak (even if that is not the
    // case syntactically) and not output any strong part.
    bool all_accepting = false;
    if (aut->acc().is_tt())
      {
	all_accepting = true;
      }
    else if (aut->acc().uses_fin_acceptance())
      {
	auto p = aut->get_acceptance().unsat_mark();
	if (p.first)
	  uacc = p.second;
	else
	  all_accepting = true;
      }
    if (all_accepting)
      {
	keep &= ~Strong;
	if (keep == Ignore)
	  return nullptr;
      }

    scc_info si(aut);
    si.determine_unknown_acceptance();

    unsigned n = si.scc_count();
    std::vector<unsigned char> want(n, Ignore);
    bool nonempty = false;
    bool strong_seen = false;

    for (unsigned i = 0; i < n; ++i) // SCC are topologically ordered
      {
	if (si.is_accepting_scc(i))
	  {
	    if (all_accepting | is_weak_scc(si, i))
	      {
		if (keep & Weak)
		  {
		    if ((keep & Weak) == Weak)
		      want[i] = Weak;
		    else
		      want[i] = keep &
			(is_complete_scc(si, i) ? Terminal : WeakStrict);
		  }
	      }
	    else
	      {
		want[i] = keep & Strong;
		strong_seen = true;
	      }
	    nonempty |= want[i];
	  }
	// An SCC is needed if one of its successor is.
	for (unsigned j: si.succ(i))
	  if (want[j])
	    {
	      want[i] |= Needed;
	      break;
	    }
      }

    if (!nonempty)
      return nullptr;

    twa_graph_ptr res = make_twa_graph(aut->get_dict());
    res->copy_ap_of(aut);
    res->prop_copy(aut, { true, false, true, false });

    if (keep & Strong)
      res->copy_acceptance_of(aut);
    else
      wacc = res->set_buchi();

    auto fun = [&si, &want, uacc, wacc, keep]
      (unsigned src, bdd& cond, acc_cond::mark_t& acc, unsigned dst)
      {
	if (want[si.scc_of(dst)] == Ignore)
	  {
	    cond = bddfalse;
	    return;
	  }
	if (want[si.scc_of(src)] == Needed)
	  {
	    acc = uacc;
	    return;
	  }
	if (keep & Strong)
	  return;
	acc = wacc;
      };

    transform_accessible(aut, res, fun);

    if (!(keep & Strong))
      {
	res->prop_weak(true);
	if (!(keep & WeakStrict))
	  {
	    assert(keep & Terminal);
	    res->prop_terminal(true);
	  }
      }
    else
      {
	res->prop_weak(!strong_seen);
      }
    return res;
  }
}