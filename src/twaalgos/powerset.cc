// -*- coding: utf-8 -*-
// Copyright (C) 2009, 2010, 2011, 2013, 2014, 2015 Laboratoire de
// Recherche et Développement de l'Epita (LRDE).
// Copyright (C) 2004 Laboratoire d'Informatique de Paris 6 (LIP6),
// département Systèmes Répartis Coopératifs (SRC), Université Pierre
// et Marie Curie.
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

#include <set>
#include <iterator>
#include <vector>
#include "powerset.hh"
#include "misc/hash.hh"
#include "twaalgos/powerset.hh"
#include "twaalgos/sccinfo.hh"
#include "twaalgos/cycles.hh"
#include "twaalgos/gtec/gtec.hh"
#include "twaalgos/product.hh"
#include "twa/bddprint.hh"
#include "twaalgos/gtec/gtec.hh"
#include "twaalgos/sccfilter.hh"
#include "twaalgos/ltl2tgba_fm.hh"
#include "twaalgos/dtgbacomp.hh"
#include "ltlast/unop.hh"
#include "misc/bitvect.hh"
#include "misc/bddlt.hh"

namespace spot
{
  namespace
  {

    static unsigned
    number_of_variables(bdd vin)
    {
      unsigned nap = 0;
      int v = vin.id();
      while (v != 1)
	{
	  v = bdd_high(v);
	  ++nap;
	}
      return nap;
    }

    static power_map::power_state
    bv_to_ps(const bitvect* in)
    {
      power_map::power_state ps;
      unsigned ns = in->size();
      for (unsigned pos = 0; pos < ns; ++pos)
	if (in->get(pos))
	  ps.insert(pos);
      return ps;
    }

    struct bv_hash
    {
      size_t operator()(const bitvect* bv) const
      {
	return bv->hash();
      }
    };

    struct bv_equal
    {
      bool operator()(const bitvect* bvl, const bitvect* bvr) const
      {
	return *bvl == *bvr;
      }
    };
  }

  twa_graph_ptr
  tgba_powerset(const const_twa_graph_ptr& aut, power_map& pm, bool merge)
  {
    bdd allap = bddtrue;
    {
      typedef std::set<bdd, bdd_less_than> sup_map;
      sup_map sup;
      // Record occurrences of all guards
      for (auto& t: aut->edges())
	sup.emplace(t.cond);
      for (auto& i: sup)
	allap &= bdd_support(i);
    }

    unsigned nap = number_of_variables(allap);

    // Call this before aut->num_states(), since it might add a state.
    unsigned init_num = aut->get_init_state_number();

    unsigned ns = aut->num_states();
    assert(ns > 0);

    if ((-1UL / ns) >> nap == 0)
      throw std::runtime_error("too many atomic propositions (or states)");

    // Build a correspondence between conjunctions of APs and unsigned
    // indexes.
    std::vector<bdd> num2bdd;
    num2bdd.reserve(1UL << nap);
    std::map<bdd, unsigned, bdd_less_than> bdd2num;
    bdd all = bddtrue;
    while (all != bddfalse)
      {
	bdd one = bdd_satoneset(all, allap, bddfalse);
	all -= one;
	bdd2num.emplace(one, num2bdd.size());
	num2bdd.push_back(one);
      }

    size_t nc = num2bdd.size();	// number of conditions
    assert(nc == (1UL << nap));

    // An array of bit vectors of size 'ns'.  Each original state is
    // represented by 'nc' bitvectors representing the possible
    // destinations for each condition.
    auto bv = std::unique_ptr<bitvect_array>(make_bitvect_array(ns, ns * nc));

    for (unsigned src = 0; src < ns; ++src)
      {
	size_t base = src * nc;
	for (auto& t: aut->out(src))
	  {
	    bdd all = t.cond;
	    while (all != bddfalse)
	      {
		bdd one = bdd_satoneset(all, allap, bddfalse);
		all -= one;
		unsigned num = bdd2num[one];
		bv->at(base + num).set(t.dst);
	      }
	  }
      }

    typedef power_map::power_state power_state;

    typedef std::unordered_map<bitvect*, int, bv_hash, bv_equal> power_set;
    power_set seen;

    std::vector<const bitvect*>toclean;

    auto res = make_twa_graph(aut->get_dict());
    res->copy_ap_of(aut);

    {
      auto bvi = make_bitvect(ns);
      bvi->set(init_num);
      power_state ps{init_num};
      unsigned num = res->new_state();
      res->set_init_state(num);
      seen[bvi] = num;
      assert(pm.map_.size() == num);
      pm.map_.emplace_back(std::move(ps));
      toclean.push_back(bvi);
    }

    // outgoing map
    auto om = std::unique_ptr<bitvect_array>(make_bitvect_array(ns, nc));

    for (unsigned src_num = 0; src_num < res->num_states(); ++src_num)
      {
	om->clear_all();

	const power_state& src = pm.states_of(src_num);

	for (auto s: src)
	  {
	    size_t base = s * nc;
	    for (unsigned c = 0; c < nc; ++c)
	      om->at(c) |= bv->at(base + c);
	  }
	for (unsigned c = 0; c < nc; ++c)
	  {
	    auto dst = &om->at(c);
	    if (dst->is_fully_clear())
	      continue;
	    auto i = seen.find(dst);
	    unsigned dst_num;
	    if (i != seen.end())
	      {
		dst_num = i->second;
	      }
	    else
	      {
		dst_num = res->new_state();
		auto dst2 = dst->clone();
		seen[dst2] = dst_num;
		toclean.push_back(dst2);
		auto ps = bv_to_ps(dst);
		assert(pm.map_.size() == dst_num);
		pm.map_.emplace_back(std::move(ps));
	      }
	    res->new_edge(src_num, dst_num, num2bdd[c]);
	  }
      }

    for (auto v: toclean)
      delete v;
    if (merge)
      res->merge_edges();
    return res;
  }

  twa_graph_ptr
  tgba_powerset(const const_twa_graph_ptr& aut)
  {
    power_map pm;
    return tgba_powerset(aut, pm);
  }


  namespace
  {

    class fix_scc_acceptance final: protected enumerate_cycles
    {
    public:
      typedef dfs_stack::const_iterator cycle_iter;
      typedef twa_graph_edge_data trans;
      typedef std::set<trans*> edge_set;
      typedef std::vector<edge_set> set_set;
    protected:
      const_twa_graph_ptr ref_;
      power_map& refmap_;
      edge_set reject_;	// set of rejecting edges
      set_set accept_;		// set of cycles that are accepting
      edge_set all_;		// all non rejecting edges
      unsigned threshold_;	// maximum count of enumerated cycles
      unsigned cycles_left_; 	// count of cycles left to explore

    public:
      fix_scc_acceptance(const scc_info& sm, const_twa_graph_ptr ref,
			 power_map& refmap, unsigned threshold)
	: enumerate_cycles(sm), ref_(ref), refmap_(refmap),
	  threshold_(threshold)
      {
      }

      bool fix_scc(const int m)
      {
	reject_.clear();
	accept_.clear();
	cycles_left_ = threshold_;
	run(m);

//	std::cerr << "SCC #" << m << '\n';
//	std::cerr << "REJECT: ";
//	print_set(std::cerr, reject_) << '\n';
//	std::cerr << "ALL: ";
//	print_set(std::cerr, all_) << '\n';
//	for (set_set::const_iterator j = accept_.begin();
//	     j != accept_.end(); ++j)
//	  {
//	    std::cerr << "ACCEPT: ";
//	    print_set(std::cerr, *j) << '\n';
//	  }

	auto acc = aut_->acc().all_sets();
	for (auto i: all_)
	  i->acc = acc;
	return threshold_ != 0 && cycles_left_ == 0;
      }

      bool is_cycle_accepting(cycle_iter begin, edge_set& ts) const
      {
	auto a = std::const_pointer_cast<twa_graph>(aut_);

	// Build an automaton representing this loop.
	auto loop_a = make_twa_graph(aut_->get_dict());
	int loop_size = std::distance(begin, dfs_.end());
	loop_a->new_states(loop_size);
	int n;
	cycle_iter i;
	for (n = 1, i = begin; n <= loop_size; ++n, ++i)
	  {
	    trans* t = &a->edge_data(i->succ);
	    loop_a->new_edge(n - 1, n % loop_size, t->cond);
	    if (reject_.find(t) == reject_.end())
	      ts.insert(t);
	  }
	assert(i == dfs_.end());

	unsigned loop_a_init = loop_a->get_init_state_number();
	assert(loop_a_init == 0);

	// Check if the loop is accepting in the original automaton.
	bool accepting = false;

	// Iterate on each original state corresponding to the
	// start of the loop in the determinized automaton.
	for (auto s: refmap_.states_of(begin->s))
	  {
	    // Check the product between LOOP_A, and ORIG_A starting
	    // in S.
	    if (!product(loop_a, ref_, loop_a_init, s)->is_empty())
	      {
		accepting = true;
		break;
	      }
	  }
	return accepting;
      }

      std::ostream&
      print_set(std::ostream& o, const edge_set& s) const
      {
	o << "{ ";
	for (auto i: s)
	  o << i << ' ';
	o << '}';
	return o;
      }

      virtual bool
      cycle_found(unsigned start) override
      {
	cycle_iter i = dfs_.begin();
	while (i->s != start)
	  ++i;
	edge_set ts;
	bool is_acc = is_cycle_accepting(i, ts);
	do
	  ++i;
	while (i != dfs_.end());
	if (is_acc)
	  {
	    accept_.push_back(ts);
	    all_.insert(ts.begin(), ts.end());
	  }
	else
	  {
	    for (auto t: ts)
	      {
		reject_.insert(t);
		for (auto& j: accept_)
		  j.erase(t);
		all_.erase(t);
	      }
	  }

	// Abort this algorithm if we have seen too much cycles, i.e.,
	// when cycle_left_ *reaches* 0.  (If cycle_left_ == 0, that
	// means we had no limit.)
	return (cycles_left_ == 0) || --cycles_left_;
      }
    };

    static bool
    fix_dba_acceptance(twa_graph_ptr det,
		       const_twa_graph_ptr ref, power_map& refmap,
		       unsigned threshold)
    {
      det->copy_acceptance_of(ref);

      scc_info sm(det);

      unsigned scc_count = sm.scc_count();

      fix_scc_acceptance fsa(sm, ref, refmap, threshold);

      for (unsigned m = 0; m < scc_count; ++m)
	if (!sm.is_trivial(m))
	  if (fsa.fix_scc(m))
	    return true;
      return false;
    }
  }

  twa_graph_ptr
  tba_determinize(const const_twa_graph_ptr& aut,
		  unsigned threshold_states, unsigned threshold_cycles)
  {
    power_map pm;
    // Do not merge edges in the deterministic automaton.  If we
    // add two self-loops labeled by "a" and "!a", we do not want
    // these to be merged as "1" before the acceptance has been fixed.
    auto det = tgba_powerset(aut, pm, false);

    if ((threshold_states > 0)
	&& (pm.map_.size() > aut->num_states() * threshold_states))
      return nullptr;
    if (fix_dba_acceptance(det, aut, pm, threshold_cycles))
      return nullptr;
    det->merge_edges();
    return det;
  }

  twa_graph_ptr
  tba_determinize_check(const twa_graph_ptr& aut,
			unsigned threshold_states,
			unsigned threshold_cycles,
			const ltl::formula* f,
			const_twa_graph_ptr neg_aut)
  {
    if (f == 0 && neg_aut == 0)
      return 0;
    if (aut->num_sets() > 1)
      return 0;

    auto det = tba_determinize(aut, threshold_states, threshold_cycles);

    if (!det)
      return nullptr;

    if (neg_aut == nullptr)
      {
	const ltl::formula* neg_f =
	  ltl::unop::instance(ltl::unop::Not, f->clone());
	neg_aut = ltl_to_tgba_fm(neg_f, aut->get_dict());
	neg_f->destroy();

	// Remove useless SCCs.
	neg_aut = scc_filter(neg_aut, true);
      }

    if (product(det, neg_aut)->is_empty())
      // Complement the DBA.
      if (product(aut, dtgba_complement(det))->is_empty())
	// Finally, we are now sure that it was safe
	// to determinize the automaton.
	return det;

    return aut;
  }
}