// -*- coding: utf-8 -*-
// Copyright (C) 2014 Laboratoire de Recherche et Développement de
// l'Epita.
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

#ifndef SPOT_TGBAALGOS_SCCINFO_HH
# define SPOT_TGBAALGOS_SCCINFO_HH

#include <vector>
#include "tgba/tgbagraph.hh"

namespace spot
{

  class SPOT_API scc_info
  {
  public:
    struct scc_trans
    {
      scc_trans(bdd cond, unsigned dst)
	: cond(cond), dst(dst)
      {
      }
      bdd cond;
      unsigned dst;
    };

    typedef std::vector<scc_trans> scc_succs;

    struct scc_node
    {
      scc_node():
	acc(0U), trivial(true), accepting(false), useful(false)
      {
      }

      scc_node(acc_cond::mark_t acc, bool trivial):
	acc(acc), trivial(trivial), accepting(false), useful(false)
      {
      }

      scc_succs succ;
      acc_cond::mark_t acc;
      std::list<unsigned> states; // States of the component
      bool trivial:1;
      bool accepting:1;
      bool useful:1;
    };

  protected:

    std::vector<unsigned> sccof_;
    std::vector<scc_node> node_;
    const_tgba_digraph_ptr aut_;


    const scc_node& node(unsigned scc) const
    {
      assert(scc < node_.size());
      return node_[scc];
    }

  public:
    scc_info(const_tgba_digraph_ptr aut);

    const_tgba_digraph_ptr get_aut() const
    {
      return aut_;
    }

    unsigned scc_count() const
    {
      return node_.size();
    }

    bool reachable_state(unsigned st) const
    {
      return scc_of(st) != -1U;
    }

    unsigned scc_of(unsigned st) const
    {
      assert(st < sccof_.size());
      return sccof_[st];
    }

    const std::list<unsigned>& states_of(unsigned scc) const
    {
      return node(scc).states;
    }

    unsigned one_state_of(unsigned scc) const
    {
      return states_of(scc).front();
    }

    /// \brief Get number of the SCC containing the initial state.
    unsigned initial() const
    {
      assert(scc_count() - 1 == scc_of(aut_->get_init_state_number()));
      return scc_count() - 1;
    }

    const scc_succs& succ(unsigned scc) const
    {
      return node(scc).succ;
    }

    bool is_trivial(unsigned scc) const
    {
      return node(scc).trivial;
    }

    acc_cond::mark_t acc(unsigned scc) const
    {
      return node(scc).acc;
    }

    bool is_accepting_scc(unsigned scc) const
    {
      return node(scc).accepting;
    }

    bool is_useful_scc(unsigned scc) const
    {
      return node(scc).useful;
    }

    bool is_useful_state(unsigned st) const
    {
      return reachable_state(st) && node(scc_of(st)).useful;
    }

    /// \brief Return the set of all used acceptance combinations, for
    /// each accepting SCC.
    std::vector<std::set<acc_cond::mark_t>> used_acc() const;

    std::set<acc_cond::mark_t> used_acc_of(unsigned scc) const;


    std::vector<bool> weak_sccs() const;

    bdd scc_ap_support(unsigned scc) const;
  };


  /// \brief Dump the SCC graph of \a aut on \a out.
  ///
  /// If \a sccinfo is not given, it will be computed.
  SPOT_API std::ostream&
  dump_scc_info_dot(std::ostream& out,
		    const_tgba_digraph_ptr aut, scc_info* sccinfo = nullptr);

}

#endif // SPOT_TGBAALGOS_SCCINFO_HH