// -*- coding: utf-8 -*-
// Copyright (C) 2009, 2010, 2012, 2013 Laboratoire de Recherche et
// Developpement de l'Epita (LRDE).
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

#ifndef SPOT_TGBAALGOS_SCCFILTER_HH
# define SPOT_TGBAALGOS_SCCFILTER_HH

#include "tgba/tgba.hh"

namespace spot
{
  class scc_map;

  /// \brief Prune unaccepting SCCs and remove superfluous acceptance
  /// conditions.
  ///
  /// This function will explore the SCCs of the automaton and remove
  /// dead SCCs (i.e. SCC that are not accepting, and those with no
  /// exit path leading to an accepting SCC).
  ///
  /// Additionally, this will try to remove useless acceptance
  /// conditions.  This operation may diminish the number of
  /// acceptance condition of the automaton (for instance when two
  /// acceptance conditions are always used together we only keep one)
  /// but it will never remove all acceptance conditions, even if it
  /// would be OK to have zero.
  ///
  /// Acceptance conditions on transitions going to non-accepting SCC
  /// are all removed.  Acceptance conditions going to an accepting
  /// SCC and coming from another SCC are only removed if \a
  /// remove_all_useless is set.  The default value of \a
  /// remove_all_useless is \c false because some algorithms (like the
  /// degeneralization) will work better if transitions going to an
  /// accepting SCC are accepting.
  ///
  /// If \a given_sm is supplied, the function will use its result
  /// without computing a map of its own.
  ///
  /// If \a susp is different from bddtrue, it should be a conjunction
  /// of (positive) variables to be removed from transitions going to
  /// non-accepting SCCs.  If early_susp is false, the previous
  /// variable are also removed from transitions entering an accepting
  /// SCC.  ignored is a conjunction of positive variables that should
  /// be removed everywhere.
  ///
  /// \warning Calling scc_filter on a TGBA that has the SBA property
  /// (i.e., transitions leaving accepting states are all marked as
  /// accepting) may destroy this property.  Use scc_filter_states()
  /// instead.
  SPOT_API tgba*
  scc_filter(const tgba* aut, bool remove_all_useless = false,
	     scc_map* given_sm = 0, bdd susp = bddtrue,
	     bool early_susp = false, bdd ignored = bddtrue);

  /// \brief Prune unaccepting SCCs.
  ///
  /// This is an abridged version of scc_filter(), that only remove
  /// useless states, without touching at the acceptance conditions.
  ///
  /// Especially, if the input TGBA has the SBA property, (i.e.,
  /// transitions leaving accepting states are all marked as
  /// accepting), then the output TGBA will also have that property.
  SPOT_API tgba*
  scc_filter_states(const tgba* aut, scc_map* given_sm = 0);
}

#endif // SPOT_TGBAALGOS_SCC_HH