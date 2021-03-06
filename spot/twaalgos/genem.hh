// -*- coding: utf-8 -*-
// Copyright (C) 2017, 2018 Laboratoire de Recherche et Developpement
// de l'Epita (LRDE).
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

#include <spot/twaalgos/sccinfo.hh>

namespace spot
{
  /// \ingroup emptiness_check_algorithms
  /// \brief Emptiness check of on automaton, for any acceptance condition.
  SPOT_API bool
  generic_emptiness_check(const const_twa_graph_ptr& aut);

  /// \ingroup emptiness_check_algorithms
  /// \brief Emptiness check of one SCC, for any acceptance condition.
  SPOT_API bool
  generic_emptiness_check_for_scc(const scc_info& si, unsigned scc);
}
