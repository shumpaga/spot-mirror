// -*- coding: utf-8 -*-
// Copyright (C) 2015 Laboratoire de Recherche et Développement
// de l'Epita.
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

#include <spot/twa/twagraph.hh>

namespace spot
{
  /// \brief Take an automaton with any acceptance condition and return
  /// an equivalent Generalized Büchi automaton.
  SPOT_API twa_graph_ptr
  to_generalized_buchi(const const_twa_graph_ptr& aut);

  /// \brief Convert Streett acceptance into generalized Büchi
  /// acceptance.
  SPOT_API twa_graph_ptr
  streett_to_generalized_buchi(const const_twa_graph_ptr& in);

  /// \brief Convert Streett acceptance into generalized Büchi
  /// only if SPOT_STREET_CONF_MIN is set to a number of pairs
  /// less than the number of pairs used by IN.
  SPOT_API twa_graph_ptr
  streett_to_generalized_buchi_maybe(const const_twa_graph_ptr& in);
}