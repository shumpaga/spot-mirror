// Copyright (C) 2003, 2004, 2013  Laboratoire d'Informatique de Paris 6 (LIP6),
// d�partement Syst�mes R�partis Coop�ratifs (SRC), Universit� Pierre
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

#ifndef SPOT_TGBA_TGBABDDCONCRETEPRODUCT_HH
# define SPOT_TGBA_TGBABDDCONCRETEPRODUCT_HH

#include "tgbabddconcrete.hh"

namespace spot
{
  /// \ingroup tgba_algorithms
  /// \brief Multiplies two spot::tgba_bdd_concrete automata.
  ///
  /// This function builds the resulting product as another
  /// spot::tgba_bdd_concrete automaton.
  SPOT_API tgba_bdd_concrete*
  product(const tgba_bdd_concrete* left, const tgba_bdd_concrete* right);
}

#endif // SPOT_TGBA_TGBABDDCONCRETEPRODUCT_HH