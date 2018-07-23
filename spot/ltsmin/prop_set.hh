// -*- coding: utf-8 -*-
// Copyright (C)  2018 Laboratoire de Recherche et
// Developpement de l'Epita (LRDE)
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

#include <vector>

namespace spot
{
  namespace kripke_cube
  {
    typedef enum {
      OP_EQ_VAR, OP_NE_VAR, OP_LT_VAR, OP_GT_VAR, OP_LE_VAR, OP_GE_VAR,
      VAR_OP_EQ, VAR_OP_NE, VAR_OP_LT, VAR_OP_GT, VAR_OP_LE, VAR_OP_GE,
      VAR_OP_EQ_VAR, VAR_OP_NE_VAR, VAR_OP_LT_VAR,
      VAR_OP_GT_VAR, VAR_OP_LE_VAR, VAR_OP_GE_VAR, VAR_DEAD
    } relop;

    // Structure for complex atomic proposition
    struct one_prop
    {
      using val_t = union { int val; unsigned var; };
      val_t lval;
      relop op;
      val_t rval;

      // For transparent transitions computing; true if the AP is positive in the formula.
      bool pos;
    };

    // Data structure to store complex atomic propositions
    typedef std::vector<one_prop> prop_set;
  }
}
