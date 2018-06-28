// -*- coding: utf-8 -*-
// Copyright (C) 2018 Laboratoire de Recherche
// et Developpement de l'Epita
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
#pragma GCC system_header

#include <mpi.h>
#include <spot/misc/common.hh>

#include <sstream>

namespace spot
{
  namespace mpi
  {
    struct attributes_
    {
      MPI_Comm comm_parent;
      MPI_Comm comm_children;
      MPI_Comm comm_everyone;
      std::stringstream ss_out;
      std::stringstream ss_err;
    };

    class SPOT_API process
    {
      process (void);
      ~process (void);
    };
  }
}
