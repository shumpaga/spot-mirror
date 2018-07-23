// -*- coding: utf-8 -*-
// Copyright (C) 2017 Laboratoire de Recherche et Développement de
// l'Epita (LRDE)
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

#include <memory>
#include <ltdl.h>

namespace spot
{
  ////////////////////////////////////////////////////////////////////////
  // spins interface

  typedef struct transition_info
  {
    int* labels; // edge labels, NULL, or pointer to the edge label(s)
    int  group;  // holds transition group or -1 if unknown
  } transition_info_t;

  typedef void (*TransitionCB)(void *ctx,
                               transition_info_t *transition_info,
                               int *dst);

  struct spins_interface
  {
    // handle to the dynamic library. The variable is of type lt_dlhandle, but
    // we need this trick since we cannot put ltdl.h in public headers
    void* handle;

    // General functions
    void (*get_initial_state)(void *to);
    int (*have_property)();
    int (*get_successors)(void* m, int *in, TransitionCB, void *arg);
    int (*get_state_size)();
    const char* (*get_state_variable_name)(int var);
    int (*get_state_variable_type)(int var);
    int (*get_type_count)();
    const char* (*get_type_name)(int type);
    int (*get_type_value_count)(int type);
    const char* (*get_type_value_name)(int type, int value);

    // Transitions functions
    int (*get_transition_count)();
    int* (*get_transition_read_dependencies)(int t);
    int* (*get_transition_write_dependencies)(int t);

    // Guards functions
    int* (*get_guards)(int t);
    int (*get_guard_count)();
    int* (*get_guard_variables_matrix)(int g);
    int* (*get_guard_nes_matrix)(int g);
    int* (*get_guard_may_be_coenabled_matrix)(int g);
    int (*get_guard)(void* model, int g, const int* src);

    // get_effects function is needed for transparent transitions computing.
    int (*get_effects)(int t, void (*callback)(const char* var_name, const char* str, void* arg), void* arg);

    ~spins_interface()
    {
      lt_dlhandle h = (lt_dlhandle) handle;
      if (h)
        lt_dlclose(h);
      lt_dlexit();
    }
  };

  using spins_interface_ptr = std::shared_ptr<const spins_interface>;
}
