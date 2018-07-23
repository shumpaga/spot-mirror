// -*- coding: utf-8 -*-
// Copyright (C) 2018 Laboratoire de Recherche
// et Développement de l'Epita (LRDE).
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

#include "config.h"
#include <iostream>
#include <cstdlib>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>
#include <spot/tl/apcollect.hh>
#include <spot/tl/formula.hh>

static void
syntax(char* prog)
{
  std::cerr << prog << " formula?" << std::endl;
  exit(2);
}

int
main(int argc, char** argv)
{
  if (argc != 2)
    syntax(argv[0]);

  auto ftmp = spot::parse_infix_psl(argv[1]);

  if (ftmp.format_errors(std::cerr))
    return 2;

  spot::formula f = ftmp.f;

  spot::atomic_prop_set s;

  literal_collect(f, &s);

  for (auto l: s)
    if (l.is(spot::op::ap))
      print_psl(std::cout, l) << ' ';
  std::cout << '\n';
  for (auto l: s)
    if (l.is(spot::op::Not))
      print_psl(std::cout, l.get_child_of(spot::op::Not)) << ' ';
  std::cout << '\n';

  return 0;
}
