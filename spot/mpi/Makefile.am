## -*- coding: utf-8 -*-
## Copyright (C) 2018 Laboratoire de
## Recherche et Developpement de l'Epita (LRDE).
##
## This file is part of Spot, a model checking library.
##
## Spot is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Spot is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
## License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) \
              -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CXXFLAGS = $(WARNING_CXXFLAGS)

mpidir = $(pkgincludedir)/mpi

mpi_HEADERS = mpi.hh

lib_LTLIBRARIES = libspotmpi.la
libspotmpi_la_DEPENDENCIES = \
  $(top_builddir)/spot/libspot.la \
  $(top_builddir)/lib/libgnu.la
libspotmpi_la_LIBADD = \
  $(top_builddir)/spot/libspot.la \
  $(top_builddir)/lib/libgnu.la
libspotmpi_la_LDFLAGS = -no-undefined $(SYMBOLIC_LDFLAGS)
libspotmpi_la_SOURCES = mpi.cc


EXTRA_DIST = $(srcdir)/libspotmpi.pc.in
pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = libspotmpi.pc

libspotmpi.pc: $(srcdir)/libspotmpi.pc.in Makefile
	sed -e 's![@]prefix[@]!$(prefix)!g' \
	    -e 's![@]exec_prefix[@]!$(exec_prefix)!g' \
	    -e 's![@]includedir[@]!$(includedir)!g' \
	    -e 's![@]libdir[@]!$(libdir)!g' \
	    -e 's![@]PACKAGE_VERSION[@]!$(PACKAGE_VERSION)!g' \
	$(srcdir)/libspotmpi.pc.in > $@.tmp && mv $@.tmp $@

CLEANFILES = libspotmpi.pc
