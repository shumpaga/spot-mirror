// -*- coding: utf-8 -*-
// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015  Laboratoire de
// Recherche et Développement de l'Epita (LRDE).
// Copyright (C) 2003, 2004, 2005, 2006  Laboratoire d'Informatique
// de Paris 6 (LIP6), département Systèmes Répartis Coopératifs (SRC),
// Université Pierre et Marie Curie.
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

%{
  // Workaround for SWIG 2.0.2 using ptrdiff_t but not including cstddef.
  // It matters with g++ 4.6.
#include <cstddef>
%}

%module(directors="1") spot_impl

%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_list.i"
%include "std_set.i"
%include "std_map.i"
%include "exception.i"
%include "typemaps.i"

 // git grep 'typedef.*std::shared_ptr' | grep -v const |
 //   sed 's/.*<\(.*\)>.*/%shared_ptr(spot::\1)/g'
%shared_ptr(spot::dstar_aut)
%shared_ptr(spot::parsed_aut)
%shared_ptr(spot::fair_kripke)
%shared_ptr(spot::kripke)
%shared_ptr(spot::kripke_explicit)
%shared_ptr(spot::kripke)
%shared_ptr(spot::ta)
%shared_ptr(spot::ta_explicit)
%shared_ptr(spot::ta_product)
%shared_ptr(spot::tgta)
%shared_ptr(spot::tgta_explicit)
%shared_ptr(spot::bdd_dict)
%shared_ptr(spot::twa)
%shared_ptr(spot::twa_graph)
%shared_ptr(spot::twa_product)
%shared_ptr(spot::twa_product_init)
%shared_ptr(spot::taa_tgba)
%shared_ptr(spot::taa_tgba_string)
%shared_ptr(spot::taa_tgba_formula)
%shared_ptr(spot::twa_safra_complement)
%shared_ptr(spot::tgba_run)
%shared_ptr(spot::emptiness_check_result)
%shared_ptr(spot::emptiness_check)
%shared_ptr(spot::emptiness_check_instantiator)
%shared_ptr(spot::tgbasl)

namespace std {
   %template(liststr) list<string>;
};

%import "buddy.i"

%{
#include <iostream>
#include <fstream>
#include <sstream>
#include <signal.h>

#include "misc/common.hh"
#include "misc/version.hh"
#include "misc/minato.hh"
#include "misc/optionmap.hh"
#include "misc/random.hh"

#include "ltlast/formula.hh"
#include "ltlast/atomic_prop.hh"
#include "ltlast/binop.hh"
#include "ltlast/constant.hh"
#include "ltlast/multop.hh"
#include "ltlast/unop.hh"
#include "ltlast/visitor.hh"

#include "ltlenv/environment.hh"
#include "ltlenv/defaultenv.hh"

#include "ltlparse/public.hh"

#include "twa/bdddict.hh"

#include "ltlvisit/apcollect.hh"
#include "ltlvisit/dot.hh"
#include "ltlvisit/dump.hh"
#include "ltlvisit/nenoform.hh"
#include "ltlvisit/print.hh"
#include "ltlvisit/simplify.hh"
#include "ltlvisit/unabbrev.hh"
#include "ltlvisit/randomltl.hh"
#include "ltlvisit/length.hh"
#include "ltlvisit/remove_x.hh"
#include "ltlvisit/relabel.hh"

#include "twa/bddprint.hh"
#include "twa/fwd.hh"
#include "twa/acc.hh"
#include "twa/twa.hh"
#include "twa/taatgba.hh"
#include "twa/twaproduct.hh"

#include "twaalgos/cleanacc.hh"
#include "twaalgos/dot.hh"
#include "twaalgos/degen.hh"
#include "twaalgos/dupexp.hh"
#include "twaalgos/emptiness.hh"
#include "twaalgos/gtec/gtec.hh"
#include "twaalgos/lbtt.hh"
#include "twaalgos/ltl2taa.hh"
#include "twaalgos/ltl2tgba_fm.hh"
#include "twaalgos/compsusp.hh"
#include "twaalgos/magic.hh"
#include "twaalgos/minimize.hh"
#include "twaalgos/neverclaim.hh"
#include "twaalgos/remfin.hh"
#include "twaalgos/safety.hh"
#include "twaalgos/sccfilter.hh"
#include "twaalgos/stats.hh"
#include "twaalgos/isdet.hh"
#include "twaalgos/simulation.hh"
#include "twaalgos/postproc.hh"
#include "twaalgos/product.hh"
#include "twaalgos/stutter.hh"
#include "twaalgos/translate.hh"
#include "twaalgos/hoa.hh"
#include "twaalgos/dtgbasat.hh"
#include "twaalgos/relabel.hh"

#include "parseaut/public.hh"

#include "ta/ta.hh"
#include "ta/tgta.hh"
#include "ta/taexplicit.hh"
#include "ta/tgtaexplicit.hh"
#include "taalgos/tgba2ta.hh"
#include "taalgos/dot.hh"
#include "taalgos/stats.hh"
#include "taalgos/minimize.hh"

using namespace spot::ltl;
using namespace spot;
%}

// For spot::emptiness_check_instantiator::construct and any other
// function that return errors via a "char **err" argument.
%typemap(in, numinputs=0) char** OUTPUT (char* temp)
  "$1 = &temp;";
%typemap(argout) char** OUTPUT {
  PyObject *obj = SWIG_FromCharPtr(*$1);
  if (!$result) {
    $result = obj;
  //# If the function returns null_ptr (i.e. Py_None), we
  //# don't want to override it with OUTPUT as in the
  //# default implementation of t_output_helper.
  // } else if ($result == Py_None) {
  //   Py_DECREF($result);
  //   $result = obj;
  } else {
    if (!PyList_Check($result)) {
      PyObject *o2 = $result;
      $result = PyList_New(1);
      PyList_SetItem($result, 0, o2);
    }
    PyList_Append($result,obj);
    Py_DECREF(obj);
  }

 };
%apply char** OUTPUT { char** err };

%exception {
  try {
    $action
  }
  catch (const spot::parse_error& e)
  {
    std::string er("\n");
    er += e.what();
    SWIG_exception(SWIG_SyntaxError, er.c_str());
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::runtime_error& e)
  {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// False and True cannot be redefined in Python3, even in a class.
// Spot uses these in an enum of the constant class.
%rename(FalseVal) False;
%rename(TrueVal) True;

%include "misc/common.hh"
%include "misc/version.hh"
%include "misc/minato.hh"
%include "misc/optionmap.hh"
%include "misc/random.hh"

%include "ltlast/formula.hh"
%include "ltlast/atomic_prop.hh"
%include "ltlast/binop.hh"
%include "ltlast/constant.hh"
%include "ltlast/multop.hh"
%include "ltlast/unop.hh"
%include "ltlast/visitor.hh"

namespace std {
  %template(atomic_prop_set)    set<const spot::ltl::atomic_prop*,
                                    spot::ltl::formula_ptr_less_than>;

  %template(mapff) map<const spot::ltl::formula*, const spot::ltl::formula*,
                       spot::ltl::formula_ptr_less_than>;
}

%{
  namespace swig {
    template <>  struct traits<atomic_prop> {
      typedef pointer_category category;
      static const char* type_name() { return "spot::ltl::atomic_prop"; }
    };
    template <>  struct traits<formula> {
      typedef pointer_category category;
      static const char* type_name() { return "spot::ltl::formula"; }
    };
  }
%}


%include "ltlenv/environment.hh"
%include "ltlenv/defaultenv.hh"

%include "ltlparse/public.hh"

 /* these must come before apcollect.hh */
%include "twa/bdddict.hh"
%include "twa/bddprint.hh"
%include "twa/fwd.hh"
%feature("flatnested") spot::acc_cond::mark_t;
%feature("flatnested") spot::acc_cond::acc_code;
%apply bool* OUTPUT {bool& max, bool& odd};
%include "twa/acc.hh"
%include "twa/twa.hh"

%include "ltlvisit/apcollect.hh"
%include "ltlvisit/dot.hh"
%include "ltlvisit/dump.hh"
%include "ltlvisit/nenoform.hh"
%include "ltlvisit/print.hh"
%include "ltlvisit/simplify.hh"
%include "ltlvisit/unabbrev.hh"
%include "ltlvisit/randomltl.hh"
%include "ltlvisit/length.hh"
%include "ltlvisit/remove_x.hh"
%include "ltlvisit/relabel.hh"

// Help SWIG with namespace lookups.
#define ltl spot::ltl
%include "twa/taatgba.hh"
%include "twa/twaproduct.hh"
%include "twa/twagraph.hh"

// Should come after the definition of twa_graph

%include "twaalgos/cleanacc.hh"
%include "twaalgos/degen.hh"
%include "twaalgos/dot.hh"
%include "twaalgos/dupexp.hh"
%include "twaalgos/emptiness.hh"
%include "twaalgos/gtec/gtec.hh"
%include "twaalgos/lbtt.hh"
%include "twaalgos/ltl2taa.hh"
%include "twaalgos/ltl2tgba_fm.hh"
%include "twaalgos/compsusp.hh"
%include "twaalgos/magic.hh"
%include "twaalgos/minimize.hh"
%include "twaalgos/neverclaim.hh"
%include "twaalgos/safety.hh"
%include "twaalgos/remfin.hh"
%include "twaalgos/sccfilter.hh"
%include "twaalgos/stats.hh"
%include "twaalgos/isdet.hh"
%include "twaalgos/simulation.hh"
%include "twaalgos/postproc.hh"
%include "twaalgos/product.hh"
%include "twaalgos/stutter.hh"
%include "twaalgos/translate.hh"
%include "twaalgos/hoa.hh"
%include "twaalgos/dtgbasat.hh"
%include "twaalgos/relabel.hh"

%include "parseaut/public.hh"

%include "ta/ta.hh"
%include "ta/tgta.hh"
%include "ta/taexplicit.hh"
%include "ta/tgtaexplicit.hh"
%include "taalgos/tgba2ta.hh"
%include "taalgos/dot.hh"
%include "taalgos/stats.hh"
%include "taalgos/minimize.hh"


#undef ltl

%extend spot::ltl::formula {

  // When comparing formula, make sure Python compare our
  // pointers, not the pointers to its wrappers.

  // __cmp__ is for Python 2.0
  int __cmp__(const spot::ltl::formula* b) { return self - b; }
  // These are for Python 2.1+ or 3.x.  They more closely match
  // the logic in Spot.
  bool __lt__(const spot::ltl::formula* b)
  { spot::ltl::formula_ptr_less_than lt; return lt(self, b); }
  bool __le__(const spot::ltl::formula* b)
  { spot::ltl::formula_ptr_less_than lt; return !lt(b, self); }
  bool __eq__(const spot::ltl::formula* b) { return self == b; }
  bool __ne__(const spot::ltl::formula* b) { return self != b; }
  bool __gt__(const spot::ltl::formula* b)
  { spot::ltl::formula_ptr_less_than lt; return lt(b, self); }
  bool __ge__(const spot::ltl::formula* b)
  { spot::ltl::formula_ptr_less_than lt; return !lt(self, b); }

  size_t __hash__() { return self->hash(); }

  std::string __repr__() { return spot::ltl::str_psl(self); }
  std::string _repr_latex_()
  {
    return std::string("$") + spot::ltl::str_sclatex_psl(self) + '$';
  }
  std::string __str__() { return spot::ltl::str_psl(self); }
}

%extend spot::acc_cond::acc_code {
  std::string __repr__()
  {
    std::ostringstream os;
    os << *self;
    return os.str();
  }

  std::string __str__()
  {
    std::ostringstream os;
    os << *self;
    return os.str();
  }
}

%nodefaultctor std::ostream;
namespace std {
  class ostream {};
  class ofstream : public ostream
  {
  public:
     ofstream(const char *fn);
     ~ofstream();
  };
  class ostringstream : public ostream
  {
  public:
     ostringstream();
     std::string str() const;
     ~ostringstream();
  };
}

%inline %{
spot::ltl::parse_error_list
empty_parse_error_list()
{
  parse_error_list l;
  return l;
}

spot::parse_aut_error_list
empty_parse_aut_error_list()
{
  parse_aut_error_list l;
  return l;
}

spot::twa_graph_ptr
ensure_digraph(const spot::twa_ptr& a)
{
  auto aa = std::dynamic_pointer_cast<spot::twa_graph>(a);
  if (aa)
    return aa;
  return spot::make_twa_graph(a, spot::twa::prop_set::all());
}

std::ostream&
get_cout()
{
  return std::cout;
}

void
nl_cout()
{
  std::cout << std::endl;
}

std::ostream&
get_cerr()
{
  return std::cerr;
}

void
nl_cerr()
{
  std::cerr << std::endl;
}

void
print_on(std::ostream& on, const std::string& what)
{
  on << what;
}

int
unblock_signal(int signum)
{
  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, signum);
  return sigprocmask(SIG_UNBLOCK, &set, 0);
}

%}

%extend spot::ltl::parse_error_list {

bool
__nonzero__()
{
  return !self->empty();
}

bool
__bool__()
{
  return !self->empty();
}

}

%extend spot::parse_aut_error_list {

bool
__nonzero__()
{
  return !self->empty();
}

bool
__bool__()
{
  return !self->empty();
}

}