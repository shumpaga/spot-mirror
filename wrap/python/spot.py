# -*- coding: utf-8 -*-
# Copyright (C) 2014, 2015  Laboratoire de
# Recherche et Développement de l'Epita (LRDE).
#
# This file is part of Spot, a model checking library.
#
# Spot is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Spot is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from spot_impl import *
import subprocess
import sys
from functools import lru_cache

def setup(**kwargs):
    """Configure Spot for fancy display.

    This is manly useful in IPython.

    Note that this function needs to be called before any automaton is
    displayed.  Afterwards it will have no effect (you should restart
    Python, or the IPython Kernel).

    Keywords arguments:
    bullets -- a Boolean indicating whether to display acceptance conditions
               as UTF8 bullets (default: True)
    fillcolor -- a string indicating the color to use for states
                 (default: '#ffffaa')
    size -- a string giving the width and height of the GraphViz output
            in inches (default: '10.2,5')
    font -- the name of the font to use in the GraphViz output
            (default: 'Lato')
    """
    import os

    s = 'size="{}" node[style=filled,fillcolor="{}"] edge[arrowhead=vee, arrowsize=.7]'
    os.environ['SPOT_DOTEXTRA'] = s.format(kwargs.get('size', '10.2,5'),
                                           kwargs.get('fillcolor', '#ffffaa'))

    bullets = 'B' if kwargs.get('bullets', True) else ''
    d = 'rf({})'.format(kwargs.get('font', 'Lato')) + bullets
    os.environ['SPOT_DOTDEFAULT'] = d

# Global BDD dict so that we do not have to create one in user code.
_bdd_dict = make_bdd_dict()

# Add a small LRU cache so that when we display automata into a
# interactive widget, we avoid some repeated calls to dot for
# identical inputs.
@lru_cache(maxsize=64)
def _str_to_svg(str):
    dotty = subprocess.Popen(['dot', '-Tsvg'],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
    dotty.stdin.write(str)
    res = dotty.communicate()
    return res[0].decode('utf-8')

def _ostream_to_svg(ostr):
    return _str_to_svg(ostr.str().encode('utf-8'))

def _render_automaton_as_svg(a, opt=None):
    ostr = ostringstream()
    print_dot(ostr, a, opt)
    return _ostream_to_svg(ostr)

twa._repr_svg_ = _render_automaton_as_svg
ta._repr_svg_ = _render_automaton_as_svg

def _render_formula_as_svg(a):
    # Load the SVG function only if we need it. This way the bindings
    # can still be used outside of IPython if IPython is not
    # installed.
    from IPython.display import SVG
    ostr = ostringstream()
    print_dot_psl(ostr, a)
    return SVG(_ostream_to_svg(ostr))

def _return_automaton_as_svg(a, opt=None):
    # Load the SVG function only if we need it. This way the bindings
    # can still be used outside of IPython if IPython is not
    # installed.
    from IPython.display import SVG
    return SVG(_render_automaton_as_svg(a, opt))
twa.show = _return_automaton_as_svg
ta.show = _return_automaton_as_svg

def _formula_str_ctor(self, str):
    self.this = parse_formula(str)

def _formula_to_str(self, format = 'spot', parenth = False):
    if format == 'spot':
        return str_psl(self, parenth)
    elif format == 'spin':
        return str_spin_ltl(self, parenth)
    elif format == 'utf8':
        return str_utf8_psl(self, parenth)
    elif format == 'lbt':
        return str_lbt_ltl(self)
    elif format == 'wring':
        return str_wring_ltl(self)
    elif format == 'latex':
        return str_latex_psl(self, parenth)
    elif format == 'sclatex':
        return str_sclatex_psl(self, parenth)
    else:
        raise ValueError("unknown string format: " + format)

formula.__init__ = _formula_str_ctor
formula.to_str = _formula_to_str
formula.show_ast = _render_formula_as_svg

def _twa_to_str(a, format='hoa', opt=None):
    format = format.lower()
    if format == 'hoa':
        ostr = ostringstream()
        print_hoa(ostr, a, opt)
        return ostr.str()
    if format == 'dot':
        ostr = ostringstream()
        print_dot(ostr, a, opt)
        return ostr.str()
    if format == 'spin':
        ostr = ostringstream()
        print_never_claim(ostr, a, opt)
        return ostr.str()
    if format == 'lbtt':
        ostr = ostringstream()
        print_lbtt(ostr, a, opt)
        return ostr.str()
    raise ValueError("unknown string format: " + format)

def _twa_save(a, filename, format='hoa', opt=None, append=False):
    with open(filename, 'a' if append else 'w') as f:
        s = a.to_str(format, opt)
        f.write(s)
        if s[-1] != '\n':
            f.write('\n')
    return a

twa.to_str = _twa_to_str
twa.save = _twa_save

def automata(*filenames):
    """Read automata from a list of filenames.

    The automata can be written in the
    [HOA format](http://adl.github.io/hoaf/), as
    [never claims](http://spinroot.com/spin/Man/never.html),
    or in [LBTT's format]
    (http://www.tcs.hut.fi/Software/lbtt/doc/html/Format-for-automata.html).

    If an argument ends with a `|`, then this argument is interpreted as
    a shell command, and the output of that command (without the `|`)
    is parsed.

    If an argument contains a newline, then it is interpreted as
    actual contents to be parsed.

    Otherwise, the argument is assumed to be a filename.
    """

    for filename in filenames:
        try:
            p = None
            if filename[-1] == '|':
                proc = subprocess.Popen(filename[:-1], shell=True,
                                        stdout=subprocess.PIPE)
                p = automaton_stream_parser(proc.stdout.fileno(),
                                            filename, True)
            elif '\n' in filename:
                proc = None
                p = automaton_stream_parser(filename, "<string>", True)
            else:
                proc = None
                p = automaton_stream_parser(filename, True)
            a = True
            while a:
                # This returns None when we reach the end of the file.
                a = p.parse_strict(_bdd_dict)
                if a:
                    yield a
        finally:
            # Make sure we destroy the parser (p) and the subprocess
            # (prop) in the correct order...
            del p
            if proc != None:
                if not a:
                    # We reached the end of the stream.  Wait for the
                    # process to finish, so that we can its exit code.
                    ret = proc.wait()
                else:
                    # if a != None, we probably got there through an
                    # exception, and the subprocess my still be
                    # running.  Check if an exit status is available
                    # just in case.
                    ret = proc.poll()
                del proc
                if ret:
                    raise RuntimeError("Command {} exited with exit status {}"
                                       .format(filename[:-1], ret))
    return

def automaton(filename):
    """Read a single automaton from a file.

    See `spot.automata()` for a list of supported formats."""
    try:
        return next(automata(filename))
    except StopIteration:
        raise RuntimeError("Failed to read automaton from {}".format(filename))

def translate(formula, *args):
    """Translate a formula into an automaton.

    Keep in mind that pref expresses just a preference that may not be
    satisfied.

    The optional arguments should be strings among the following:
    - at most one in 'TGBA', 'BA', or 'Monitor'
      (type of automaton to build)
    - at most one in 'Small', 'Deterministic', 'Any'
      (preferred characteristics of the produced automaton)
    - at most one in 'Low', 'Medium', 'High'
      (optimization level)
    - any combination of 'Complete', 'Unambiguous', and
      'StateBasedAcceptance' (or 'SBAcc' for short)

    The default correspond to 'tgba', 'small' and 'high'.
    """

    type_ = None
    pref_ = None
    optm_ = None
    comp_ = 0
    unam_ = 0
    sbac_ = 0

    def type_set(val):
        nonlocal type_
        if type_ != None and type_ != val:
            raise ValueError("type cannot be both {} and {}"
                             .format(type_, val))
        elif val == 'tgba':
            type_ = postprocessor.TGBA
        elif val == 'ba':
            type_ = postprocessor.BA
        else:
            assert(val == 'monitor')
            type_ = postprocessor.Monitor

    def pref_set(val):
        nonlocal pref_
        if pref_ != None and pref_ != val:
            raise ValueError("preference cannot be both {} and {}"
                             .format(pref_, val))
        elif val == 'small':
            pref_ = postprocessor.Small
        elif val == 'deterministic':
            pref_ = postprocessor.Deterministic
        else:
            assert(val == 'any')
            pref_ = postprocessor.Any

    def optm_set(val):
        nonlocal optm_
        if optm_ != None and optm_ != val:
            raise ValueError("optimization level cannot be both {} and {}"
                             .format(optm_, val))
        if val == 'high':
            optm_ = postprocessor.High
        elif val.startswith('med'):
            optm_ = postprocessor.Medium
        else:
            assert(val == 'low')
            optm_ = postprocessor.Low

    def misc_set(val):
        nonlocal comp_, unam_, sbac_
        if val == 'complete':
            comp_ = postprocessor.Complete
        elif val == 'sbacc' or val == 'state-based-acceptance':
            sbac_ = postprocessor.SBAcc
        else:
            assert(val == 'unambiguous')
            unam_ = postprocessor.Unambiguous

    options = {
        'tgba': type_set,
        'ba': type_set,
        'monitor': type_set,
        'small': pref_set,
        'deterministic': pref_set,
        'any': pref_set,
        'high': optm_set,
        'medium': optm_set,
        'low': optm_set,
        'complete': misc_set,
        'unambiguous': misc_set,
        'statebasedacceptance': misc_set,
        'sbacc': misc_set,
    }

    for arg in args:
        arg = arg.lower()
        fn = options.get(arg)
        if fn:
            fn(arg)
        else:
            # arg is not an know option, but maybe it is a prefix of
            # one of them
            compat = []
            f = None
            for key, fn in options.items():
                if key.startswith(arg):
                    compat.append(key)
                    f = fn
            lc = len(compat)
            if lc == 1:
                f(compat[0])
            elif lc < 1:
                raise ValueError("unknown option '{}'".format(arg))
            else:
                raise ValueError("ambiguous option '{}' is prefix of {}"
                                 .format(arg, str(compat)))

    if type_ == None:
        type_ = postprocessor.TGBA
    if pref_ == None:
        pref_ = postprocessor.Small
    if optm_ == None:
        optm_ = postprocessor.High

    if type(formula) == str:
        formula = parse_formula(formula)


    a = translator(_bdd_dict)
    a.set_type(type_)
    a.set_pref(pref_ | comp_ | unam_ | sbac_)
    a.set_level(optm_)

    return a.run(formula)

formula.translate = translate

# Wrapper around a formula iterator to which we add some methods of formula
# (using _addfilter and _addmap), so that we can write things like
# formulas.simplify().is_X_free().
class formulaiterator:
    def __init__(self, formulas):
        self._formulas = formulas

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._formulas)

# fun shoud be a predicate and should be a method of formula.
# _addfilter adds this predicate as a filter whith the same name in
# formulaiterator.
def _addfilter(fun):
    def filtf(self, *args, **kwargs):
        it = filter(lambda f: getattr(f, fun)(*args, **kwargs), self)
        return formulaiterator(it)
    def nfiltf(self, *args, **kwargs):
        it = filter(lambda f: not getattr(f, fun)(*args, **kwargs), self)
        return formulaiterator(it)
    setattr(formulaiterator, fun, filtf)
    if fun[:3] == 'is_':
        notfun = fun[:3] + 'not_' + fun[3:]
    elif fun[:4] == 'has_':
        notfun = fun[:4] + 'no_' + fun[4:]
    else:
        notfun = 'not_' + fun
    setattr(formulaiterator, fun, filtf)
    setattr(formulaiterator, notfun, nfiltf)

# fun should be a function taking a formula as its first parameter and returning
# a formula.
# _addmap adds this function as a method of formula and formulaiterator.
def _addmap(fun):
    def mapf(self, *args, **kwargs):
        return formulaiterator(map(lambda f: getattr(f, fun)(*args, **kwargs),
self))
    setattr(formula, fun, lambda self, *args, **kwargs: globals()[fun](self,
    *args, **kwargs))
    setattr(formulaiterator, fun, mapf)

def randltl(ap, n = -1, **kwargs):
    """Generate random formulas.

    Returns a random formula iterator.

    ap: the number of atomic propositions used to generate random formulas.

    n: number of formulas to generate, or unbounded if n < 0.

    **kwargs:
    seed: seed for the random number generator (0).
    output: can be 'ltl', 'psl', 'bool' or 'sere' ('ltl').
    allow_dups: allow duplicate formulas (False).
    tree_size: tree size of the formulas generated, before mandatory
    simplifications (15)
    boolean_priorities: set priorities for Boolean formulas.
    ltl_priorities: set priorities for LTL formulas.
    sere_priorities: set priorities for SERE formulas.
    dump_priorities: show current priorities, do not generate any formula.
    simplify:
      0           No rewriting
      1           basic rewritings and eventual/universal rules
      2           additional syntactic implication rules
      3 (default) better implications using containment
    """
    opts = option_map()
    output_map = {
        "ltl" : OUTPUTLTL,
        "psl" : OUTPUTPSL,
        "bool" : OUTPUTBOOL,
        "sere" : OUTPUTSERE
    }

    if isinstance(ap, list):
        aprops = atomic_prop_set()
        e = default_environment.instance()
        for elt in ap:
            aprops.insert(is_atomic_prop(e.require(elt)))
        ap = aprops
    ltl_priorities = kwargs.get("ltl_priorities", None)
    sere_priorities = kwargs.get("sere_priorities", None)
    boolean_priorities = kwargs.get("boolean_priorities", None)
    output = output_map[kwargs.get("output", "ltl")]
    opts.set("output", output)
    opts.set("seed", kwargs.get("seed", 0))
    tree_size = kwargs.get("tree_size", 15)
    if isinstance(tree_size, tuple):
        tree_size_min, tree_size_max = tree_size
    else:
        tree_size_min = tree_size_max = tree_size
    opts.set("tree_size_min", tree_size_min)
    opts.set("tree_size_max", tree_size_max)
    opts.set("unique", not kwargs.get("allow_dups", False))
    opts.set("wf", kwargs.get("weak_fairness", False))
    simpl_level = kwargs.get("simplify", 0)
    if simpl_level > 3 or simpl_level < 0:
        sys.stderr.write('invalid simplification level: ' + simpl_level)
        return
    opts.set("simplification_level", simpl_level)

    rg = randltlgenerator(ap, opts, ltl_priorities, sere_priorities,
        boolean_priorities)

    dump_priorities = kwargs.get("dump_priorities", False)
    if dump_priorities:
        dumpstream = ostringstream()
        if output == OUTPUTLTL:
            print('Use argument ltl_priorities=STRING to set the following ' \
                    'LTL priorities:\n')
            rg.dump_ltl_priorities(dumpstream)
            print(dumpstream.str())
        elif output == OUTPUTBOOL:
            print('Use argument boolean_priorities=STRING to set the ' \
                    'following Boolean formula priorities:\n')
            rg.dump_bool_priorities(dumpstream)
            print(dumpstream.str())
        elif output == OUTPUTPSL or output == OUTPUTSERE:
            if output != OUTPUTSERE:
                print('Use argument ltl_priorities=STRING to set the following ' \
                        'LTL priorities:\n')
                rg.dump_psl_priorities(dumpstream)
                print(dumpstream.str())
            print('Use argument sere_priorities=STRING to set the following ' \
                    'SERE priorities:\n')
            rg.dump_sere_priorities(dumpstream)
            print(dumpstream.str())
            print('Use argument boolean_priorities=STRING to set the ' \
                    'following Boolean formula priorities:\n')
            rg.dump_sere_bool_priorities(dumpstream)
            print(dumpstream.str())
        else:
            sys.stderr.write("internal error: unknown type of output")
        return

    def _randltlgenerator(rg):
        i = 0
        while i != n:
            f = rg.next()
            if f is None:
                sys.stderr.write("Warning: could not generate a new unique formula " \
                "after " + str(MAX_TRIALS) + " trials.\n")
                yield None
            else:
                yield f
                i += 1
    return formulaiterator(_randltlgenerator(rg))

def simplify(f, **kwargs):
    level = kwargs.get('level', None)
    if level is not None:
        return ltl_simplifier(ltl_simplifier_options(level)).simplify(f)

    basics = kwargs.get('basics', True)
    synt_impl = kwargs.get('synt_impl', True)
    event_univ = kwargs.get('event_univ', True)
    containment_checks = kwargs.get('containment_checks', False)
    containment_checks_stronger = kwargs.get('containment_checks_stronger', False)
    nenoform_stop_on_boolean = kwargs.get('nenoform_stop_on_boolean', False)
    reduce_size_strictly = kwargs.get('reduce_size_strictly', False)
    boolean_to_isop = kwargs.get('boolean_to_isop', False)
    favor_event_univ = kwargs.get('favor_event_univ', False)

    simp_opts = ltl_simplifier_options(basics,
                                       synt_impl,
                                       event_univ,
                                       containment_checks,
                                       containment_checks_stronger,
                                       nenoform_stop_on_boolean,
                                       reduce_size_strictly,
                                       boolean_to_isop,
                                       favor_event_univ)
    return ltl_simplifier(simp_opts).simplify(f)

for fun in dir(formula):
    if (callable(getattr(formula, fun)) and
        (fun[:3] == 'is_' or fun[:4] == 'has_')):
        _addfilter(fun)

for fun in ['remove_x', 'get_literal', 'relabel', 'relabel_bse',
            'simplify', 'unabbreviate_ltl']:
    _addmap(fun)