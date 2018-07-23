// -*- coding: utf-8 -*-
// Copyright (C) 2009, 2010, 2013, 2014, 2016, 2017 Laboratoire de Recherche
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

#include <spot/kripke/fairkripke.hh>
#include <spot/twacube/cube.hh>
#include <memory>

namespace spot
{
  /// \brief This class is a template representation of a Kripke
  /// structure. It is composed of two template parameters: State
  /// represents a state of the Kripke structure, SuccIterator is
  /// an iterator over the (possible) successors of a state.
  ///
  /// Do not delete by hand any states and/or iterator that
  /// are provided by this template class. Specialisations
  /// will handle it.
  template<typename State, typename SuccIterator>
  class SPOT_API kripkecube:
    public std::enable_shared_from_this<kripkecube<State, SuccIterator>>
  {
  public:
    /// \brief Returns the initial State of the System. The \a tid parameter
    /// is used internally for sharing this structure among threads.
    State initial(unsigned tid);

    /// \brief Provides a string representation of the parameter state
    std::string to_string(const State, unsigned tid) const;

    /// \brief Returns an iterator over the successors of the parameter state.
    SuccIterator* succ(const State, unsigned tid);

    /// \brief Allocation and deallocation of iterator is costly. This
    /// method allows to reuse old iterators.
    void recycle(SuccIterator*, unsigned tid);

    /// \brief This method allow to deallocate a given state.
    const std::vector<std::string> get_ap();
  };

  /// \brief This method allows to ensure (at compile time) if
  /// a given parameter is of type kripkecube. It also check
  /// if the iterator has the good interface.
  template <typename State, typename SuccIter>
  bool is_a_kripkecube(kripkecube<State, SuccIter>&)
  {
    // Hardly waiting C++17 concepts...
    State (kripkecube<State, SuccIter>::*test_initial)(unsigned) =
      &kripkecube<State, SuccIter>::initial;
    std::string (kripkecube<State, SuccIter>::*test_to_string)
      (const State, unsigned) const = &kripkecube<State, SuccIter>::to_string;
    void (kripkecube<State, SuccIter>::*test_recycle)(SuccIter*, unsigned) =
      &kripkecube<State, SuccIter>::recycle;
    const std::vector<std::string>
      (kripkecube<State, SuccIter>::*test_get_ap)() =
      &kripkecube<State, SuccIter>::get_ap_str;
    void (SuccIter::*test_next)() = &SuccIter::next;
    State (SuccIter::*test_state)() const= &SuccIter::state;
    bool (SuccIter::*test_done)() const= &SuccIter::done;
    cube (SuccIter::*test_condition)() const = &SuccIter::condition;
    void (SuccIter::*test_fireall)() = &SuccIter::fireall;
    bool (SuccIter::*test_naturally_expanded)() const =
      &SuccIter::naturally_expanded;

    // suppress warnings about unused variables
    (void) test_initial;
    (void) test_to_string;
    (void) test_recycle;
    (void) test_get_ap;
    (void) test_next;
    (void) test_state;
    (void) test_done;
    (void) test_condition;
    (void) test_fireall;
    (void) test_naturally_expanded;

    // Always return true since otherwise a compile-time error will be raised.
    return true;
  }

  /// \ingroup kripke
  /// \brief Iterator code for Kripke structure
  ///
  /// This iterator can be used to simplify the writing
  /// of an iterator on a Kripke structure (or lookalike).
  ///
  /// If you inherit from this iterator, you should only
  /// redefine
  ///
  ///   - kripke_succ_iterator::first()
  ///   - kripke_succ_iterator::next()
  ///   - kripke_succ_iterator::done()
  ///   - kripke_succ_iterator::dst()
  ///
  /// This class implements kripke_succ_iterator::cond(),
  /// and kripke_succ_iterator::acc().
  class SPOT_API kripke_succ_iterator : public twa_succ_iterator
  {
  public:
    /// \brief Constructor
    ///
    /// The \a cond argument will be the one returned
    /// by kripke_succ_iterator::cond().
    kripke_succ_iterator(const bdd& cond)
      : cond_(cond)
    {
    }

    void recycle(const bdd& cond)
    {
      cond_ = cond;
    }

    virtual ~kripke_succ_iterator();

    virtual bdd cond() const override;
    virtual acc_cond::mark_t acc() const override;
  protected:
    bdd cond_;
  };

  /// \ingroup kripke
  /// \brief Interface for a Kripke structure
  ///
  /// A Kripke structure is a graph in which each node (=state) is
  /// labeled by a conjunction of atomic proposition.
  ///
  /// Such a structure can be seen as spot::tgba without
  /// any acceptance condition.
  ///
  /// A programmer that develops an instance of Kripke structure needs
  /// just provide an implementation for the following methods:
  ///
  ///   - kripke::get_init_state()
  ///   - kripke::succ_iter()
  ///   - kripke::state_condition()
  ///   - kripke::format_state()
  ///
  /// The other methods of the tgba interface (like those dealing with
  /// acceptance conditions) are supplied by this kripke class and
  /// need not be defined.
  ///
  /// See also spot::kripke_succ_iterator.
  class SPOT_API kripke: public fair_kripke
  {
  public:
    kripke(const bdd_dict_ptr& d)
      : fair_kripke(d)
      {
      }

    virtual ~kripke();

    virtual
      acc_cond::mark_t state_acceptance_mark(const state*) const override;
  };

  typedef std::shared_ptr<kripke> kripke_ptr;
  typedef std::shared_ptr<const kripke> const_kripke_ptr;
}
