// Copyright (C) 2012 Laboratoire de Recherche et Developpement de
// l'Epita (LRDE).
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


#ifndef SPOT_FASTTGBA_MARKSET_HH
# define SPOT_FASTTGBA_MARKSET_HH

#include <vector>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include "acc_dict.hh"

namespace spot
{

  /// \brief the type that represents a mark
  typedef int mark;

  /// \brief This class represents a set of acceptance marks
  class markset
  {
  private:
    /// \brief Internal Constructor
    markset(boost::dynamic_bitset<>, acc_dict& acc);

  public:
    /// \brief Initialize a mark of size \a size
    markset(acc_dict& acc);

    /// \brief destructor for a markset
    virtual ~markset(){ }

    /// \brief A copy constructor
    markset(const markset& b);

    /// \brief perform a logic AND with \a b
    ///
    /// Assign the result to this
    markset& operator&=(const markset& b);

    /// \brief a comparator for two markset
    bool operator==(const markset& b);

    /// \brief perform a logic OR with \a b
    ///
    /// Assign the result to this
    markset& operator|=(const markset& b);

    /// \brief perform a logic XOR with \a b
    ///
    /// Assign the result to this
    markset& operator^=(const markset & b);

    /// \brief Compute the difference with \a b
    ///
    /// Assign the result to this
    markset& operator-=(const markset& b);

    /// \brief Compute the difference with \a b
    ///
    /// Assign the result to this
    markset& operator-=(const mark b);

    /// \brief Retrun true if the mark is set
    bool is_set(const mark b);

    /// \brief Affect the value of  \a b in this
    markset& operator=(const markset& b);

    /// \brief Perform a logical AND
    markset operator&(const markset& b) const;

    /// \brief Perform a logical OR
    markset operator|(const markset& b) const;

    /// \brief Set a mark in the markset
    ///
    /// \param m the mark to be set
    /// Return true if the mark was not already set, false otherwise
    void set_mark(mark m);

    /// \brief test is a least one mark is set
    ///
    /// Return a boolean that indicates if there is
    /// at least a mark which is set
    bool empty();

    /// \brief test if all mark are set
    bool all() const;

    /// \brief Access to the first mark on the set
    ///
    /// To grab all mark, a loop which remove all
    /// mark must be realized.
    mark one();

    /// \brief Return the number of marks in the set
    size_t size();

    /// \brief Perform the negation of the mark
    ///
    /// Return a new mark which is the negation of thiss
    markset operator~() const;

    unsigned long to_ulong() const;
    markset operator|(const unsigned long) const;
    markset operator&(const unsigned long) const;

    /// \brief Display the content of the marking
    virtual std::string dump() const;

  protected:
    boost::dynamic_bitset<> markset_;   ///< the set of acceptance
    acc_dict& accs_;			///< the reference over the acc

  };
}

#endif  // SPOT_FASTTGBA_MARKSET_HH