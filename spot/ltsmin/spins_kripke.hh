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

#include <bricks/brick-hash>
#include <bricks/brick-hashset>
#include <spot/kripke/kripke.hh>
#include <spot/ltsmin/spins_interface.hh>
#include <spot/misc/fixpool.hh>
#include <spot/misc/mspool.hh>
#include <spot/misc/intvcomp.hh>
#include <spot/misc/intvcmp2.hh>
#include <spot/twacube/cube.hh>

/// This file aggregates all classes and typedefs necessary
/// to build a kripke that is thread safe
namespace spot
{

  /// \brief A Spins state is represented as an array of integer
  /// Note that this array has two reserved slots (position 0 an 1).
  ///
  /// At position 0 we store the hash associated to the state to avoid
  /// multiple computations.
  ///
  /// At position 1 we store the size of the state: keeping this information
  /// allows to compress the state
  typedef int* cspins_state;

  /// \brief This class provides the ability to compare two states
  struct cspins_state_equal
  {
    bool operator()(const cspins_state lhs, const cspins_state rhs) const
    {
      return 0 == memcmp(lhs, rhs, (2+rhs[1])* sizeof(int));
    }
  };

  /// \brief This class provides the ability to hash a state
  struct cspins_state_hash
  {
    size_t operator()(const cspins_state that) const
    {
      return that[0];
    }
  };

  /// \brief This class provides a hasher as required by the bricks classes
  struct cspins_state_hasher
  {
    cspins_state_hasher(cspins_state&) { }
    cspins_state_hasher() = default;
    brick::hash::hash128_t hash(cspins_state t) const
    {
      // FIXME we should compute a better hash value for this particular
      // case. Shall we use two differents hash functions?
      return std::make_pair(t[0], t[0]);
    }
    bool equal(cspins_state lhs, cspins_state rhs) const
    {
      return 0 == memcmp(lhs, rhs, (2+rhs[1])* sizeof(int));
    }
  };

  /// \brief The management of states (i.e. allocation/deallocation) can
  /// be painless since every time we have to consider wether the state will
  /// be compressed or not. This class aims to simplify this management
  class cspins_state_manager final
  {
  public:
    cspins_state_manager(unsigned int state_size, int compress);
    int* unbox_state(cspins_state s) const;
    // cmp is the area we can use to compute the compressed rep of dst.
    cspins_state alloc_setup(int *dst, int* cmp, size_t cmpsize);
    void decompress(cspins_state s, int* uncompressed, unsigned size) const;
    void  dealloc(cspins_state s);
    unsigned int size() const;

  private:
    fixed_size_pool p_;
    multiple_size_pool msp_;
    bool compress_;
    const unsigned int state_size_;
    void (*fn_compress_)(const int*, size_t, int*, size_t&);
    void (*fn_decompress_)(const int*, size_t, int*, size_t);
  };


  // This structure is used as a parameter during callback when
  // generating states from the shared librarie produced by LTSmin
  struct inner_callback_parameters
  {
    cspins_state_manager* manager;   // The state manager
    std::vector<cspins_state>* succ; // The successors of a state
    int* compressed_;
    int* uncompressed_;
    bool compress;
    bool selfloopize;
  };

  // This class provides an iterator over the successors of a state.
  // All successors are computed once when an iterator is recycled or
  // created.
  class cspins_iterator final
  {
  public:
    cspins_iterator(cspins_state s,
                    const spot::spins_interface* d,
                    cspins_state_manager& manager,
                    inner_callback_parameters& inner,
                    cube cond,
                    bool compress,
                    bool selfloopize,
                    cubeset& cubeset,
                    int dead_idx, unsigned tid);

    void recycle(cspins_state s,
                 const spot::spins_interface* d,
                 cspins_state_manager& manager,
                 inner_callback_parameters& inner,
                 cube cond,
                 bool compress,
                 bool selfloopize,
                 cubeset& cubeset, int dead_idx, unsigned tid);

    ~cspins_iterator();
    void next();
    bool done() const;
    cspins_state state() const;
    cube condition() const;


  private:
    std::vector<cspins_state> successors_;
    unsigned int current_;
    cube cond_;
    unsigned tid_;
  };


  // A specialisation of the template class kripke that is thread safe.
  template<>
  class kripkecube<cspins_state, cspins_iterator> final
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
      int lval;
      relop op;
      int rval;
    };

    // Data structure to store complex atomic propositions
    typedef std::vector<one_prop> prop_set;
    prop_set pset_;

  public:
    kripkecube(spins_interface_ptr sip, bool compress,
               std::vector<std::string> visible_aps,
               bool selfloopize, std::string dead_prop,
               unsigned int nb_threads);
    ~kripkecube();
    cspins_state initial(unsigned tid);
    std::string to_string(const cspins_state s, unsigned tid = 0) const;
    cspins_iterator* succ(const cspins_state s, unsigned tid);
    void recycle(cspins_iterator* it, unsigned tid);
    const std::vector<std::string> get_ap();
    unsigned get_threads();

  private:
    /// Parse the set of atomic proposition to have a more
    /// efficient data strucure for computation
    void match_aps(std::vector<std::string>& aps, std::string dead_prop);

    /// Compute the cube associated to each state. The cube
    /// will then be given to all iterators.
    void compute_condition(cube c, cspins_state s, unsigned tid = 0);

    spins_interface_ptr sip_;
    const spot::spins_interface* d_; // To avoid numerous sip_.get()
    cspins_state_manager* manager_;
    bool compress_;
    std::vector<std::vector<cspins_iterator*>> recycle_;
    inner_callback_parameters* inner_;
    cubeset cubeset_;
    bool selfloopize_;
    int dead_idx_;
    std::vector<std::string> aps_;
    unsigned int nb_threads_;
  };

  /// \brief shortcut to manipulate the kripke below
  typedef std::shared_ptr<spot::kripkecube<spot::cspins_state,
                                           spot::cspins_iterator>>
  ltsmin_kripkecube_ptr;

}

#include <spot/ltsmin/spins_kripke.hxx>
