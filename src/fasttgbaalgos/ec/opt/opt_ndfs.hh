// Copyright (C) 2012 Laboratoire de Recherche et Développement
// de l'Epita (LRDE).
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

#ifndef SPOT_FASTTGBAALGOS_EC_OPT_OPT_NDFS_HH
# define SPOT_FASTTGBAALGOS_EC_OPT_OPT_NDFS_HH

#include <tuple>

#include <stack>
#include <map>

#include "fasttgbaalgos/ec/ec.hh"
#include "fasttgbaalgos/ec/deadstore.hh"
#include "fasttgbaalgos/ec/lowlink_stack.hh"
#include "fasttgbaalgos/ec/root_stack.hh"

namespace spot
{

  /// \brief Compute the SCCs of a TGBA.
  ///
  /// This implementation use many optimisations
  ///   - Keep Live stack as small as possible [Nuutila]
  ///   - Use (or not!) a compressed lowlink stack
  class SPOT_API opt_ndfs : public ec
  {
  protected:
    enum rcol  {Blue, Cyan, Red};

    /// The map of visited states
    typedef Sgi::hash_map<const fasttgba_state*, rcol,
			  fasttgba_state_ptr_hash,
			  fasttgba_state_ptr_equal> seen_map;

  public:

    /// A constructor taking the automaton to check
    opt_ndfs(instanciator* i, std::string option = "",
		   bool swarm = false, int tid = 0);

    /// A destructor
    virtual ~opt_ndfs();

    /// The implementation of the interface
    bool check();

    /// Supply more information in a CSV way
    /// Informations are : Number of merge, number of states mark as dead.
    std::string extra_info_csv();

  protected:

    /// \brief Fix set ups for the algo
    void init();

    /// \brief Push a new state to explore
    virtual void dfs_push(fasttgba_state*);

    /// \brief  Pop states already explored
    virtual void dfs_pop();

    /// \brief the update for backedges
    virtual bool dfs_update (fasttgba_state* s);

    /// \brief the main procedure
    virtual void main();

    virtual void nested_dfs(const spot::fasttgba_state* a);

    /// \brief The color for a new State
    enum color {Alive, Dead, Unknown};

    // An element in Todo stack
    struct pair_state_iter
    {
      const spot::fasttgba_state* state;
      fasttgba_succ_iterator* lasttr;
      bool allred;
    };

    /// \brief the todo stack
    std::vector<pair_state_iter> todo;

    /// \brief Use a "generic" lowlink stack
    // stack_of_lowlink* dstack_;
    // generic_stack* stack_;

    /// \brief Access  the color of a state
    virtual opt_ndfs::color get_color(const fasttgba_state*);

    ///< \brief Storage for counterexample found or not
    bool counterexample_found;

    /// \brief the automaton that will be used for the Emptiness check
    const fasttgba* a_;

    // /// \brief The live stack
    // std::vector<const spot::fasttgba_state*> live;

    /// \brief the HashMap Live
    seen_map H;

    /// \brief The store for Deads states
    deadstore* deadstore_;

    /// \brief The instance automaton
    const instance_automaton* inst;

    unsigned int dfs_size_;	 ///< \brief keep dfs size
    unsigned int max_live_size_; ///< \brief keep peack size
    unsigned int max_dfs_size_;	 ///< \brief keep peack size
    int update_cpt_;		 ///< \brief count UPDATE calls
    int roots_poped_cpt_;	 ///< \brief count UPDATE loop iterations
    int states_cpt_;		 ///< \brief count states
    int transitions_cpt_;	 ///< \brief count transitions
    int memory_cost_;		 ///< \brief evaluates memory
    int trivial_scc_;            ///< \brief count trivial SCCs
    bool swarm_;		 ///< \brief shall use swarming?
    int tid_;                    ///< \brief the thread id
  };
}

#endif // SPOT_FASTTGBAALGOS_EC_OPT_OPT_NDFS_HH