// Copyright (C) 2013 Laboratoire de Recherche et Développement
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

// #define UNION_SCCTRACE
#ifdef UNION_SCCTRACE
#define trace std::cerr
#else
#define trace while (0) std::cerr
#endif

#include "fasttgba/fasttgba_product.hh"
#include <iostream>
#include "union_scc.hh"
#include <assert.h>

namespace spot
{
  union_scc::union_scc(instanciator* i, std::string option) :
    counterexample_found(false),
    inst(i->new_instance()),
    max_dfs_size_(0),
    update_cpt_(0),
    update_loop_cpt_(0),
    roots_poped_cpt_(0),
    states_cpt_(0),
    transitions_cpt_(0),
    memory_cost_(0),
    trivial_scc_(0)
  {
    a_ = inst->get_automaton ();

    if (!option.compare("-cs-ds"))
      {
	K = 3;
	uf  = new setOfDisjointSetsIPC_LRPC_MS (a_->get_acc());
	roots_stack_ = new stack_of_roots (a_->get_acc());
      }
    else if (!option.compare("-cs+ds"))
      {
	K = 4;
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	roots_stack_ = new stack_of_roots (a_->get_acc());
      }
    else if (!option.compare("+cs-ds"))
      {
	K = 3;
	uf  = new setOfDisjointSetsIPC_LRPC_MS (a_->get_acc());
	roots_stack_ = new compressed_stack_of_roots (a_->get_acc());
      }
    else
      {
	K = 4;
	assert(!option.compare("+cs+ds") || !option.compare(""));
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	roots_stack_ = new compressed_stack_of_roots (a_->get_acc());
      }
    empty_ = new markset(a_->get_acc());
  }

  union_scc::~union_scc()
  {
    delete empty_;
    delete roots_stack_;
    delete uf;
    while (!todo.empty())
      {
    	delete todo.back().lasttr;
    	todo.pop_back();
      }

    delete inst;
  }

  bool
  union_scc::check()
  {
    init();
    main();
    return counterexample_found;
  }

  void union_scc::init()
  {
    trace << "Union_Scc::Init" << std::endl;
    fasttgba_state* init = a_->get_init_state();
    dfs_push(init);
  }

  void union_scc::dfs_push(fasttgba_state* s)
  {
    trace << "Union_Scc::DFS_push "
    	  << std::endl;
    ++states_cpt_;

    int p;
    uf->add (s, &p);

    todo.push_back ({s, 0});
    // Count !
    max_dfs_size_ = max_dfs_size_ > todo.size() ?
      max_dfs_size_ : todo.size();
    roots_stack_->push_trivial(todo.size()-1);

    int tmp_cost = 1*roots_stack_->size() + K*uf->size()
      + 1*uf->dead_size();
    if (tmp_cost > memory_cost_)
      memory_cost_ = tmp_cost;

  }

  void union_scc::dfs_pop()
  {
    trace << "Union_Scc::DFS_pop " << std::endl;
    pair_state_iter pair = todo.back();
    delete pair.lasttr;
    todo.pop_back();

    if (todo.size() == roots_stack_->root_of_the_top())
      {
	++roots_poped_cpt_;
	uf->make_dead(pair.state);
	roots_stack_->pop();
      }
  }

  bool union_scc::merge(fasttgba_state* d)
  {
    trace << "Union_Scc::merge " << std::endl;
    ++update_cpt_;
    int i = roots_stack_->root_of_the_top();
    assert(todo[i].lasttr);
    roots_stack_->pop();

    while (!uf->same_partition(d, todo[i].state))
      {
	++update_loop_cpt_;
 	uf->unite(d, todo[i].state);
	assert(todo[i].lasttr);
	i = roots_stack_->root_of_the_top();
	roots_stack_->pop();
      }
    roots_stack_->push_non_trivial(i, *empty_, todo.size() -1);

    return false;
  }

  void union_scc::main()
  {
    union_find::color c;
    while (!todo.empty())
      {
	trace << "Main " << std::endl;
	assert(!uf->is_dead(todo.back().state));

	if (!todo.back().lasttr)
	  {
	    todo.back().lasttr = a_->succ_iter(todo.back().state);
	    todo.back().lasttr->first();
	  }
	else
	  {
	    assert(todo.back().lasttr);
	    todo.back().lasttr->next();
	  }

    	if (todo.back().lasttr->done())
    	  {
	    dfs_pop ();
    	  }
    	else
    	  {
	    ++transitions_cpt_;
	    assert(todo.back().lasttr);
    	    fasttgba_state* d = todo.back().lasttr->current_state();
	    c = uf->get_color(d);
    	    if (c == union_find::Unknown)
    	      {
		dfs_push (d);
    	    	continue;
    	      }
    	    else if (c == union_find::Alive)
    	      {
    	    	if (merge (d))
    	    	  {
    	    	    counterexample_found = true;
    	    	    d->destroy();
    	    	    return;
    	    	  }
    	      }
    	    d->destroy();
    	  }
      }
  }

  std::string
  union_scc::extra_info_csv()
  {
    // dfs max size
    // root max size
    // live max size
    // deadstore max size
    // number of UPDATE calls
    // number of loop inside UPDATE
    // Number of Roots poped
    // visited states
    // visited transitions

    return
      std::to_string(max_dfs_size_)
      + ","
      + std::to_string(roots_stack_->max_size())
      + ","
      + std::to_string(uf->max_alive())
      + ","
      + std::to_string(uf->max_dead())
      + ","
      + std::to_string(update_cpt_)
      + ","
      + std::to_string(update_loop_cpt_)
      + ","
      + std::to_string(roots_poped_cpt_)
      + ","
      + std::to_string(transitions_cpt_)
      + ","
      + std::to_string(states_cpt_)
      + ","
      + std::to_string(memory_cost_)
      + ","
      + std::to_string(trivial_scc_);
  }
}