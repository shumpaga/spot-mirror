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

// #define UNIONCHECKTRACE
#ifdef UNIONCHECKTRACE
#define trace std::cerr
#else
#define trace while (0) std::cerr
#endif

#include "fasttgba/fasttgba_product.hh"
#include <iostream>
#include "unioncheck.hh"
#include <assert.h>

namespace spot
{
  unioncheck::unioncheck(instanciator* i, std::string option) :
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

    if (!option.compare("-cs"))
      {
	K = 4;
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	stack_ = new generic_stack (a_->get_acc());
      }
    else
      {
	K = 4;
	assert(!option.compare("+cs") || !option.compare(""));
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	stack_ = new compressed_generic_stack (a_->get_acc());
      }
  }

  unioncheck::~unioncheck()
  {
    delete stack_;
    delete uf;
    while (!todo.empty())
      {
    	delete todo.back().lasttr;
    	todo.pop_back();
      }

    delete inst;
  }

  bool
  unioncheck::check()
  {
    init();
    main();
    return counterexample_found;
  }

  void unioncheck::init()
  {
    trace << "Unioncheck::Init" << std::endl;
    fasttgba_state* init = a_->get_init_state();
    dfs_push(init);
  }

  void unioncheck::dfs_push(fasttgba_state* s)
  {
    trace << "Unioncheck::DFS_push "
	  << a_->format_state(s)
    	  << std::endl;
    ++states_cpt_;

    int p;
    uf->add (s, &p);
    stack_->push_transient(todo.size());

    todo.push_back ({s, 0});
    // Count !
    max_dfs_size_ = max_dfs_size_ > todo.size() ?
      max_dfs_size_ : todo.size();

    int tmp_cost = 2*stack_->size() + K*uf->size()
      + 1*uf->dead_size();
    if (tmp_cost > memory_cost_)
      memory_cost_ = tmp_cost;

  }

  void unioncheck::dfs_pop()
  {
    trace << "Unioncheck::DFS_pop " << std::endl;
    pair_state_iter pair = todo.back();
    delete pair.lasttr;
    todo.pop_back();

    // todo is already popped (no need -1)
    if (todo.size() == stack_->top(todo.size()).pos)
      {
	++roots_poped_cpt_;
	uf->make_dead(pair.state);
	stack_->pop(todo.size());
      }
  }

  bool unioncheck::merge(fasttgba_state* d)
  {
    trace << "Unioncheck::merge "
	  << a_->format_state(todo.back().state)
	  << std::endl;
    ++update_cpt_;

    auto top = stack_->pop(todo.size()-1);
    //    markset a //= top.acc |
    top.acc |= todo.back().lasttr->current_acceptance_marks();
    int r = top.pos;
    assert(todo[r].state);

    while (!uf->same_partition(d, todo[r].state))
      {
	++update_loop_cpt_;
 	uf->unite(d, todo[r].state);
	assert(todo[r].lasttr);
	auto newtop = stack_->pop(r-1);

	// [r-1] Because acceptances are stored in the predecessor!
	top.acc |= newtop.acc | todo[r-1].lasttr->current_acceptance_marks();
	r = newtop.pos;
      }
    stack_->push_non_transient(r, top.acc);

    return top.acc.all();
  }

  void unioncheck::main()
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
  unioncheck::extra_info_csv()
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
      + std::to_string(stack_->max_size())
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


  tarjanunioncheck::tarjanunioncheck(instanciator* i, std::string option) :
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

    if (!option.compare("-cs"))
      {
	K = 4;
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	stack_ = new generic_stack(a_->get_acc());
      }
    else
      {
	K = 4;
	assert(!option.compare("+cs") || !option.compare(""));
	uf  = new setOfDisjointSetsIPC_LRPC_MS_Dead (a_->get_acc());
	stack_ = new compressed_generic_stack (a_->get_acc());
      }
  }

  tarjanunioncheck::~tarjanunioncheck()
  {
    delete stack_;
    delete uf;
    while (!todo.empty())
      {
    	delete todo.back().lasttr;
    	todo.pop_back();
      }

    delete inst;
  }

  bool
  tarjanunioncheck::check()
  {
    init();
    main();
    return counterexample_found;
  }

  void tarjanunioncheck::init()
  {
    trace << "Tarjanunioncheck::Init" << std::endl;
    fasttgba_state* init = a_->get_init_state();
    dfs_push(init);
  }

  void tarjanunioncheck::dfs_push(fasttgba_state* s)
  {
    trace << "Tarjanunioncheck::DFS_push "
	  << a_->format_state(s)
    	  << std::endl;
    ++states_cpt_;

    int p = 0;
    uf->add (s, &p);
    todo.push_back ({s, 0, p});

    // Count !
    max_dfs_size_ = max_dfs_size_ > todo.size() ?
      max_dfs_size_ : todo.size();
    stack_->push_transient(p);

    int tmp_cost = 2*stack_->size() + K*uf->size() + 1*uf->dead_size();
    if (tmp_cost > memory_cost_)
      memory_cost_ = tmp_cost;

  }

  void tarjanunioncheck::dfs_pop()
  {
    trace << "Tarjanunioncheck::DFS_pop "
	  << a_->format_state(todo.back().state)
	  << std::endl;
    dfs_element pair = todo.back();
    delete pair.lasttr;
    todo.pop_back();

    auto top = stack_->pop(pair.pos);

    if (pair.pos == (int) top.pos)
      {
	++roots_poped_cpt_;
	uf->make_dead(pair.state);
      }
    else
      {
 	uf->unite(pair.state, todo.back().state);

	auto newtop = stack_->pop(todo.back().pos);
	newtop.acc |= top.acc |
	  todo.back().lasttr->current_acceptance_marks();

	if (top.pos <= newtop.pos)
	  stack_->push_non_transient(top.pos, newtop.acc);
	else
	  stack_->push_non_transient(newtop.pos, newtop.acc);

	if (newtop.acc.all())
	  counterexample_found = true;
      }
  }

  bool tarjanunioncheck::dfs_update(fasttgba_state* d)
  {
    trace << "Tarjanunioncheck::merge "
	  << a_->format_state(todo.back().state)
	  << std::endl;
    ++update_cpt_;
    uf->unite(d, todo.back().state);

    auto top = stack_->pop(todo.back().pos);
    top.acc |= todo.back().lasttr->current_acceptance_marks();

    int livenum = uf->live_get(d);
    if (livenum <= (int) top.pos)
      stack_->push_non_transient(livenum, top.acc);
    else
      stack_->push_non_transient(top.pos, top.acc);

    return top.acc.all();
  }

  void tarjanunioncheck::main()
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
	    if (counterexample_found)
	      return;
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
    	    	if (dfs_update(d))
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
  tarjanunioncheck::extra_info_csv()
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
      + std::to_string(stack_->max_size())
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