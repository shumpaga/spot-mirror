// -*- coding: utf-8 -*-
// Copyright (C)  2015, 2016, 2017, 2018 Laboratoire de Recherche et
// Developpement de l'Epita (LRDE)
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
#include <vector>
#include <cstring>
#include <unordered_set>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <spot/ltsmin/porinfos.hh>

namespace spot
{
  porinfos_stats::porinfos_stats()
  { }

  void porinfos_stats::cumul(unsigned red, unsigned en)
  {
    red_ += red;
    en_ += en;
    ++cumul_;
    min_red_ = min_red_ < red ? min_red_ : red;
    max_red_ = max_red_ > red ? max_red_ : red;
    min_en_ = min_en_ < en ? min_en_ : en;
    max_en_ = max_en_ > en ? max_en_ : en;
  }

  void porinfos_stats::dump()
  {
    assert(cumul_);
    std::cout << "reduced_set_min: " << min_red_ << std::endl;
    std::cout << "reduced_set_avg: " << ((int) (red_/cumul_)) << std::endl;
    std::cout << "reduced_set_max: " << max_red_ << std::endl;
    std::cout << "enabled_set_min: " << min_en_ << std::endl;
    std::cout << "enabled_set_avg: " << ((int) (en_/cumul_)) << std::endl;
    std::cout << "enabled_set_max: " << max_en_ << std::endl;
  }

  porinfos::porinfos(const spins_interface* si)
    : porinfos(si, porinfos_options())
  {
  }

  porinfos::porinfos(const spins_interface* si, const porinfos_options& opt)
  {
    d_ = si;
    transitions_ = si->get_transition_count();
    variables_ = si->get_state_size();
    guards_ = si->get_guard_count();
    m_read.resize(transitions_);
    m_write.resize(transitions_);
    m_nes.resize(guards_);
    m_mbc.resize(guards_);
    m_guards.resize(transitions_);
    m_guard_variables.resize(guards_);

    switch (opt.method)
      {
        case porinfos_options::porinfos_method::stubborn_set:
          f_not_enabled_transition = &porinfos::stubborn_set;
        break;
        case porinfos_options::porinfos_method::base:
          f_not_enabled_transition = &porinfos::base;
        break;
      }

    transparent_ =
      opt.ltl_condition ==
        porinfos_options::porinfos_ltl_condition::transparent;
    invisible_ =
      transparent_ ||
        opt.ltl_condition == porinfos_options::porinfos_ltl_condition::invisible;

    for (unsigned i = 0; i < transitions_; ++i)
      m_read[i].resize(variables_);
    for (unsigned i = 0; i < transitions_; ++i)
      m_write[i].resize(variables_);
    for (unsigned i = 0; i < guards_; ++i)
      m_nes[i].resize(transitions_);
    for (unsigned i = 0; i < guards_; ++i)
      m_mbc[i].resize(guards_);
    for (unsigned i = 0; i < guards_; ++i)
      m_guard_variables[i].resize(variables_);

    // Grab Read dependency Matrix
    for (unsigned i = 0; i < transitions_; ++i)
      {
        int* read_dep = d_->get_transition_read_dependencies(i);
        for (unsigned j = 0; j < variables_; ++j)
          m_read[i][j] = read_dep[j];
      }

    // Grab Write depency Matrix
    for (unsigned i = 0; i < transitions_; ++i)
      {
        int* write_dep = d_->get_transition_write_dependencies(i);
        for (unsigned j = 0; j < variables_; ++j)
          m_write[i][j] = write_dep[j];
      }

    // Grab NES matrix
    for (unsigned i = 0; i < guards_; ++i)
      {
        int* nes = d_->get_guard_nes_matrix(i);
        for (unsigned j = 0; j < transitions_; ++j)
          m_nes[i][j] = nes[j];
      }

    // Grab "may-be-coenabled" matrix
    for (unsigned i = 0; i < guards_; ++i)
      {
        int* mbc = d_->get_guard_may_be_coenabled_matrix(i);
        for (unsigned j = 0; j < guards_; ++j)
          m_mbc[i][j] = mbc[j];
      }

    // Grab guards matrix
    for (unsigned i = 0; i < transitions_; ++i)
      {
        int* array_guards = d_->get_guards(i);
        m_guards[i].resize(array_guards[0]);
        for (unsigned j = 0; j < (unsigned) array_guards[0]; ++j)
          m_guards[i][j] = array_guards[j+1];
      }

    // Grab guard variables matrix
    for (unsigned i = 0; i < guards_; ++i)
      {
        int* guard_variables = d_->get_guard_variables_matrix(i);
        for (unsigned j = 0; j < variables_; ++j)
          m_guard_variables[i][j] = guard_variables[j];
      }

    // Compute the list of processes with their id.
    // Primitives types and processes are stored in the same place, so
    // we just have to keep only entries which are not a type.
    const std::vector<std::string> primitive_types = { "int", "byte" };
    std::vector<std::pair<std::string, int>> processes;
    unsigned type_count = d_->get_type_count();
    for (unsigned i = 0; i < type_count; ++i)
      {
        auto type = std::string(d_->get_type_name(i));
        if (std::find(primitive_types.begin(), primitive_types.end(), type)
           == primitive_types.end())
          processes.emplace_back(type, i);
      }

    // Compute the list of PC (Process Counter) and their processes.
    // -1 means that the variable is not a PC.
    std::vector<int> pc_proc(d_->get_state_size());
    std::vector<unsigned> pc_numbers(0); // Contains the number of the variable
    std::vector<unsigned> pc_type(0); // Contains the number of the type

    for (unsigned i = 0; i < variables_; ++i)
      {
        auto name = std::string(d_->get_state_variable_name(i));
        auto pred = [&name](std::pair<std::string, int> p)
          { return p.first == name; };
        auto res = std::find_if(processes.begin(), processes.end(), pred);

        if (res != processes.end())
        {
          pc_proc[i] = res->second;
          pc_numbers.push_back(i);
          pc_type.push_back(res->second);
        }
        else
          pc_proc[i] = -1;
      }

    // Compute all processes of each transition, based on the tested
    // PC. A transition can be synchronized with several
    // processes, so a transition may have more than one process.
    t_processes.resize(transitions_);
    for (unsigned t = 0; t < transitions_; ++t)
      {
        t_processes[t].clear();
        for (unsigned g = 0; g < m_guards[t].size(); g++)
        {
          for (auto pc : pc_numbers)
            if (m_guard_variables[m_guards[t][g]][pc])
              t_processes[t].insert(pc_proc[pc]);
        }
       /* for (unsigned v = 0; v < variables_; ++v)
          {
            if (m_read[t][v])
              {
                if (pc_proc[v] != -1)
                  t_processes[t].insert(pc_proc[v]);
              }
          } */
      }

    // Compute the matrix of state dependancy (true iff two transitions start
    // from the same local state in a process).
    m_dep_state.resize(transitions_);
    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      m_dep_state[t1].resize(transitions_);

    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      for (unsigned t2 = t1; t2 < transitions_; ++t2)
        {
          bool res = false;
          for (unsigned i = 0; i < m_guards[t1].size() && !res; i++)
            for (unsigned j = 0; j < m_guards[t2].size() && !res; j++)
              if (m_guards[t1][i] == m_guards[t2][j])
                for (auto pc : pc_numbers)
                  if (m_guard_variables[m_guards[t1][i]][pc])
                    res = true;

          m_dep_state[t1][t2] = res;
          m_dep_state[t2][t1] = res;
        }

    // Compute the matrix of processes dependency (true iff two transitions
    // belong to one same process).
    m_dep_process.resize(transitions_);
    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      m_dep_process[t1].resize(transitions_);

    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      {
        for (unsigned t2 = t1; t2 < transitions_; ++t2)
          {
            std::vector<unsigned> res(t_processes[t1].size());

            auto it = std::set_intersection(t_processes[t1].begin(),
                                            t_processes[t1].end(),
                                            t_processes[t2].begin(),
                                            t_processes[t2].end(), res.begin());
            res.resize(it - res.begin());

            m_dep_process[t1][t2] = !res.empty();
            m_dep_process[t2][t1] = !res.empty();
          }
      }

    // Setup dependency between transitions (true iff two transitions read
    // or write one same variable).
    m_dep_tr.resize(transitions_);
    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      m_dep_tr[t1].resize(transitions_);

    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      {
        for (unsigned t2 = t1; t2 < transitions_; ++t2)
          {
            bool res = false;
            for (unsigned i = 0; i < variables_; ++i)
              if ((m_read[t1][i] || m_write[t1][i])
                && (m_read[t2][i] || m_write[t2][i]))
                {
                  res = true;
                  break;
                }
            m_dep_tr[t1][t2] = res;
            m_dep_tr[t2][t1] = res;
          }
      }

    // Setup conflicting transitions matrix
    m_conflict_tr.resize(transitions_);
    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      m_conflict_tr[t1].resize(transitions_);

    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      {
        for (unsigned t2 = t1; t2 < transitions_; ++t2)
          {
            bool res = (m_dep_state[t1][t2])
                       || (!m_dep_process[t1][t2] && m_dep_tr[t1][t2]);
            m_conflict_tr[t1][t2] = res;
            m_conflict_tr[t2][t1] = res;
          }
      }

    // Setup non maybe coenabled
    non_mbc_tr.resize(transitions_);
    for (unsigned i = 0; i < transitions_; ++i)
      non_mbc_tr[i].resize(transitions_);

    for (unsigned t1 = 0; t1 < transitions_; ++t1)
      for (unsigned t2 = t1; t2 < transitions_; ++t2)
        {
          non_mbc_tr[t1][t2] = non_maybecoenabled(t1, t2);
          non_mbc_tr[t2][t1] = non_mbc_tr[t1][t2];
        }

    // Speedup Stubborn set computation
    std::vector<int> tmp;
    for (unsigned i = 0; i < guards_; ++i)
      {
        unsigned guard_to_look = i;
        for (unsigned j = 0; j < transitions_; ++j)
          {
            if (m_nes[guard_to_look][j])
              tmp.emplace_back(j);
          }
        gnes_cache_.insert({guard_to_look, tmp});
        tmp.clear();
      }


    for (unsigned beta = 0; beta < transitions_; ++beta)
      {
        for (unsigned i = 0; i < transitions_; ++i)
          {
            // Beta is already in t_s or transitions are independants ...
            // ... so we avoid extra computation.
            if (i == (unsigned) beta || !m_dep_tr[i][beta])
              continue;

            // The transition must be processed!
            if (non_mbc_tr[beta][i])
              {
                tmp.emplace_back(i);
              }
          }
        gnes_cache_e_.insert({beta, tmp});
        tmp.clear();
      }

    if (invisible_)
      {
        if (transparent_)
          {
            // Get the patterns of transparent transitions from Divine.
            // This needs the function get_effect in Divine.
            t_simple_effects.resize(transitions_);
            for (unsigned t = 0; t < transitions_; ++t)
              {
                get_effect_data data = { t, &t_simple_effects, d_, variables_ };
                d_->get_effects(t, [] (const char* var_name, const char* str, void* arg) -> void
                  {
                    get_effect_data* data = reinterpret_cast<get_effect_data*>(arg);
                    std::vector<std::vector<simple_op>>* t_simple_effects =
                      data->t_simple_effects;
                    const spins_interface* d = data->d;
                    unsigned variables_ = data->variables;
                    unsigned t = data->t;
                    unsigned var_num = 0;
                    for (unsigned i = 0; i < variables_; ++i)
                      if (strcmp(var_name, d->get_state_variable_name(i)) == 0)
                        {
                          var_num = i;
                          break;
                        }
                    if (strcmp(str, "incr") == 0)
                      (*t_simple_effects)[t].emplace_back(simple_op{ simple_ops_types::incr, var_num });
                    if (strcmp(str, "decr") == 0)
                      (*t_simple_effects)[t].emplace_back(simple_op{ simple_ops_types::decr, var_num });
                  }, (void*) &data);
              }
          }

        // Compute the list of invisible/transparent transitions.
        m_invisible.resize(transitions_);
        for (unsigned t = 0; t < transitions_; ++t)
          {
            m_invisible[t] = true;
            for (const auto& op: *(opt.aps))
              {
                bool transparent_ok_r = false;
                bool transparent_ok_l = false;
                if (transparent_)
                  {
                    for (auto simple_eff : t_simple_effects[t])
                      {
                        if (simple_eff.op_type == simple_ops_types::decr)
                          {
                            if ((op.op == kripke_cube::VAR_OP_LT ||
                                 op.op == kripke_cube::VAR_OP_LE) &&
                                 op.lval.var == simple_eff.variable && op.pos)
                              transparent_ok_l = true;
                            else if ((op.op == kripke_cube::OP_LT_VAR ||
                                 op.op == kripke_cube::OP_LE_VAR) &&
                                 op.rval.var == simple_eff.variable && !op.pos)
                              transparent_ok_r = true;
                            else if ((op.op == kripke_cube::VAR_OP_GT ||
                                 op.op == kripke_cube::VAR_OP_GE) &&
                                 op.lval.var == simple_eff.variable && !op.pos)
                              transparent_ok_l = true;
                            else if ((op.op == kripke_cube::OP_GT_VAR ||
                                 op.op == kripke_cube::OP_GE_VAR) &&
                                 op.rval.var == simple_eff.variable && op.pos)
                              transparent_ok_r = true;
                          }
                        else if (simple_eff.op_type == simple_ops_types::incr)
                          {
                            if ((op.op == kripke_cube::VAR_OP_LT ||
                                 op.op == kripke_cube::VAR_OP_LE) &&
                                 op.lval.var == simple_eff.variable && !op.pos)
                              transparent_ok_l = true;
                            else if ((op.op == kripke_cube::OP_LT_VAR ||
                                 op.op == kripke_cube::OP_LE_VAR) &&
                                 op.rval.var == simple_eff.variable && op.pos)
                              transparent_ok_r = true;
                            else if ((op.op == kripke_cube::VAR_OP_GT ||
                                 op.op == kripke_cube::VAR_OP_GE) &&
                                 op.lval.var == simple_eff.variable && op.pos)
                              transparent_ok_l = true;
                            else if ((op.op == kripke_cube::OP_GT_VAR ||
                                 op.op == kripke_cube::OP_GE_VAR) &&
                                 op.rval.var == simple_eff.variable && !op.pos)
                              transparent_ok_r = true;
                          }
                      }
                  }

                if (!transparent_ok_l && (op.op == kripke_cube::VAR_OP_EQ
                    || op.op == kripke_cube::VAR_OP_NE
                    || op.op == kripke_cube::VAR_OP_LT
                    || op.op == kripke_cube::VAR_OP_GT
                    || op.op == kripke_cube::VAR_OP_LE
                    || op.op == kripke_cube::VAR_OP_GE
                    || op.op == kripke_cube::VAR_OP_EQ_VAR
                    || op.op == kripke_cube::VAR_OP_NE_VAR
                    || op.op == kripke_cube::VAR_OP_LT_VAR
                    || op.op == kripke_cube::VAR_OP_GT_VAR
                    || op.op == kripke_cube::VAR_OP_LE_VAR
                    || op.op == kripke_cube::VAR_OP_GE_VAR))
                  {
                    if (m_write[t][op.lval.var])
                      m_invisible[t] = false;
                  }
                if (!transparent_ok_r && (op.op == kripke_cube::OP_EQ_VAR
                    || op.op == kripke_cube::OP_NE_VAR
                    || op.op == kripke_cube::OP_LT_VAR
                    || op.op == kripke_cube::OP_GT_VAR
                    || op.op == kripke_cube::OP_LE_VAR
                    || op.op == kripke_cube::OP_GE_VAR
                    || op.op == kripke_cube::VAR_OP_EQ_VAR
                    || op.op == kripke_cube::VAR_OP_NE_VAR
                    || op.op == kripke_cube::VAR_OP_LT_VAR
                    || op.op == kripke_cube::VAR_OP_GT_VAR
                    || op.op == kripke_cube::VAR_OP_LE_VAR
                    || op.op == kripke_cube::VAR_OP_GE_VAR))
                  {
                    if (m_write[t][op.rval.var])
                      m_invisible[t] = false;
                  }
              }
          }
      }
  }

  bool
  porinfos::base(int t, std::vector<int>& t_work,
                 const std::vector<int>& t_s,
                 const std::vector<int>& enabled,
                 const int* for_spins_state)
  {
    (void) t;
    (void) t_work;
    (void) t_s;
    (void) enabled;
    (void) for_spins_state;
    return false;
  }

  bool
  porinfos::stubborn_set(int t, std::vector<int>& t_work,
                         const std::vector<int>& t_s,
                         const std::vector<int>& enabled,
                         const int* for_spins_state)
  {
    (void) enabled;
    int beta = t;
    unsigned beta_guards_size = m_guards[beta].size();
    bool goon = true;
    for (unsigned i = 0; i < beta_guards_size && goon; ++i)
      {
        unsigned guard_to_look = m_guards[beta][i];
        if (!(d_->get_guard(nullptr, guard_to_look, for_spins_state)))
          {
            goon = false;
            for (unsigned j = 0; j < transitions_; ++j)
              {
                if (m_nes[guard_to_look][j]
                   && std::find(t_work.begin(), t_work.end(), j) == t_work.end()
                   && std::find(t_s.begin(), t_s.end(), j) == t_s.end())
                  {
                    t_work.push_back(j);
                  }
              }
          }
      }
    return true;
  }

  std::vector<bool>
  porinfos::compute_reduced_set(const std::vector<int>& enabled,
                                const int* for_spins_state)
  {
    (void) for_spins_state;
    std::vector<bool> res(enabled.size(), true);

    if (enabled.empty())
      {
        stats_.cumul(0, enabled.size());
        return res;
      }

    // Compute the stubborn set algorithm as described by Elwin Pater
    // in Partial Order Reduction for PINS [2011] (page 21)

    // Declare usefull variables
    thread_local std::vector<int> t_work;
    thread_local std::vector<int> t_s;
    t_work.clear();
    t_s.clear();

    // Randomly take one enabled transition
    {
      int alpha = enabled[0];
      t_work.push_back(alpha);
    }

    // Compute the conflicting transitions algorithm as described by
    // Elwin Pater in Partial Order Reduction for PINS [2011] (page 15)

    while (!t_work.empty())
      {
        unsigned alpha = t_work.back();
        t_work.pop_back();
        t_s.push_back(alpha);

        if (std::find(enabled.begin(), enabled.end(), alpha) !=
            enabled.end())
          {
            for (unsigned beta = 0; beta < transitions_; ++beta)
              if (alpha != beta && m_conflict_tr[alpha][beta]
                && std::find(t_work.begin(), t_work.end(), beta) == t_work.end()
                && std::find(t_s.begin(), t_s.end(), beta) == t_s.end())
                t_work.push_back(beta);
          }
        else
          {
            if (!((this->*f_not_enabled_transition)(alpha, t_work, t_s, enabled,
                                                  for_spins_state)))
              return res;
          }
      }

    // If there is one or more no-invisible transition in the reduced set,
    // all enabled transition are returned.
    if (invisible_)
      for (unsigned i = 0; i < enabled.size(); ++i)
        {
          if (std::find(t_s.begin(), t_s.end(), enabled[i]) != t_s.end()
              && !m_invisible[enabled[i]])
            return res;
        }

    // Compute intersection between t_s and enabled
    for (unsigned i = 0; i < enabled.size(); ++i)
      if (std::find(t_s.begin(), t_s.end(), enabled[i]) == t_s.end())
        res[i] = false;

    return res;
  }

  bool porinfos::non_maybecoenabled(int t1, int t2)
  {
    unsigned t1_guards_size = m_guards[t1].size();
    unsigned t2_guards_size = m_guards[t2].size();
    for (unsigned i_t2 = 0; i_t2 < t2_guards_size; ++i_t2)
      for (unsigned i_t1 = 0; i_t1 < t1_guards_size; ++i_t1)
        if (!m_mbc[m_guards[t2][i_t2]][m_guards[t1][i_t1]])
          return false;
    return true;
  }

  unsigned porinfos::transitions()
  {
    return transitions_;
  }

  unsigned porinfos::variables()
  {
    return variables_;
  }

  void porinfos::dump_read_dependency()
  {
    std::cout << "Read Dependency Matrix" << std::endl;
    for (unsigned i = 0; i < transitions_; ++i)
      {
        std::cout << "   ";
        for (unsigned j = 0; j < variables_; ++j)
          std::cout << (m_read[i][j]? "+" : "-");
        std::cout << std::endl;
      }
    std::cout << std::endl;
  }

  void porinfos::dump_write_dependency()
  {
    std::cout << "Write Dependency Matrix" << std::endl;
    for (unsigned i = 0; i < transitions_; ++i)
      {
        std::cout << "   ";
        for (unsigned j = 0; j < variables_; ++j)
          std::cout << (m_write[i][j]? "+" : "-");
        std::cout << std::endl;
      }
    std::cout << std::endl;
  }

  void porinfos::dump_nes_guards()
  {
    std::cout << "NES Matrix" << std::endl;
    for (unsigned i = 0; i < guards_; ++i)
      {
        std::cout << "   ";
        for (unsigned j = 0; j < transitions_; ++j)
          std::cout << (m_nes[i][j]? "+" : "-");
        std::cout << std::endl;
      }
    std::cout << std::endl;
  }

  void porinfos::dump_mbc_guards()
  {
    std::cout << "May be coenabled Matrix" << std::endl;
    for (unsigned i = 0; i < guards_; ++i)
      {
        std::cout << "   ";
        for (unsigned j = 0; j < guards_; ++j)
          std::cout << (m_mbc[i][j]? "+" : "-");
        std::cout << std::endl;
      }
    std::cout << std::endl;
  }

  porinfos_stats& porinfos::stats()
  {
    return stats_;
  }
}
