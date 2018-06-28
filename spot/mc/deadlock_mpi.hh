// -*- coding: utf-8 -*-
// Copyright (C) 2015, 2016, 2017 Laboratoire de Recherche et
// Developpement de l'Epita
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

#include <atomic>
#include <chrono>
#include <bricks/brick-hashset>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <spot/misc/common.hh>
#include <spot/kripke/kripke.hh>
#include <spot/misc/fixpool.hh>
#include <spot/misc/timer.hh>

#include <spot/mpi/mpi.hh>

#define SEUIL 32

namespace spot
{
  /// \brief This object is returned by the algorithm below
  struct SPOT_API deadlock_stats
  {
    unsigned states;            ///< \brief Number of states visited
    unsigned transitions;       ///< \brief Number of transitions visited
    unsigned instack_dfs;       ///< \brief Maximum DFS stack
    bool has_deadlock;          ///< \brief Does the model contains a deadlock
    unsigned walltime;          ///< \brief Walltime for this thread in ms
  };

  /// \brief This class aims to explore a model to detect wether it
  /// contains a deadlock. This deadlock detection performs a DFS traversal
  /// sharing information shared among multiple threads.
  template<typename State, typename SuccIterator, typename StateHash,
      typename StateEqual>
    class swarmed_deadlock
    {
      /// \brief Describes the status of a state
      enum st_status
      {
        UNKNOWN = 1,    // First time this state is discoverd by this thread
        OPEN = 2,       // The state is currently processed by this thread
        CLOSED = 4,     // All the successors of this state have been visited
      };

      /// \brief Describes the structure of a shared state
      struct deadlock_pair
      {
        State st;                 ///< \brief the effective state
        int* colors;              ///< \brief the colors (one per thread)
      };

      /// \brief The haser for the previous state.
      struct pair_hasher
      {
        pair_hasher(const deadlock_pair&)
        {
        }

        pair_hasher() = default;

        brick::hash::hash128_t
        hash(const deadlock_pair& lhs) const
        {
          StateHash hash;
          // Not modulo 31 according to brick::hashset specifications.
          unsigned u = hash(lhs.st) % (1 << 30);
          return
          { u, u};
        }

        bool
        equal(const deadlock_pair& lhs, const deadlock_pair& rhs) const
        {
          StateEqual equal;
          return equal(lhs.st, rhs.st);
        }
      };

    public:

      ///< \brief Shortcut to ease shared map manipulation
      using shared_map = brick::hashset::FastConcurrent <deadlock_pair,
      pair_hasher>;

      swarmed_deadlock(kripkecube<State, SuccIterator>& sys, shared_map& map,
                       unsigned tid, std::atomic<bool>& stop) :
          sys_(sys), tid_(tid), map_(map), nb_th_(
              std::thread::hardware_concurrency()), p_(
              sizeof(int) * std::thread::hardware_concurrency()), stop_(stop)
      {
        SPOT_ASSERT(is_a_kripkecube(sys));
      }

      virtual
      ~swarmed_deadlock()
      {
      }

      void
      setup()
      {
        tm_.start("DFS thread " + std::to_string(tid_));
      }

      bool
      push(State s)
      {
        // Prepare data for a newer allocation
        int* ref = (int*)p_.allocate();
        for (unsigned i = 0; i < nb_th_; ++i)
          ref[i] = UNKNOWN;

        // Try to insert the new state in the shared map.
        auto it = map_.insert(
        { s, ref });
        bool b = it.isnew();

        // Insertion failed, delete element
        // FIXME Should we add a local cache to avoid useless allocations?
        if (!b)
          p_.deallocate(ref);

        // The state has been mark dead by another thread
        for (unsigned i = 0; !b && i < nb_th_; ++i)
          if (it->colors[i] == static_cast<int>(CLOSED))
            return false;

        // The state has already been visited by the current thread
        if (it->colors[tid_] == static_cast<int>(OPEN))
          return false;

        // Keep a ptr over the array of colors
        refs_.push_back(it->colors);

        // Mark state as visited.
        it->colors[tid_] = OPEN;
        ++states_;
        return true;
      }

      bool
      pop()
      {
        // Track maximum dfs size
        dfs_ = todo_.size() > dfs_ ? todo_.size() : dfs_;

        // Don't avoid pop but modify the status of the state
        // during backtrack
        refs_.back()[tid_] = CLOSED;
        refs_.pop_back();
        return true;
      }

      void
      finalize()
      {
        stop_ = true;
        tm_.stop("DFS thread " + std::to_string(tid_));
      }

      unsigned
      states()
      {
        return states_;
      }

      unsigned
      transitions()
      {
        return transitions_;
      }

      void
      run()
      {
        setup();
        State initial = sys_.initial(tid_);
        if (SPOT_LIKELY(push(initial)))
        {
          todo_.push_back(
          { initial, sys_.succ(initial, tid_), transitions_ });
        }
        while (!todo_.empty() && !stop_.load(std::memory_order_relaxed))
        {
          if (todo_.back().it->done())
          {
            if (SPOT_LIKELY(pop()))
            {
              deadlock_ = todo_.back().current_tr == transitions_;
              if (deadlock_)
                break;
              sys_.recycle(todo_.back().it, tid_);
              todo_.pop_back();
            }
          }
          else
          {
            ++transitions_;
            State dst = todo_.back().it->state();

            if (SPOT_LIKELY(push(dst)))
            {
              todo_.back().it->next();
              todo_.push_back(
              { dst, sys_.succ(dst, tid_), transitions_ });
            }
            else
            {
              todo_.back().it->next();
            }
          }
        }
        finalize();
      }

      bool
      has_deadlock()
      {
        return deadlock_;
      }

      unsigned
      walltime()
      {
        return tm_.timer("DFS thread " + std::to_string(tid_)).walltime();
      }

      deadlock_stats
      stats()
      {
        return
        { states(), transitions(), dfs_, has_deadlock(), walltime()};
      }

      /************************************************************************
       *                                   MPI                                *
       ************************************************************************/

      bool
      push_mpi(State s)
      {
        int* ref = (int*)p_.allocate();

        for (unsigned i = 0; i < nb_th_; ++i)
          ref[i] = UNKNOWN;

        auto it = map_.insert(
        { s, ref });
        bool b = it.isnew();

        if (!b)
          p_.deallocate(ref);

        // The state has been mark dead by another thread
        for (unsigned i = 0; !b && i < nb_th_; ++i)
          if (it->colors[i] == static_cast<int>(CLOSED))
            return false;

        // The state has already been visited by the current thread
        if (it->colors[tid_] == static_cast<int>(OPEN))
          return false;

        // Mark state as visited.
        it->colors[tid_] = CLOSED;
        return true;
      }

      void
      run_mpi(struct spot::mpi::attributes_& process_attributes)
      {
        int rank = 0;
        int size = 1;
        /* // must be different from 0 ! 0 is reserved for display */
        int deadlock_tag = 512;
        int state_tag = 1024;
        int end_tag = 2048;
        MPI_Request* array_of_deadlock_request;
        int* cpt_state_send;
        int* cpt_state_recv;
        std::vector<State> state_recv;
        int deadlock_flag = 0;
        int state_flag = 0;
        MPI_Status deadlock_status;
        MPI_Status state_status;
        MPI_Message deadlock_handle;
        MPI_Message state_handle;
        bool break_while = false;
        unsigned cpt_pop = 0;

        MPI_Comm_rank(process_attributes.comm_everyone, &rank);
        MPI_Comm_size(process_attributes.comm_everyone, &size);

        array_of_deadlock_request = new MPI_Request[size];
        cpt_state_send = new int[size];
        cpt_state_recv = new int[size];

        for (int i = 0; i < size; i++)
        {
          array_of_deadlock_request[i] = MPI_REQUEST_NULL;
          cpt_state_send[i] = 0;
          cpt_state_recv[i] = 0;
        }

        setup();
        State initial = sys_.initial(tid_);
        if (SPOT_LIKELY(push(initial)))
        {
          todo_.push_back(
          { initial, sys_.succ(initial, tid_), transitions_ });
        }

        while (!todo_.empty() && !stop_.load(std::memory_order_relaxed))
        {
          for (int i = 0; i < size; i++)
          {
            MPI_Improbe(i, deadlock_tag + i + tid_,
                        process_attributes.comm_everyone, &deadlock_flag,
                        &deadlock_handle, &deadlock_status);

            if (deadlock_flag)
            {
              char deadlock_message = '0';

              MPI_Mrecv(&deadlock_message, 1, MPI_CHAR, &deadlock_handle,
                        &deadlock_status);
              deadlock_ = true;
              break_while = true;
              deadlock_flag = 0;
              break;
            }

            MPI_Improbe(i, state_tag + i + tid_,
                        process_attributes.comm_everyone, &state_flag,
                        &state_handle, &state_status);

            if (state_flag)
            {
              std::vector<todo__element> tmp;
              int* s;
              int state_size = 0;

              MPI_Get_count(&state_status, MPI_INT, &state_size);
              s = new int[state_size];
              MPI_Mrecv(s, state_size, MPI_INT, &state_handle, &state_status);
              cpt_state_recv[i]++;

              if (SPOT_LIKELY(push_mpi((State)s)))
              {
                state_recv.push_back((State)s);
                tmp.push_back(
                { s, sys_.succ(s, tid_), 0 });

                while (!tmp.empty())
                {
                  if (tmp.back().it->done())
                  {
                    sys_.recycle(tmp.back().it, tid_);
                    tmp.pop_back();
                  }

                  else
                  {
                    s = tmp.back().it->state();

                    if (SPOT_LIKELY(push_mpi((State)s)))
                    {
                      tmp.back().it->next();
                      tmp.push_back(
                      { s, sys_.succ(s, tid_), 0 });
                    }

                    else
                    {
                      tmp.back().it->next();
                    }
                  }
                }
              }

              state_flag = 0;
            }
          }

          if (break_while)
            break;

          if (todo_.back().it->done())
          {
            int* current = (int*)todo_.back().s;

            cpt_pop++;

            if (SPOT_LIKELY(pop()))
            {
              deadlock_ = todo_.back().current_tr == transitions_;

              if (rank == 0)
                deadlock_ = true;

              if (deadlock_)
              {
                char deadlock_message = 'd';

                for (int i = (rank + 1) % size; i != rank; i = (i + 1) % size)
                {
                  MPI_Isend(&deadlock_message, 1,
                  MPI_CHAR,
                            i, deadlock_tag + rank + tid_,
                            process_attributes.comm_everyone,
                            &array_of_deadlock_request[i]);
                }
                break;
              }

              if (cpt_pop >= SEUIL)
              {
                for (int i = (rank + 1) % size; i != rank; i = (i + 1) % size)
                {
                  MPI_Request request = MPI_REQUEST_NULL;

                  MPI_Isend(current, current[1] + 2,
                  MPI_INT,
                            i, state_tag + rank + tid_,
                            process_attributes.comm_everyone, &request);
                  cpt_state_send[i]++;
                }
                cpt_pop = 0;
              }

              sys_.recycle(todo_.back().it, tid_);
              todo_.pop_back();
            }
          }

          else
          {
            ++transitions_;
            State dst = todo_.back().it->state();

            if (SPOT_LIKELY(push(dst)))
            {
              todo_.back().it->next();
              todo_.push_back(
              { dst, sys_.succ(dst, tid_), transitions_ });
            }
            else
            {
              todo_.back().it->next();
            }
          }
        }

        finalize();

        MPI_Barrier(process_attributes.comm_everyone);

        for (int i = (rank + 1) % size; i != rank; i = (i + 1) % size)
        {
          MPI_Request request = MPI_REQUEST_NULL;

          MPI_Isend(&cpt_state_send[i], 1,
          MPI_INT,
                    i, end_tag + rank + tid_, process_attributes.comm_everyone,
                    &request);
        }

        for (int i = (rank + 1) % size; i != rank; i = (i + 1) % size)
        {
          MPI_Status status;
          int max = 0;

          MPI_Recv(&max, 1,
          MPI_INT,
                   i, end_tag + i + tid_, process_attributes.comm_everyone,
                   &status);

          while (cpt_state_recv[i] < max)
          {
            MPI_Improbe(i, state_tag + i + tid_,
                        process_attributes.comm_everyone, &state_flag,
                        &state_handle, &state_status);

            if (state_flag)
            {
              int* s;
              int state_size = 0;

              MPI_Get_count(&state_status, MPI_INT, &state_size);
              s = new int[state_size];
              MPI_Mrecv(s, state_size, MPI_INT, &state_handle, &state_status);
              cpt_state_recv[i]++;
              state_flag = 0;
              delete[] s;
            }
          }

        }

        for (State& s : state_recv)
        {
          delete[] s;
        }

        for (int i = 0; i < size; i++)
        {
          int test_flag = 0;
          MPI_Status status;

          MPI_Test(&array_of_deadlock_request[i], &test_flag, &status);

          if (!test_flag)
          {
            MPI_Cancel(&array_of_deadlock_request[i]);
            MPI_Request_free(&array_of_deadlock_request[i]);
          }

          test_flag = 0;
        }

        delete[] array_of_deadlock_request;
      }

      /************************************************************************/

    private:
      struct todo__element
      {
        State s;
        SuccIterator* it;
        unsigned current_tr;
      };
      kripkecube<State, SuccIterator>& sys_; ///< \brief The system to check
      std::vector<todo__element> todo_;      ///< \brief The DFS stack
      unsigned transitions_ = 0;         ///< \brief Number of transitions
      unsigned tid_;                     ///< \brief Thread's current ID
      shared_map map_;                       ///< \brief Map shared by threads
      spot::timer_map tm_;                   ///< \brief Time execution
      unsigned states_ = 0;                  ///< \brief Number of states
      unsigned dfs_ = 0;                     ///< \brief Maximum DFS stack size
      /// \brief Maximum number of threads that can be handled by this algorithm
      unsigned nb_th_ = 0;
      fixed_size_pool p_;                    ///< \brief State Allocator
      bool deadlock_ = false;                ///< \brief Deadlock detected?
      std::atomic<bool>& stop_;              ///< \brief Stop-the-world boolean
      /// \brief Stack that grows according to the todo stack. It avoid multiple
      /// concurent access to the shared map.
      std::vector<int*> refs_;
    };
}
