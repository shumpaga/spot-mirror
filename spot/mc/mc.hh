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

#include <functional>
#include <string>
#include <thread>
#include <tuple>
#include <vector>
#include <utility>
#include <spot/kripke/kripke.hh>
#include <spot/mc/ec.hh>
#include <spot/mc/deadlock.hh>
#include <spot/mc/bloemen.hh>
#include <spot/misc/common.hh>
#include <spot/misc/timer.hh>

namespace spot
{
  /// \brief Check for the emptiness between a system and a twa.
  /// Return a pair containing a boolean indicating wether a counterexample
  /// has been found and a string representing the counterexample if the
  /// computation have been required
  template<typename kripke_ptr, typename State, typename Iterator,
      typename Hash, typename Equal>
    static std::tuple<bool, std::string, std::vector<istats>>
    modelcheck(kripke_ptr sys, spot::twacube_ptr twa, bool compute_ctrx = false)
    {
      // Must ensure that the two automata are working on the same
      // set of atomic propositions.
      SPOT_ASSERT(sys->get_ap().size() == twa->get_ap().size());
      for (unsigned int i = 0; i < sys->get_ap().size(); ++i)
        SPOT_ASSERT(sys->get_ap()[i].compare(twa->get_ap()[i]) == 0);

      bool stop = false;
      std::vector<ec_renault13lpar<State, Iterator, Hash, Equal>> ecs;
      for (unsigned i = 0; i < sys->get_threads(); ++i)
        ecs.push_back(
        { *sys, twa, i, stop });

      std::vector<std::thread> threads;
      for (unsigned i = 0; i < sys->get_threads(); ++i)
        threads.push_back(
            std::thread(&ec_renault13lpar<State, Iterator, Hash, Equal>::run,
                        &ecs[i]));

      for (unsigned i = 0; i < sys->get_threads(); ++i)
        threads[i].join();

      bool has_ctrx = false;
      std::string trace = "";
      std::vector<istats> stats;
      for (unsigned i = 0; i < sys->get_threads(); ++i)
      {
        has_ctrx |= ecs[i].counterexample_found();
        if (compute_ctrx && ecs[i].counterexample_found() &&
            trace.compare("") == 0)
          trace = ecs[i].trace(); // Pick randomly one !
        stats.push_back(ecs[i].stats());
      }
      return std::make_tuple(has_ctrx, trace, stats);
    }

  /// \bief Check wether the system contains a deadlock. The algorithm
  /// spawns multiple threads performing a classical swarming DFS. As
  /// soon one thread detects a deadlock all the other threads are stopped.
  template<typename kripke_ptr, typename State, typename Iterator,
      typename Hash, typename Equal>
    static std::tuple<bool, std::vector<deadlock_stats>, spot::timer_map>
    has_deadlock(kripke_ptr sys)
    {
      spot::timer_map tm;
      using algo_name = spot::swarmed_deadlock<State, Iterator, Hash, Equal>;

      unsigned nbth = sys->get_threads();
      typename algo_name::shared_map map;
      std::atomic<bool> stop(false);

      tm.start("Initialisation");
      std::vector<algo_name*> swarmed(nbth);
      for (unsigned i = 0; i < nbth; ++i)
        swarmed[i] = new algo_name(*sys, map, i, stop);
      tm.stop("Initialisation");

      std::mutex iomutex;
      std::atomic<bool> barrier(true);
      std::vector<std::thread> threads(nbth);
      for (unsigned i = 0; i < nbth; ++i)
      {
        threads[i] = std::thread([&swarmed, &iomutex, i, & barrier]
        {
#if defined(unix) || defined(__unix__) || defined(__unix)
                                 {
                                   std::lock_guard<std::mutex> iolock(iomutex);
                                   std::cout << "Thread #" << i
                                   << ": on CPU " << sched_getcpu() << '\n';
                                 }
#endif

                                 // Wait all threads to be instanciated.
            while (barrier)
            continue;
            swarmed[i]->run();
          });

#if defined(unix) || defined(__unix__) || defined(__unix)
        //  Pins threads to a dedicated core.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0)
        {
          std::lock_guard<std::mutex> iolock(iomutex);
          std::cerr << "Error calling pthread_setaffinity_np: " << rc << '\n';
        }
#endif
      }

      tm.start("Run");
      barrier.store(false);

      for (auto& t : threads)
        t.join();
      tm.stop("Run");

      std::vector<deadlock_stats> stats;
      bool has_deadlock = false;
      for (unsigned i = 0; i < sys->get_threads(); ++i)
      {
        has_deadlock |= swarmed[i]->has_deadlock();
        stats.push_back(swarmed[i]->stats());
      }

      for (unsigned i = 0; i < nbth; ++i)
        delete swarmed[i];

      return std::make_tuple(has_deadlock, stats, tm);
    }

  /// \brief Perform the SCC computation algorithm of bloemen.16.ppopp
  template<typename kripke_ptr, typename State, typename Iterator,
      typename Hash, typename Equal>
    static std::pair<std::vector<bloemen_stats>, spot::timer_map>
    bloemen(kripke_ptr sys)
    {
      spot::timer_map tm;
      using algo_name = spot::swarmed_bloemen<State, Iterator, Hash, Equal>;
      using uf_name = spot::iterable_uf<State, Hash, Equal>;

      unsigned nbth = sys->get_threads();
      typename uf_name::shared_map map;

      tm.start("Initialisation");
      std::vector<algo_name*> swarmed(nbth);
      std::vector<uf_name*> ufs(nbth);
      for (unsigned i = 0; i < nbth; ++i)
      {
        ufs[i] = new uf_name(map, i);
        swarmed[i] = new algo_name(*sys, *ufs[i], i);
      }
      tm.stop("Initialisation");

      std::mutex iomutex;
      std::atomic<bool> barrier(true);
      std::vector<std::thread> threads(nbth);
      for (unsigned i = 0; i < nbth; ++i)
      {
        threads[i] = std::thread([&swarmed, &iomutex, i, & barrier]
        {
#if defined(unix) || defined(__unix__) || defined(__unix)
                                 {
                                   std::lock_guard<std::mutex> iolock(iomutex);
                                   std::cout << "Thread #" << i
                                   << ": on CPU " << sched_getcpu() << '\n';
                                 }
#endif

                                 // Wait all threads to be instanciated.
            while (barrier)
            continue;
            swarmed[i]->run();
          });

#if defined(unix) || defined(__unix__) || defined(__unix)
        //  Pins threads to a dedicated core.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0)
        {
          std::lock_guard<std::mutex> iolock(iomutex);
          std::cerr << "Error calling pthread_setaffinity_np: " << rc << '\n';
        }
#endif
      }

      tm.start("Run");
      barrier.store(false);

      for (auto& t : threads)
        t.join();
      tm.stop("Run");

      std::vector<bloemen_stats> stats;
      for (unsigned i = 0; i < sys->get_threads(); ++i)
        stats.push_back(swarmed[i]->stats());

      for (unsigned i = 0; i < nbth; ++i)
        delete swarmed[i];

      return std::make_pair(stats, tm);
    }
}
