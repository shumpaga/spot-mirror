import spot

a1 = spot.automaton('''
HOA: v1  name: "aut"  States: 4  Start: 0  AP: 0
Acceptance: 4 Fin(0) & (Inf(1) | (Fin(2) & Inf(3)))
--BODY--
State: 0  [t] 0 {0}  [t] 1 {0 2}
State: 1  [t] 2
State: 2  [t] 1      [t] 0 {0}    [t] 3 {3}
State: 3  [t] 2 {1}  [t] 0
--END--''')

a2 = spot.automaton('''
HOA: v1 States: 7 Start: 0 AP: 3 "a" "b" "c" Acceptance: 3 Fin(0) & Fin(2)
& Inf(1) properties: trans-labels explicit-labels trans-acc --BODY--
State: 0 [2] 1 [!0&!2] 2 [0&!2] 3 [0&!2] 4 State: 1 [!0 | 1] 1 {1 2}
[0&!1] 5 {1} State: 2 [t] 2 {0} State: 3 [!0&!2] 2 [0&!2] 3 State: 4
[0&!1&2] 1 [!0&2 | 1&2] 1 {2} [!0&!2] 6 {2} [0&!1&!2] 4 [0&1&!2] 4 {2}
State: 5 [!0 | 1] 1 {2} [0&!1] 5 State: 6 [0&!1] 6 {0} [!0 | 1] 6 {0
2} --END--''')

a3 = spot.automaton('''
HOA: v1 States: 11 Start: 0 AP: 3 "a" "b" "c" Acceptance: 5 (Fin(0)
| Inf(1)) & (Fin(2) | Inf(3)) & Fin(4) properties: trans-labels
explicit-labels trans-acc --BODY-- State: 0 [0&!1&!2] 1 {1} [0&!1&!2]
2 {1} [!0&!2] 3 {1} [0&1&!2] 4 {1} [0&1&!2] 5 {1} [0&!1&2] 6 {1} [!0&2 |
1&2] 7 {1} State: 1 [0&!1&!2] 1 {0 2} [!0&!2] 3 {0 2} [0&1&!2] 4 {0 2}
State: 2 [0&!1&!2] 2 {0 2} [!0&!2] 8 {0 2 4} [0&1&!2] 5 {0 2 4} [0&!1&2]
6 {0 2} [!0&2 | 1&2] 7 {0 2 4} State: 3 [0&!1] 9 {1 2} [!0 | 1] 3 {1 2}
State: 4 [0&!1&!2] 1 {1 2} [!0&!2] 3 {1 2} [0&1&!2] 4 {1 2} State: 5
[0&!1&!2] 2 {1 2} [!0&!2] 8 {1 2 4} [0&1&!2] 5 {1 2 4} [0&!1&2] 6 {1 2}
[!0&2 | 1&2] 7 {1 2 4} State: 6 [0&!1] 6 {0 3} [!0 | 1] 7 {0 3 4} State:
7 [0&!1] 6 {1 3} [!0 | 1] 7 {1 3 4} State: 8 [0&!1] 10 {1 2} [!0 | 1]
8 {1 2 4} State: 9 [0&!1] 9 {0 2} [!0 | 1] 3 {0 2} State: 10 [0&!1]
10 {0 2} [!0 | 1] 8 {0 2 4} --END--''')

a4 = spot.automaton('''
HOA: v1 States: 8 Start: 0 AP: 3 "a" "b" "c" Acceptance: 6 ((Fin(0) &
Inf(1)) | (Fin(2) & Inf(3))) & Fin(4) & Inf(5) properties: trans-labels
explicit-labels state-acc complete properties: deterministic --BODY--
State: 0 {2 4} [0&2] 1 [!0&2] 2 [!0&!1&!2] 3 [!0&1&!2] 4 [0&!2] 5 State:
1 {1 2 4} [!0] 6 [0] 1 State: 2 {1 2 5} [!0] 6 [0] 1 State: 3 {1 2 4}
[t] 3 State: 4 {1 2 4} [!0&2] 6 [0&2] 1 [!1&!2] 3 [0&1&!2] 4 [!0&1&!2]
7 State: 5 {1 2 4} [!0&2] 6 [0&2] 1 [!0&!1&!2] 3 [0&!2] 5 [!0&1&!2]
7 State: 6 {2 5} [!0] 6 [0] 1 State: 7 {3 4} [!0&2] 6 [0&2] 1 [!1&!2]
3 [0&1&!2] 4 [!0&1&!2] 7 --END--''')

a5 = spot.automaton("""
HOA: v1 States: 10 Start: 0 AP: 2 "a" "b" Acceptance: 5 ((Inf(4)
| Fin(0)) & Fin(1)) | (Inf(2)&Inf(3)) properties: trans-labels
explicit-labels trans-acc --BODY-- State: 0 [0&!1] 2 State: 1 [!0&!1]
9 {3} [!0&1] 7 {0 1 2} State: 2 [!0&1] 1 [0&!1] 3 {1 2} [!0&1] 2 {3}
[0&!1] 5 {1 4} State: 3 [0&!1] 5 {0 3} [!0&!1] 9 [!0&!1] 2 [0&!1] 4 State:
4 [0&!1] 3 [!0&!1] 6 {0 1 2} State: 5 [!0&!1] 7 {0} [0&!1] 4 [0&!1] 1
[0&!1] 8 {0} State: 6 [!0&!1] 5 {2} [0&!1] 4 {3} State: 7 [!0&1] 5 [0&!1]
2 {2} [!0&1] 4 State: 8 [0&!1] 1 {1} [!0&!1] 6 {0} State: 9 [0&1] 0 {1}
[0&1] 3 {2} [0&1] 7 {1} [0&!1] 8 {2 3 4} --END--""")

a6 = spot.automaton("""
HOA: v1 States: 10 Start: 0 AP: 2 "a" "b" Acceptance: 5 (Inf(1) | (Fin(0)
& Inf(4)) | Fin(2)) & Fin(3) properties: trans-labels explicit-labels
trans-acc --BODY-- State: 0 [0&1] 9 {3} [!0&!1] 0 [0&!1] 5 {0 1} State:
1 [0&!1] 9 {4} [0&1] 8 {3} State: 2 [!0&!1] 8 {0} [!0&1] 6 {2 4} [0&1]
2 [!0&1] 7 State: 3 [0&!1] 2 {0 4} [!0&!1] 3 {1} [!0&1] 4 {0} State: 4
[0&!1] 5 {2} [0&1] 0 [!0&1] 1 {0} State: 5 [!0&!1] 0 [!0&!1] 6 State: 6
[0&1] 3 {2} [!0&1] 1 [0&1] 2 {0 1 3 4} State: 7 [0&1] 1 [!0&1] 7 {0 2}
State: 8 [!0&1] 7 [!0&!1] 9 {0} State: 9 [0&1] 8 {0} [0&!1] 5 [0&!1]
1 --END--""")

a7 = spot.automaton("""
HOA: v1 States: 10 Start: 0 AP: 2 "a" "b" acc-name: Rabin 3 Acceptance:
6 (Fin(0) & Inf(1)) | (Fin(2) & Inf(3)) | (Fin(4) & Inf(5)) properties:
trans-labels explicit-labels trans-acc --BODY-- State: 0 [0&!1] 8 {0}
[!0&!1] 6 {0 1} State: 1 [!0&1] 4 {2} [0&1] 8 {2} State: 2 [0&1] 6 {1 4}
[!0&!1] 3 {1 4} State: 3 [!0&!1] 8 {2 4} [0&1] 4 State: 4 [!0&!1] 8 {4}
[!0&!1] 7 State: 5 [!0&!1] 2 {0 5} [!0&!1] 8 {0 4} [!0&!1] 9 {4} State:
6 [!0&1] 1 {2 3 4} State: 7 [0&!1] 5 {0} [0&!1] 7 State: 8 [!0&1] 4 {0 2}
State: 9 [0&1] 3 {4} [!0&1] 5 {4} --END--""")

a8 = spot.automaton('''
HOA: v1 States: 10 Start: 0 AP: 2 "a" "b" Acceptance: 6 Fin(5) &
((Fin(1) & (Inf(3) | Inf(4))) | Fin(0) | Fin(2)) properties: trans-labels
explicit-labels trans-acc --BODY-- State: 0 [0&1] 8 {0} [0&!1] 6 {2}
State: 1 [!0&1] 9 {0 4 5} State: 2 [!0&1] 1 State: 3 [0&!1] 3 {2}
[0&1] 4 {3 5} State: 4 [0&1] 7 {5} [0&!1] 9 {2} [!0&1] 0 {0 2} State:
5 [!0&1] 1 [!0&1] 3 {2 3} State: 6 [0&!1] 8 {1 2 5} [!0&1] 7 {3} State:
7 [0&1] 2 {0} [!0&1] 5 State: 8 [0&!1] 3 {4 5} State: 9 [!0&1] 3 {1 2}
[0&1] 1 {4} [0&!1] 5 {2} --END--''')

a9 = spot.automaton("""
HOA: v1 States: 10 Start: 0 AP: 2 "a" "b" acc-name: Streett 3 Acceptance:
6 (Fin(0) | Inf(1)) & (Fin(2) | Inf(3)) & (Fin(4) | Inf(5)) properties:
trans-labels explicit-labels trans-acc deterministic --BODY-- State: 0
[0&1] 1 [0&!1] 9 {0 5} State: 1 [0&!1] 5 State: 2 [!0&!1] 4 {1} State: 3
[!0&!1] 8 {0} State: 4 [0&1] 6 {0 3} State: 5 [!0&!1] 7 State: 6 [!0&1]
4 State: 7 [!0&!1] 3 {2 5} State: 8 [0&!1] 1 {2} [!0&!1] 2 {2} State:
9 [!0&1] 6 {2 4} --END--""")

a10 = spot.automaton("""
HOA: v1 States: 2 Acceptance: 4 (Fin(0)|Fin(1))&(Fin(2)|Fin(3)) Start: 0
AP: 0 --BODY-- State: 0 [t] 0 {0 2 3} [t] 1 {1} State: 1 [t] 0 {2}
[t] 1 {3 0 1} --END--""")

a11 = spot.automaton("""
HOA: v1 States: 2 Acceptance: 6 (Fin(0)|Fin(1))&(Fin(2)|Fin(3))&
(Fin(4)|Fin(5)) Start: 0 AP: 0 --BODY-- State: 0 [t] 0 {0 2 3} [t]
1 {1 4} State: 1 [t] 0 {2 5} [t] 1 {3 0 1} --END--""")

# From issue #360.
a360 = spot.automaton("""HOA: v1
States: 2
Start: 0
AP: 2 "a" "b"
Acceptance: 8 Fin(5) & (Inf(4) | (Fin(3) & (Inf(2) | (Fin(1) & Inf(0))))) &
              (Inf(6) | Inf(7)) & (Fin(6)|Fin(7))
properties: trans-labels explicit-labels trans-acc complete
properties: deterministic
--BODY--
State: 0
[0&1] 0 {4 6 7}
[0&!1] 1 {0 6}
[!0&1] 0 {3 7}
[!0&!1] 0 {0}
State: 1
[0&1] 0 {4 6 7}
[0&!1] 1 {3 6}
[!0&1] 0 {4 7}
[!0&!1] 1 {0}
--END--""")




def generic_emptiness2_rec(aut):
    spot.cleanup_acceptance_here(aut, False)
    # Catching 'false' acceptance here is an optimization that could be removed.
    if aut.acc().is_f():
        return True
    # Catching Fin-less acceptance here to use a regular emptiness check is an
    # optimization.  This "if" block could be removed without breaking the
    # algorithm.
    if not aut.acc().uses_fin_acceptance():
        return aut.is_empty()
    si = spot.scc_info(aut, True)
    acc_scc = si.one_accepting_scc()
    if acc_scc >= 0:
        return False
    nscc = si.scc_count()
    # Now recurse in all non-rejecting SCC
    for scc in range(nscc):
        if not si.is_rejecting_scc(scc):
            acc = aut.acc()
            sets = si.acc_sets_of(scc)
            acc = acc.restrict_to(sets)
            # Do we have any unit Fin?
            fu = acc.fin_unit()
            if fu:
                for part in si.split_on_sets(scc, fu):
                    if not generic_emptiness2(part):
                        return False
            else:
                # Find some Fin set, we necessarily have one, otherwise the SCC
                # would have been found to be either rejecting or accepting.
                fo = acc.fin_one()
                assert fo >= 0, acc
                for part in si.split_on_sets(scc, [fo]):
                    if not generic_emptiness2(part):
                        return False
                whole = si.split_on_sets(scc, [])[0]
                whole.set_acceptance(acc.force_inf([fo]))
                if not generic_emptiness2(whole):
                    return False
    return True

# The python version of spot.generic_emptiness_check()
def generic_emptiness2(aut):
    old_a = spot.acc_cond(aut.acc())
    res = generic_emptiness2_rec(aut)
    # Restore the original acceptance condition
    aut.set_acceptance(old_a)
    return res

def run_bench(automata):
    for aut in automata:
        # Make sure our three implementation behave identically
        res3 = spot.generic_emptiness_check(aut)
        res2 = spot.remove_fin(aut).is_empty()
        res1 = generic_emptiness2(aut)
        res = str(res1)[0] + str(res2)[0] + str(res3)[0]
        print(res)
        assert res in ('TTT', 'FFF')

run_bench([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a11, a360])
