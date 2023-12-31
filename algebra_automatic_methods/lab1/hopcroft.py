from finite_automation import FiniteAutomaton


def hopcroft_step(P: list, L: list, state_transitions: list, alpabet: list):
    C = L[0]
    L.remove(C)

    P1 = P.copy()

    for i in range(len(alpabet)):
        for cls in P:
            states = []

            for j in range(len(state_transitions)):
                if state_transitions[j][i] in C and j in cls:
                    states.append(j)

            if len(states) > 0:
                B2 = states
                B1 = [x for x in cls if x not in B2]

                if cls in P1:
                    P1.remove(cls)

                if len(B1) > 0 and B1 not in P1:
                    P1.append(B1)
                if len(B2) > 0 and B2 not in P1:
                    P1.append(B2)

                if cls in L:
                    L.remove(cls)
                    if len(B1) > 0:
                        L.append(B1)
                    if len(B2) > 0:
                        L.append(B2)
                else:
                    if len(B1) <= len(B2) and len(B1) > 0:
                        L.append(B1)
                    elif len(B1) > len(B2) and len(B2) > 0:
                        L.append(B2)
        P = P1.copy()

    return P1, L


def find_state_index(P: list, transitions: list) -> int | None:
    for i, s in enumerate(P):
        if set(s).union(set(transitions)) == set(s):
            return i
    return None


def form_transition_table(P: list, state_transitions: list, alphabet: list):
    table = []

    for state in P:
        row = []

        for i in range(len(alphabet)):
            transitions = [state_transitions[j][i] for j in state]
            idx = find_state_index(P, transitions)
            row.append(idx if idx is not None else -1)

        table.append(row)

    return table


def update_states(P: list, states: list, begin):
    final_states = set()
    begin_state = begin

    for i, s in enumerate(P):
        if set(s).intersection(states):
            final_states.add(i)
        if begin in s:
            begin_state = i

    return list(final_states), begin_state


def hopcroft(states_count, alphabet, state_transitions, final_states):
    P = [[x for x in range(states_count) if x not in final_states], final_states]
    L = []

    if len(P[0]) < len(P[1]):
        L.append(P[0])
    else:
        L.append(P[1])

    counter = 1
    while len(L) > 0:
        P, L = hopcroft_step(P, L, state_transitions, alphabet)
        L.sort(key=lambda x: len(x))
        counter += 1

    print(f"Minimized P: {P}\n")

    new_state_transitions = form_transition_table(P, state_transitions, alphabet)

    return P, new_state_transitions


def main():
    automation = FiniteAutomaton.load_from(file=FILE_NAME)

    states_count, alphabet, state_transitions, final_states = len(
        automation.states), automation.alphabet, automation.transition_matrix, automation.finals
    print(automation, "\n")

    P, table = hopcroft(states_count, alphabet, state_transitions, final_states)
    new_final_states, begin = update_states(P, final_states, automation.begin)
    min_automation = FiniteAutomaton.from_matrix(alphabet, table, new_final_states, begin)

    print(min_automation)
    min_automation.save(f"min_{FILE_NAME}")


FILE_NAME = "data2.json"
main()
