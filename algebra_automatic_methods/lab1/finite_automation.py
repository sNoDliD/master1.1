import ast
import json


class FiniteAutomaton:
    def __init__(self, states, alphabet, transitions, begin, finals):
        """
        :param states: set of states A
        :param alphabet: set of signs X
        :param transitions: transition function f
        :param begin: begin state
        :param finals: set of final stales F
        """
        self.states = states  # A
        self.alphabet = alphabet  # X
        self.transitions = transitions  # f
        self.begin = begin  # a0
        self.finals = finals  # F

    @property
    def transition_matrix(self) -> list[list]:
        return [[self(state, sign, default=-1) for sign in self.alphabet] for state in self.states]

    @classmethod
    def from_matrix(cls, alphabet, matrix: list[list], finals, begin):
        states = [i for i in range(len(matrix))]
        return cls(
            states,
            alphabet,
            {(i, alphabet[j]): item
             for i, row in enumerate(matrix)
             for j, item in enumerate(row)
             if item != -1},
            begin,
            finals,
        )

    def __str__(self):
        transitions = "\n".join(f"{key} -> {self.transitions[key]}" for key in self.transitions)
        return f"A={self.states}\nX={self.alphabet}\nf=\n{transitions}\na0={self.begin}\nF={self.finals}"

    def __call__(self, state, sign, default=None):
        return self.transitions.get((state, sign), default)

    def to_json(self):
        automaton_data = {
            "states": list(self.states),
            "alphabet": list(self.alphabet),
            "transitions": {str(key): value for key, value in self.transitions.items()},
            "begin": self.begin,
            "finals": list(self.finals)
        }
        return json.dumps(automaton_data, indent=4)

    @classmethod
    def from_json(cls, json_str):
        automaton_data = json.loads(json_str)
        states = list(automaton_data["states"])
        alphabet = list(automaton_data["alphabet"])

        transitions = {ast.literal_eval(key): value for key, value in automaton_data["transitions"].items()}
        begin = automaton_data["begin"]
        finals = list(automaton_data["finals"])
        return cls(states, alphabet, transitions, begin, finals)

    @classmethod
    def load_from(cls, file: str):
        with open(file) as file:
            return cls.from_json('\n'.join(file.readlines()))

    def save(self, file: str):
        with open(file, 'wt') as file:
            return file.write(self.to_json())
