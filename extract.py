import json
from config import *
from typing import List, Dict, Union, TypeVar, Optional, Tuple, NamedTuple
from pyzx import Circuit, Graph, extract_circuit
from pyzx.circuit import CNOT
from pyzx.extract import clean_frontier, neighbors_of_frontier, remove_gadget, bi_adj, greedy_reduction, column_optimal_swap, apply_cnots, graph_to_swaps, id_simp, filter_duplicate_cnots
from pyzx.linalg import Mat2, Z2
from queue import Queue
import random

VT = TypeVar('VT', bound=int) # The type that is used for representing vertices (e.g. an integer)
ET = TypeVar('ET') # The type used for representing edges (e.g. a pair of integers)

class State:
    def __init__(self, cost: int, index: Tuple[int, ...], rows: List[Z2]):
        self.cost = cost
        self.index = index
        self.rows = rows

    def __hash__(self):
        return hash(sum((v<<i for i, v in enumerate(self.rows))))

    def __repr__(self):
        return f"State({self.cost}, {self.index}, {self.rows})"


def bfs_reduction(m: Mat2) -> Optional[List[Tuple[int, int]]]:
    r = m.rows()
    c = m.cols()
    d = m.data
    if any(sum(r) == 1 for r in d):
        return tuple()
    
    q: Queue[State] = Queue()
    for i in range(r): q.put(State(0, (i,), d[i]))
    visited = {}
    n_iter: int = 0
    while not q.empty():
        s: State = q.get()
        if s in visited: continue
        n_iter += 1
        if n_iter > 1e6: return None
        visited[s] = True
        if sum(s.rows) == 1:
            # row = max([[i, d[i]] for i in s.index], key=lambda x: sum(x[1]))[0]
            cnots = []
            indices = list(s.index[:])
            while len(indices) > 1:
                control = random.choice(indices)
                indices.remove(control)
                target = random.choice(indices)
                cnots.append((control, target))
            return cnots
        for i in range(r):
            if i in s.index: continue
            new_rows = [d[i][j] ^ s.rows[j] for j in range(c)]
            q.put(State(s.cost + 1, s.index + (i,), new_rows))        

class Extractor:
    def __init__(self, graph: dict):
        self.load_graph(graph)

    def load_graph(self, graph: dict):
        self.pis = [n for n, v in graph['wire_vertices'].items() if "input" in v["annotation"]]
        self.pos = [n for n, v in graph['wire_vertices'].items() if "output" in v["annotation"]]
        self.nodes = list(graph["wire_vertices"].keys()) + list(graph["node_vertices"].keys())
        self.node_data = {**{n: v["data"] for n, v in graph["node_vertices"].items()}, **{n: None for n, _ in graph["wire_vertices"].items()}}
        self.edges = {k: [] for k in self.nodes}
        [self.edges[e["src"]].append(e["tgt"]) or self.edges[e["tgt"]].append(e["src"]) for _, e in graph["undir_edges"].items()]

    def run():
        raise NotImplementedError

class BaselineExtractor:
    def __init__(self, graph_json: str):
        self.graph = Graph.from_json(graph_json)
    
    def run(self):
        return extract_circuit(self.graph, up_to_perm=UP_TO_PERM)

class ZXExtractor:
    def __init__(self, graph_json: str):
        self.graph: Graph = Graph.from_json(graph_json)
    
    def run(self):
        optimize_czs = True
        optimize_cnots = 3
        g = self.graph.copy()
        gadgets = {}
        inputs, outputs = g.inputs(), g.outputs()
        c = Circuit(len(outputs))
        for v in g.vertices():
                if g.vertex_degree(v) == 1 and v not in inputs and v not in outputs:
                    n = list(g.neighbors(v))[0]
                    gadgets[n] = v

        qubit_map: Dict[VT,int] = {}
        frontier = []
        for i, o in enumerate(outputs):
            v = list(g.neighbors(o))[0]
            if v in inputs:
                continue
            frontier.append(v)
            qubit_map[v] = i
        while True:
            # preprocessing
            clean_frontier(g, c, frontier, qubit_map, optimize_czs)
            
            # Now we can proceed with the actual extraction
            # First make sure that frontier is connected in correct way to inputs
            neighbor_set = neighbors_of_frontier(g, frontier)
            
            if not frontier:
                break  # No more vertices to be processed. We are done.
            
            # First we check if there is a phase gadget in the way
            if remove_gadget(g, frontier, qubit_map, neighbor_set, gadgets):
                # There was a gadget in the way. Go back to the top
                continue
                
            neighbors = list(neighbor_set)
            m = bi_adj(g, neighbors, frontier)
            if all(sum(row) != 1 for row in m.data):  # No easy vertex
                # print the matrix
                greedy_operations = greedy_reduction(m)
                bfs_operations = bfs_reduction(m)
                
                if bfs_operations is not None:
                    bfs_cnots = [CNOT(target, control) for control, target in bfs_operations]
                    cnots = bfs_cnots
                    greedy_cnots = [CNOT(target, control) for control, target in greedy_operations]
                    # print("-"*80 + "\nBFS reduction\n" + f"{m}\n" + f"gates: {cnots}")
                    # print(f"Greedy reduction: {greedy_cnots}")

                if greedy_operations is None or len(bfs_cnots) > 1:
                    perm = column_optimal_swap(m)
                    perm = {v: k for k, v in perm.items()}
                    neighbors2 = [neighbors[perm[i]] for i in range(len(neighbors))]
                    m2 = bi_adj(g, neighbors2, frontier)
                    if optimize_cnots > 0:
                        cnots = m2.to_cnots(optimize=True)
                    else:
                        cnots = m2.to_cnots(optimize=False)
                    # Since the matrix is not square, the algorithm sometimes introduces duplicates
                    cnots = filter_duplicate_cnots(cnots)
                    # print("-"*80 + f"\nColumn optimal swap\n bfs_cnots: {bfs_cnots}\n" + f"{m}\n" + f"gates: {cnots}")

                    if greedy_operations is not None:
                        m3 = m2.copy()
                        for cnot in cnots:
                            m3.row_add(cnot.target,cnot.control)
                        reductions = sum(1 for row in m3.data if sum(row) == 1)
                        if bfs_cnots and (len(cnots)/reductions > len(bfs_cnots)-0.1):
                            cnots = bfs_cnots
                        else:
                            neighbors = neighbors2
                            m = m2
            else:
                cnots = []

            extracted = apply_cnots(g, c, frontier, qubit_map, cnots, m, neighbors)

        # Outside of loop. Finish up the permutation
        id_simp(g, quiet=True)  # Now the graph should only contain inputs and outputs
        # Since we were extracting from right to left, we reverse the order of the gates
        c.gates = list(reversed(c.gates))
        return graph_to_swaps(g, UP_TO_PERM) + c

def extract_main(filename: str, mode: str = "baseline") -> Circuit:
    if mode == "baseline":
        extractor = BaselineExtractor(open(filename).read())
    elif mode == "zx":
        extractor = ZXExtractor(open(filename).read())
    else:
        extractor = Extractor(json.load(open(filename)))
    return extractor.run().to_basic_gates()

