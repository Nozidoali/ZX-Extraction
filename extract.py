import json
from config import *
from typing import List, Dict, Union, TypeVar
from pyzx import Circuit, Graph, extract_circuit
from pyzx.circuit import CNOT
from pyzx.extract import clean_frontier, neighbors_of_frontier, remove_gadget, bi_adj, greedy_reduction, column_optimal_swap, apply_cnots, graph_to_swaps, id_simp, filter_duplicate_cnots

VT = TypeVar('VT', bound=int) # The type that is used for representing vertices (e.g. an integer)
ET = TypeVar('ET') # The type used for representing edges (e.g. a pair of integers)

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
        return extract_circuit(self.graph)

class ZXExtractor:
    def __init__(self, graph_json: str):
        self.graph: Graph = Graph.from_json(graph_json)
    
    def run(self):
        optimize_czs = True
        optimize_cnots = 2
        g = self.graph.copy()
        gadgets = {}
        inputs, outputs = g.inputs(), g.outputs()
        c = Circuit(len(outputs))
        for v in g.vertices():
                if g.vertex_degree(v) == 1 and v not in inputs and v not in outputs:
                    n = list(g.neighbors(v))[0]
                    gadgets[n] = v

        qubit_map: Dict[VT,int] = dict()
        frontier = []
        for i, o in enumerate(outputs):
            v = list(g.neighbors(o))[0]
            if v in inputs:
                continue
            frontier.append(v)
            qubit_map[v] = i

        czs_saved = 0
        q: Union[float, int]
        
        while True:
            # preprocessing
            czs_saved += clean_frontier(g, c, frontier, qubit_map, optimize_czs)
            
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
                if optimize_cnots > 1:
                    greedy_operations = greedy_reduction(m)
                else:
                    greedy_operations = None

                if greedy_operations is not None:
                    greedy = [CNOT(target, control) for control, target in greedy_operations]
                    cnots = greedy

                if greedy_operations is None or (optimize_cnots == 3 and len(greedy) > 1):
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

                    if greedy_operations is not None:
                        m3 = m2.copy()
                        for cnot in cnots:
                            m3.row_add(cnot.target,cnot.control)
                        reductions = sum(1 for row in m3.data if sum(row) == 1)
                        if greedy and (len(cnots)/reductions > len(greedy)-0.1):
                            cnots = greedy
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
        return graph_to_swaps(g, False) + c

def extract_qc(mode: str = "baseline"):
    if mode == "baseline":
        extractor = BaselineExtractor(open(BENCHMARK_PATH).read())
    elif mode == "zx":
        extractor = ZXExtractor(open(BENCHMARK_PATH).read())
    else:
        extractor = Extractor(json.load(open(BENCHMARK_PATH)))
    return extractor.run().to_basic_gates()
