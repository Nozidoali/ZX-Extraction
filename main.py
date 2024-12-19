import json
from config import *
from typing import List
from dataclasses import dataclass

@dataclass
class QGate:
    name: str
    qubits: List[int]
    params: List[float]

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


class BaselineExtractor:
    def __init__(self, graph_json: str):
        from pyzx import Graph
        self.graph = Graph.from_json(graph_json)
    
    def run(self):
        from pyzx import extract_circuit
        circuit = extract_circuit(self.graph).to_basic_gates()
        return sum(1 for g in circuit.gates if g.name in ('CNOT','CZ')) # 2-qubit gates

def extract_qc():
    if MODE == "baseline":
        extractor = BaselineExtractor(open(BENCHMARK_PATH).read())
    else:
        graph = json.load(open(BENCHMARK_PATH))
        extractor = Extractor(graph)
    n_cnots = extractor.run()
    print(n_cnots)