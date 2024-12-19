import json
from config import *
from typing import Dict, List

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
        

graph = json.load(open(BENCHMARK_PATH))
extractor = Extractor(graph)
print(extractor.pis)
print(extractor.pos)
print(extractor.nodes)
print(extractor.edges)
