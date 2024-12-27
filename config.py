import os, glob
QASM_DIR = './benchmarks/circuits'
GRAPH_DIR = './benchmarks/graphs'
QASMS = glob.glob(os.path.join(QASM_DIR, "*.qasm"))
GRAPHS = glob.glob(os.path.join(GRAPH_DIR, "*.graph"))

UP_TO_PERM = True
TIMEOUT = 60