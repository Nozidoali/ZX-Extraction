from config import *
from extract import *

if __name__ == "__main__":
    circuit = extract_qc(MODE)
    n_cnots = sum(1 for g in circuit.gates if g.name in ('CNOT','CZ')) # 2-qubit gates
    print("Extracted circuit with", len(circuit.gates), "gates (", n_cnots, "CNOTs/CZs )")