from config import *
from extract import *
from multiprocessing import Process, Queue
import os
import pandas as pd

def num_2q_gates(circ: Circuit) -> int:
    return sum(1 for g in circ.gates if g.name in ['CNOT', 'CZ'])

def num_t_gates(circ: Circuit) -> int:
    return circ.tcount()

def full_reduce_step(g: Graph) -> None:
    from pyzx.simplify import clifford_simp, gadget_simp, pivot_gadget_simp, interior_clifford_simp
    from pyzx.simplify import id_simp, lcomp_simp, pivot_simp, spider_simp, to_gh, pivot_boundary_simp
    quiet = True
    stats = None
    record = []
    def record_stats(method_name: str) -> None:
        _g = g.copy()
        try:
            circ = extract_circuit(_g, quiet=True)
            record.append({
                'i_step': len(record),
                'num_2q': num_2q_gates(circ),
                'num_t': num_t_gates(circ),
                'method': method_name,
            })
        except Exception as e:
            pass
    spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
    to_gh(g)
    while True:
        i1 = id_simp(g, quiet=quiet, stats=stats); record_stats('id_simp')
        i2 = spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
        i3 = pivot_simp(g, quiet=quiet, stats=stats); record_stats('pivot_simp')
        i4 = lcomp_simp(g, quiet=quiet, stats=stats); record_stats('lcomp_simp')
        if i1+i2+i3+i4==0: break
    pivot_gadget_simp(g, quiet=quiet, stats=stats); record_stats('pivot_gadget_simp')
    while True:
        while True:
            spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
            to_gh(g)
            while True:
                i1 = id_simp(g, quiet=quiet, stats=stats); record_stats('id_simp')
                i2 = spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
                i3 = pivot_simp(g, quiet=quiet, stats=stats); record_stats('pivot_simp')
                i4 = lcomp_simp(g, quiet=quiet, stats=stats); record_stats('lcomp_simp')
                if i1+i2+i3+i4==0: break
            i2 = pivot_boundary_simp(g, quiet=quiet, stats=stats); record_stats('pivot_boundary_simp')
            if i2 == 0:
                break
        i = gadget_simp(g, quiet=quiet, stats=stats); record_stats('gadget_simp')
        spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
        to_gh(g)
        while True:
            i1 = id_simp(g, quiet=quiet, stats=stats); record_stats('id_simp')
            i2 = spider_simp(g, quiet=quiet, stats=stats); record_stats('spider_simp')
            i3 = pivot_simp(g, quiet=quiet, stats=stats); record_stats('pivot_simp')
            i4 = lcomp_simp(g, quiet=quiet, stats=stats); record_stats('lcomp_simp')
            if i1+i2+i3+i4==0: break
        j = pivot_gadget_simp(g, quiet=quiet, stats=stats); record_stats('pivot_gadget_simp')
        if i+j == 0:
            break
    return record

def plot_record(record: List[Dict[str, Union[int, str]]], png_file: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    df = pd.DataFrame(record)
    if df.empty: return
    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    df['method_cate'] = pd.Categorical(df['method'])
    df_2q = df[df['num_2q'] > 0]
    df_t = df[df['num_t'] > 0]
    for i, (df_, title) in enumerate([(df_2q, '2qubit gates'), (df_t, 'T gates')]):
        ax[i].scatter(df_['i_step'], df_['num_2q' if title == '2qubit gates' else 'num_t'], c=df_['method_cate'].cat.codes, cmap='tab20')
        ax[i].set_title(title)
        ax[i].set_xlabel('Step')
        ax[i].set_ylabel('Number of gates')
        ax[i].grid()
        handles, labels = ax[i].get_legend_handles_labels()
        ax[i].legend(handles, labels, title='Method', loc='upper right')
    plt.savefig(png_file)

def run_opt(qasm_file: str, q: Queue):
    import time
    from pyzx import full_reduce
    try:
        print(f"Processing {qasm_file}")
        time_init = time.time()
        circ_init = Circuit.from_qasm_file(qasm_file).to_basic_gates()
        graph = circ_init.to_graph()
        record = full_reduce_step(graph)
        plot_record(record, f"plots/{os.path.basename(qasm_file).split('.')[0]}.png")
        circ_final = extract_circuit(graph, up_to_perm=UP_TO_PERM).to_basic_gates()
        num_2q_init, num_t_init = num_2q_gates(circ_init), num_t_gates(circ_init)
        num_2q_final, num_t_final = num_2q_gates(circ_final), num_t_gates(circ_final)
        q.put({
            'bmark': os.path.basename(qasm_file).split('.')[0],
            'num_2q_init': num_2q_init,
            'num_t_init': num_t_init,
            'num_2q_final': num_2q_final,
            'num_t_final': num_t_final,
            'time': time.time() - time_init
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        q.put({'error': str(e), 'bmark': os.path.basename(qasm_file).split('.')[0]})

def run_with_timeout(qasm_file, timeout):
    q = Queue()
    process = Process(target=run_opt, args=(qasm_file, q))
    process.start()
    process.join(timeout)
    if process.is_alive():
        process.terminate()  # Kill the process if it exceeds the timeout
        process.join()  # Ensure the process is fully cleaned up
        return {'bmark': os.path.basename(qasm_file).split('.')[0], 'error': 'Timeout'}
    else:
        return q.get()  # Retrieve the result from the queue

datas = {}
for qasm in QASMS:
    result = run_with_timeout(qasm, timeout=TIMEOUT)
    if "error" not in result:
        datas[result["bmark"]] = result
    else:
        print(f"Error for {result['bmark']}: {result['error']}")
pd.DataFrame(datas).T.to_csv('results.csv')
