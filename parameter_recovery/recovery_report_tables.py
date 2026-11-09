"""Tables of the single-edge recovery for the parameter-recovery report sent with
the reviewer response (NTM_Revision/parameter_recovery_report.tex).

Reads (all in this repository):
  recovery_all_models/{flux,w1}/recovery_summary.csv   synthetic_parameter_recovery_all_models.m
  ../network_validation/results/single_edge_error_locked.csv   single_edge_error_locked.m
Writes recovery_all_models/report/:
  recovery_error_single.csv, recovery_error_pairs.csv   median recovery error per output,
      group (interior, edge, ood), model and parameter, as a fraction of the training range
      (log10 coordinates for gamma1 and lambda1; for ood sets the truth is replaced by the
      nearest bound), with the number of sets / recoveries, the number not identifiable
      (flat loss) and the number on a bound
  table_pairs_mlp.tex        body of the pairwise (joint) recovery table, MLP and GP
  table_locked_error.tex     body of the prediction-error table on the locked sets
Run:  python -I recovery_report_tables.py   (from parameter_recovery/)
"""
import csv
import math
import os
import statistics

here = os.path.dirname(os.path.abspath(__file__))
rec_dir = os.path.join(here, 'recovery_all_models')
out_dir = os.path.join(rec_dir, 'report')
os.makedirs(out_dir, exist_ok=True)

RANGE = {'gamma1': (5e-5, 5e-2, True), 'lambda1': (0.01, 0.1, True),
         'delta': (10.0, 100.0, False), 'epsilon': (10.0, 100.0, False)}
PARAMS = ['gamma1', 'lambda1', 'delta', 'epsilon']
MODELS = ['MLP', 'GP', 'SR', 'LR']
GROUPS = ['interior', 'edge', 'ood']
TEX = {'gamma1': r'$\gamma_1$', 'lambda1': r'$\lambda_1$', 'delta': r'$\delta$', 'epsilon': r'$\epsilon$'}


def error(p, true, rec):
    lo, hi, is_log = RANGE[p]
    true = min(max(true, lo), hi)
    if is_log:
        return abs(math.log10(rec) - math.log10(true)) / math.log10(hi / lo)
    return abs(rec - true) / (hi - lo)


def summarize(output, typ):
    rows = list(csv.DictReader(open(os.path.join(rec_dir, output, 'recovery_summary.csv'))))
    out = []
    for g in GROUPS:
        for m in MODELS:
            for p in PARAMS:
                sel = [r for r in rows if r['model'] == m and r['type'] == typ and r['param'] == p
                       and r['set_id'].startswith(g + '_')]
                flat = [r for r in sel if r['flat'] == '1' or r['recovered_value'] in ('', 'NaN')]
                ok = [r for r in sel if r not in flat]
                errs = [error(p, float(r['true_value']), float(r['recovered_value'])) for r in ok]
                out.append({'output': output, 'type': typ, 'group': g, 'model': m, 'param': p,
                            'n': len(sel), 'n_flat': len(flat),
                            'n_bound': sum(r['at_bound'] == '1' for r in ok),
                            'median_error': statistics.median(errs) if errs else float('nan'),
                            'max_error': max(errs) if errs else float('nan')})
    return out


def write_csv(name, rows):
    with open(os.path.join(out_dir, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


single = summarize('flux', 'single') + summarize('w1', 'single')
pairs = summarize('flux', 'joint') + summarize('w1', 'joint')
write_csv('recovery_error_single.csv', single)
write_csv('recovery_error_pairs.csv', pairs)


def cell(r):
    if r['n_flat'] == r['n']:
        return 'n.i.'
    s = f"{100 * r['median_error']:.1f}"
    if r['n_bound'] > r['n'] / 2:
        s += r'\,b'
    return s


# pairwise table: rows output x model (MLP, GP), columns group x parameter
lines = []
for output, lab in [('flux', r'$\phi_\tau$'), ('w1', r'$\Delta_\tau$')]:
    for m in ['MLP', 'GP']:
        cells = []
        for g in GROUPS:
            for p in PARAMS:
                r = next(x for x in pairs if x['output'] == output and x['group'] == g and x['model'] == m and x['param'] == p)
                cells.append(cell(r))
        lines.append(f"{lab} & {m} & " + ' & '.join(cells) + r' \\')
    if output == 'flux':
        lines.append(r'\midrule')
open(os.path.join(out_dir, 'table_pairs_mlp.tex'), 'w').write('\n'.join(lines) + '\n')

# prediction error on the locked sets: median (max) over the 8 sets of each group
err = list(csv.DictReader(open(os.path.join(here, '..', 'network_validation', 'results', 'single_edge_error_locked.csv'))))


def fmt(x):
    if 0.1 <= x < 100:
        return f'{x:.2g}'
    e = math.floor(math.log10(x))
    return f'{x / 10 ** e:.1f}' + r'\times10^{' + str(e) + '}'


lines = []
for m, lab in [('mlp', 'MLP'), ('gp', 'GP'), ('sr', 'SR'), ('lr', 'LR')]:
    cells = []
    for target in ['flux_nmse', 'w1_nmse']:
        for g in GROUPS:
            v = [float(r[target]) for r in err if r['method'] == m and r['group'] == g]
            cells.append(f'${fmt(statistics.median(v))}$ (${fmt(max(v))}$)')
    lines.append(f'{lab} & ' + ' & '.join(cells) + r' \\')
open(os.path.join(out_dir, 'table_locked_error.tex'), 'w').write('\n'.join(lines) + '\n')

for r in single:
    if r['model'] == 'MLP':
        print(r['output'], r['group'], r['param'], f"{r['median_error']:.4f}", r['n_bound'], r['n_flat'])
print(open(os.path.join(out_dir, 'table_pairs_mlp.tex')).read())
print(open(os.path.join(out_dir, 'table_locked_error.tex')).read())
