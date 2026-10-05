from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .io import csv_write, save_json

def response(t, concentration_M, assoc_end, kon, koff, rmax, offset):
    rate = kon * concentration_M + koff
    req = rmax * kon * concentration_M / rate
    during = req * -np.expm1(-rate * np.minimum(t, assoc_end))
    return offset + during * np.exp(-koff * np.maximum(t - assoc_end, 0))

def read_traces(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        required = {'trace_id', 'time_s', 'response_nm', 'concentration_nM', 'association_end_s'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f'BLI CSV requires columns: {sorted(required)}')
        data = list(reader)
    grouped = {}
    for r in data:
        grouped.setdefault(r['trace_id'], []).append(r)
    traces = []
    for name, rows in grouped.items():
        times = np.array([float(r['time_s']) for r in rows])
        values = np.array([float(r['response_nm']) for r in rows])
        conc = np.array([float(r['concentration_nM']) for r in rows])
        ends = np.array([float(r['association_end_s']) for r in rows])
        if not all((np.isfinite(a).all() for a in [times, values, conc, ends])):
            raise ValueError(f'Nonfinite BLI values: {name}')
        if not np.all(conc == conc[0]) or not np.all(ends == ends[0]) or conc[0] <= 0 or (ends[0] <= 0):
            raise ValueError('Each trace needs one positive concentration and association duration; subtract blank/reference traces first')
        ix = np.argsort(times)
        t, y = (times[ix], values[ix])
        if len(t) < 12 or t[0] < 0 or np.any(np.diff(t) <= 0):
            raise ValueError(f'Trace {name}: need >=12 distinct nonnegative times')
        if not np.isclose(t[0], 0, atol=1e-06):
            raise ValueError(f'Trace {name}: first point must be t=0 at association start')
        if sum(t < ends[0]) < 5 or sum(t > ends[0]) < 5:
            raise ValueError(f'Trace {name}: insufficient association or dissociation points')
        traces.append({'name': name, 't': t, 'y': y, 'c': float(conc[0] * 1e-09), 'end': float(ends[0])})
    if len({r['c'] for r in traces}) < 3:
        raise ValueError('Need at least three distinct positive analyte concentrations for this global fitting workflow')
    return traces

def fit_arrays(traces, shared_rmax=False):
    n = len(traces)
    nr = 1 if shared_rmax else n
    scale = max(float(np.ptp(np.concatenate([r['y'] for r in traces]))), 0.0001)
    amp = np.array([max(float(np.ptp(r['y'])), 0.0001) for r in traces])
    offsets = np.array([r['y'][0] for r in traces])
    lo = np.r_[1.0, -7.0, np.full(nr, -8.0), np.full(n, -np.inf)]
    hi = np.r_[9.0, 0.0, np.full(nr, np.log10(scale * 1000)), np.full(n, np.inf)]

    def unpack(p):
        rm = 10 ** p[2:2 + nr]
        if shared_rmax:
            rm = np.repeat(rm, n)
        return (10 ** p[0], 10 ** p[1], rm, p[2 + nr:])

    def predict(p):
        kon, koff, rm, off = unpack(p)
        return [response(r['t'], r['c'], r['end'], kon, koff, rm[i], off[i]) for i, r in enumerate(traces)]

    def residual(p):
        return np.concatenate([pred - r['y'] for pred, r in zip(predict(p), traces)]) / scale
    best = None
    for lk in [4.0, 5.5, 7.0]:
        for lf in [-4.0, -2.0]:
            rm = np.array([amp.max() * 1.5]) if shared_rmax else amp * 1.5
            x = np.r_[lk, lf, np.log10(rm), offsets]
            result = least_squares(residual, x, bounds=(lo, hi), max_nfev=4000, xtol=1e-10, ftol=1e-10, gtol=1e-10)
            if best is None or result.cost < best.cost:
                best = result
    kon, koff, rm, off = unpack(best.x)
    pred = predict(best.x)
    res = np.concatenate([p - r['y'] for p, r in zip(pred, traces)])
    singular = np.linalg.svd(best.jac, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else None
    dof = len(res) - len(best.x)
    issues = []
    noise = float(np.sqrt(np.mean(res ** 2)))
    signal = max((float(np.ptp(v)) for v in pred))
    if signal < 5 * max(noise, 1e-12):
        issues.append('low_fitted_signal_relative_to_residual_noise')
    if not best.success:
        issues.append('optimizer_not_converged')
    if np.any(best.active_mask):
        issues.append('parameter_at_bound')
    if condition is None or condition > 100000000.0:
        issues.append('poor_parameter_identifiability')
    if koff * max((r['t'][-1] - r['end'] for r in traces)) < 0.1:
        issues.append('little_dissociation_observed')
    if max((r['c'] for r in traces)) < koff / kon or min((r['c'] for r in traces)) > koff / kon:
        issues.append('concentration_range_does_not_bracket_KD')
    ci = None
    if dof > 0 and condition is not None and (condition < 100000000.0) and (not np.any(best.active_mask)):
        cov = np.linalg.pinv(best.jac.T @ best.jac) * (2 * best.cost / dof)
        v = max(float(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]), 0.0)
        z = best.x[1] - best.x[0]
        if 1.96 * np.sqrt(v) < 10:
            ci = [float(10 ** (z - 1.96 * np.sqrt(v))), float(10 ** (z + 1.96 * np.sqrt(v)))]
    summary = {'kon_M_inverse_s_inverse': float(kon), 'koff_s_inverse': float(koff), 'KD_M': float(koff / kon), 'KD_nM': float(koff / kon * 1000000000.0), 'KD_95pct_local_CI_M': ci, 'rmse_nm': float(np.sqrt(np.mean(res ** 2))), 'optimizer_converged': bool(best.success), 'jacobian_condition_number': condition, 'issues': issues, 'interpretation': 'review_required' if issues else 'model_fit_completed_not_a_binding_verdict', 'rmax_mode': 'shared' if shared_rmax else 'per_trace', 'model': 'global 1:1 Langmuir; zero initial occupancy; per-trace constant offsets', 'uncertainty_note': 'Local Jacobian CI assumes independent homoscedastic errors; serial correlation and systematic assay errors are not covered.', 'traces': [{'trace_id': r['name'], 'Rmax_nm': float(rm[i]), 'offset_nm': float(off[i])} for i, r in enumerate(traces)]}
    return (summary, pred)

def fit_file(path, out, shared_rmax=False):
    traces = read_traces(path)
    summary, pred = fit_arrays(traces, shared_rmax)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    save_json(out / 'fit.json', summary)
    rows = [{'trace_id': r['name'], 'time_s': float(t), 'observed_nm': float(y), 'predicted_nm': float(v), 'residual_nm': float(y - v)} for r, p in zip(traces, pred) for t, y, v in zip(r['t'], r['y'], p)]
    csv_write(out / 'fitted_traces.csv', rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    for r, p in zip(traces, pred):
        line, = ax[0].plot(r['t'], p, label=f"{r['c'] * 1000000000.0:g} nM • {r['name']}")
        ax[0].scatter(r['t'], r['y'], s=7, color=line.get_color(), alpha=0.5)
        ax[1].plot(r['t'], r['y'] - p, color=line.get_color(), lw=0.8)
    ax[0].set(ylabel='Response (nm)', title=f"Global 1:1 fit • KD = {summary['KD_nM']:.3g} nM")
    ax[0].legend(fontsize=8)
    ax[1].axhline(0, color='black', lw=0.6)
    ax[1].set(xlabel='Time from association start (s)', ylabel='Residual (nm)')
    fig.tight_layout()
    fig.savefig(out / 'kinetics.png', dpi=180)
    plt.close(fig)
    return summary
