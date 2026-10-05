import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from binderflow.bli import fit_arrays, response, read_traces
from binderflow.cli import run
from binderflow.design import load_config, prepare, design_spec, mpnn_input
from binderflow.io import Stage, save_json
from binderflow.scoring import validate_prediction, filter_reasons, identity
from binderflow.structure import ca, fit_transform, rmsd, read_pdb
from scripts.demo import setup

class PipelineTests(unittest.TestCase):

    def test_end_to_end_and_resume_and_tamper(self):
        with tempfile.TemporaryDirectory() as td:
            c = load_config(setup(td))
            run(c)
            out = Path(c['run']['output'])
            summary = json.loads((out / '06_results/summary.json').read_text())
            self.assertEqual(summary['unique_evaluated'], 4)
            self.assertEqual(summary['passed'], 2)
            self.assertEqual(summary['selected'], 2)
            log = out / '02_backbones/bb0000/rfdiffusion.log'
            before = log.stat().st_mtime_ns
            run(c)
            self.assertEqual(log.stat().st_mtime_ns, before)
            f = out / '02_backbones/bb0000/complex.pdb'
            f.write_text(f.read_text() + 'REMARK tampered\n')
            with self.assertRaisesRegex(RuntimeError, 'changed'):
                run(c)

    def test_anchors_and_chain_assignment(self):
        with tempfile.TemporaryDirectory() as td:
            c = load_config(setup(td))
            m = prepare(c)
            spec = design_spec(c, m, 0)
            self.assertEqual(spec['anchors'], {'4': {'donor_position': 4, 'aa': 'K'}, '9': {'donor_position': 9, 'aa': 'D'}})
            run(c, through='sequences')
            candidates = json.loads((Path(c['run']['output']) / 'candidates.json').read_text())
            for row in candidates:
                self.assertEqual((row['sequence'][3], row['sequence'][8]), ('K', 'D'))
            obj = mpnn_input(candidates[0]['backbone'], 'test')
            self.assertEqual(obj['seq_chain_B'], m['target_sequence'])

    def test_prediction_chain_mismatch_and_pae_shape(self):
        with tempfile.TemporaryDirectory() as td:
            c = load_config(setup(td))
            run(c)
            pred = next((Path(c['run']['output']) / '04_predictions').glob('*/models/*rank_001*.pdb'))
            js = pred.with_name(pred.name.replace('_unrelaxed_', '_scores_').replace('.pdb', '.json'))
            with self.assertRaisesRegex(ValueError, 'sequence'):
                validate_prediction(pred, js, 'AAAA', 'BBBB')
            rs = read_pdb(pred)
            from binderflow.structure import sequence
            b = sequence([r for r in rs if r.chain == 'A'])
            t = sequence([r for r in rs if r.chain == 'B'])
            data = json.loads(js.read_text())
            data['pae'] = [[0.0]]
            save_json(js, data)
            with self.assertRaisesRegex(ValueError, 'invalid size'):
                validate_prediction(pred, js, b, t)

    def test_rigid_transform(self):
        rng = np.random.default_rng(3)
        a = rng.normal(size=(15, 3))
        rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        b = a @ rot + [4, -2, 8]
        r, t = fit_transform(a, b)
        self.assertLess(rmsd(a @ r + t, b), 1e-10)

    def test_noisy_bli_parameter_recovery(self):
        rng = np.random.default_rng(5)
        traces = []
        for i, c in enumerate([5e-09, 2e-08, 8e-08, 3.2e-07]):
            t = np.arange(0, 601, 3, dtype=float)
            y = response(t, c, 180, 200000.0, 0.004, 1 + i * 0.1, 0.01) + rng.normal(0, 0.001, len(t))
            traces.append({'name': str(i), 't': t, 'y': y, 'c': c, 'end': 180})
        fit, _ = fit_arrays(traces)
        self.assertLess(abs(fit['KD_nM'] - 20), 0.5)
        self.assertLess(abs(fit['kon_M_inverse_s_inverse'] / 200000.0 - 1), 0.03)
        self.assertLess(abs(fit['koff_s_inverse'] / 0.004 - 1), 0.03)
        self.assertIsNotNone(fit['KD_95pct_local_CI_M'])

    def test_bli_continuity(self):
        x = response(np.array([179.999999, 180.0, 180.000001]), 2e-08, 180.0, 200000.0, 0.004, 1.0, 0.0)
        self.assertLess(np.max(abs(np.diff(x))), 1e-07)

    def test_identity(self):
        self.assertEqual(identity('AAAA', 'AAAA'), 1)
        self.assertEqual(identity('AAAA', 'AAAT'), 0.75)
        self.assertEqual(identity('AAAA', 'AAA'), 0.75)

    def test_changed_config_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            c = load_config(setup(td))
            run(c, through='prepare')
            c['design']['seed'] += 1
            with self.assertRaisesRegex(RuntimeError, 'changed'):
                run(c)
if __name__ == '__main__':
    unittest.main()
