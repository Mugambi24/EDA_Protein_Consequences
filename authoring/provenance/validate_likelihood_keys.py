#!/usr/bin/env python3
"""Reproduce sealed likelihood expected keys from seeded observations and controls.

Reproducible key producer: tests/independent_calibration.py:calculate
 (derivative-based optimization).
Crosschecks: tests/independent_evidence.py (UMI sorting/grouping) and
 solution/inference.py:infer (independently implemented golden-section optimizer).
Neither the fixed answer keys nor the source audit are exposed to the agent.
All files are read from the authoring bundle; no external data access is used.
"""
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / 'tests'
sys.path.insert(0, str(T))
from independent_evidence import independently_verify
from independent_calibration import calculate


def load(rel):
    return json.loads((ROOT / rel).read_text())


def digest(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def named_module(path):
    spec = importlib.util.spec_from_file_location('likelihood_reference_oracle', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def numeric_differences(left, right):
    a = []
    for x, y in zip(left['records'], right['records']):
        assert x['case_id'] == y['case_id']
        assert x['mixture_class'] == y['mixture_class']
        for prop in ('delta_a1_mle',):
            va, vb = x[prop], y[prop]
            assert (va is None) == (vb is None)
            if va is not None: a.append(abs(va - vb))
        for group in ('case', 'control'):
            assert x[group]['molecules'] == y[group]['molecules']
            for prop in ('a1_mle', 'log_likelihood'):
                va, vb = x[group][prop], y[group][prop]
                assert (va is None) == (vb is None)
                if va is not None:a.append(abs(va - vb))
            va, vb = x[group]['profile_95'], y[group]['profile_95']
            assert (va is None) == (vb is None)
            if va is not None:a.extend(abs(i - j) for i, j in zip(va, vb))
    return a


def main():
    calibration = load('tests/assets/spikein_calibration.json')
    reference = named_module(ROOT / 'solution/inference.py')
    report = {'calibration_input': 'tests/assets/spikein_calibration.json',
              'reproducible_key_producer': 'tests/independent_calibration.py:calculate',
              'reference_crosscheck': 'solution/inference.py:infer',
              'independent_evidence': 'tests/independent_evidence.py:independently_verify',
              'independent_optimizer': 'tests/independent_calibration.py:calculate',
              'numeric_gate_absolute_tolerance': 0.000003,
              'output_rounding_decimal_places': 6,
              'seed_for_sensitivity': 20261109, 'fixtures': {}}
    rng = random.Random(report['seed_for_sensitivity'])
    for tag in ('sealed', 'hidden'):
        paths = {role: f'tests/{tag}_{role}.json' for role in
                 ('observations', 'evidence_expected', 'likelihood_expected')}
        evidence = independently_verify(load(paths['observations']))
        assert evidence == load(paths['evidence_expected']), f'{tag} evidence keys mismatch'
        expected = load(paths['likelihood_expected'])
        golden = reference.infer(evidence['records'], calibration)
        derivative = calculate(evidence['records'], calibration)
        assert derivative == expected, f'{tag} reproducible likelihood keys mismatch'
        independent_differences = numeric_differences(golden, expected)
        assert max(independent_differences, default=0) <= 0.000003
        sensitivity = []
        for original in evidence['records']:
            for _ in range(4):
                row = json.loads(json.dumps(original))
                for group in ('case', 'control'):
                    counts = row[group]['counts']
                    for signal in ('A1_unique', 'A2_unique', 'shared'):
                        counts[signal] = max(0, counts[signal] + rng.randint(-7, 7))
                trial_golden = reference.infer([row], calibration)
                trial_derivative = calculate([row], calibration)
                sensitivity.extend(numeric_differences(trial_golden, trial_derivative))
        max_delta = max(sensitivity, default=0)
        assert max_delta <= 0.000003, f'{tag} independent optimizer mismatch {max_delta}'
        report['fixtures'][tag] = {
            'source_files_sha256': {path: digest(path) for path in paths.values()},
            'records': len(expected['records']),
            'source_evidence_exact_match': True,
            'reproducible_keys_exact_match': True,
            'independent_optimizer_max_absolute_difference': max(independent_differences, default=0),
            'perturbed_panels': 4 * len(expected['records']) * 2,
            'perturbed_optimizer_max_absolute_difference': max_delta,
        }
    report['interpretation'] = (
        'The absolute 3e-6 threshold is three six-decimal output units. Fixed-key '
        'matches and deterministic seeded count perturbations compare two '
        'independently implemented convex optimizers (golden section versus '
        'derivative bisection). The bound permits minor numerical convergence '
        'differences without substituting biological accuracy criteria; all '
        'output-range, rounding, molecule-count and delta-identity constraints '
        'are enforced separately. This is algorithmic reproducibility evidence, '
        'not clinical measurement calibration.'
    )
    outfile = ROOT / 'authoring/evidence/likelihood_key_provenance_2026-11-09.json'
    outfile.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'fixtures': report['fixtures'],
                      'tolerance': report['numeric_gate_absolute_tolerance']}, indent=2))


if __name__ == '__main__':
    main()
