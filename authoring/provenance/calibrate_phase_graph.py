#!/usr/bin/env python3
"""Reproducible, independent-method numeric tolerance and Brier calibration.

All data are synthetic.  This script reads sealed truth only in the authoring tree.
It does not modify any sealed/hidden tests or agent-accessible files.
"""
import hashlib
import importlib.util
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'authoring/provenance'))
from generate_phase_graph import build


def module(filename, name):
    spec = importlib.util.spec_from_file_location(name, filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

reference = module(ROOT/'solution/inference.py','reference_phase_calibration')
alternative = module(ROOT/'authoring/provenance/independent_submission/inference.py','independent_phase_calibration')


def measure(assay, truth):
    truth_by_id = {r['case_id']:r['phase_bits'] for r in truth['records']}
    ref = reference.reconcile_phase_graph(assay)['records']
    alt = {r['case_id']:r for r in alternative.reconcile_phase_graph(assay)['records']}
    sites = []
    map_differences=[]
    site_differences=[]
    for row in ref:
        row_alt = alt[row['case_id']]
        map_differences.append(abs(row['map_posterior']-row_alt['map_posterior']))
        for site, alt_site, planted_bit in zip(row['loci'], row_alt['loci'], truth_by_id[row['case_id']]):
            sites.append((site['p_opposite_anchor'],planted_bit))
            site_differences.append(abs(site['p_opposite_anchor']-alt_site['p_opposite_anchor']))
    return {'cases':len(ref),'sites':len(sites),
        'mean_brier':round(statistics.mean((p-bit)**2 for p,bit in sites),9),
        'fixed_zero_brier':round(statistics.mean(bit for _,bit in sites),9),
        'fixed_half_brier':0.25,
        'max_cross_method_site_difference':round(max(site_differences),9),
        'max_cross_method_map_difference':round(max(map_differences),9)}


def main():
    existing={}
    inputs=[]
    for kind in ('sealed','hidden'):
        assay_path=ROOT/f'tests/{kind}_phase_graph_assay.json'
        truth_path=ROOT/f'tests/{kind}_phase_graph_truth.json'
        assay=json.loads(assay_path.read_text());truth=json.loads(truth_path.read_text())
        existing[kind]=measure(assay,truth)
        inputs.extend([assay_path,truth_path])
    simulations=[]
    for seed in range(1200,1232):
        assay,truth=build(seed,[f'calibration_{i:02d}' for i in range(16)])
        m=measure(assay,truth)
        simulations.append({'seed':seed,**m})
    assert all(x['mean_brier']<.08 for x in list(existing.values())+simulations)
    assert all(x['fixed_zero_brier']>.08 and x['fixed_half_brier']>.08 for x in list(existing.values()))
    assert all(x['max_cross_method_site_difference']<=.000001 and x['max_cross_method_map_difference']<=.000001 for x in list(existing.values())+simulations)
    record={
      'description':'Fully synthetic method-neutral phase-graph accuracy calibration; no patient evidence',
      'calibration_algorithm_a':'solution/inference.py exact tuple-enumeration log-sum-exp',
      'calibration_algorithm_b':'authoring/provenance/independent_submission/inference.py independent integer-bitmask enumeration',
      'requirements':{'per_site_absolute_posterior_error_strict':.03,'per_case_map_posterior_error_strict':.025,'mean_absolute_posterior_error_strict':.01,'mean_brier_strict':.08},
      'definition':'Brier score is mean squared error (reported posterior opposite-anchor phase probability - planted binary phase bit)^2 over all sites in each withheld dataset',
      'six_decimal_rounding':True,
      'sealed_hidden':existing,
      'simulation':{'seeds':'1200-1231 inclusive','cases_per_seed':16,'datasets':32,
       'maximum_brier':max(x['mean_brier'] for x in simulations),
       'mean_brier':round(statistics.mean(x['mean_brier'] for x in simulations),9),
       'datasets_passing_brier':sum(x['mean_brier']<.08 for x in simulations),
       'results':simulations},
      'interpretation':'The 0.03 marginal and 0.025 MAP tolerances are loose agreement ceilings, not measured computational error. The 0.08 Brier gate is supported by simulations and separates fixed guesses, not clinical calibration.',
      'sha256':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
    }
    out=ROOT/'authoring/evidence/phase_graph_calibration_2026-11-09.json'
    out.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'sealed_hidden':existing,'simulated_pass':record['simulation']['datasets_passing_brier'],'simulated_max_brier':record['simulation']['maximum_brier']},indent=2))
if __name__=='__main__':main()
