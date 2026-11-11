#!/usr/bin/env python3
"""Reproduce verifier perturbation magnitudes using distinct synthetic estimators.

Does not import the verifier, alter agent inputs, or use clinical data. Both
implementations receive identical fixture panels and perturbed copies.
"""
from pathlib import Path
import importlib.util, json, math, copy, hashlib, statistics
ROOT=Path(__file__).resolve().parents[2]
def load_module(filepath,name):
    spec=importlib.util.spec_from_file_location(name,str(filepath)); module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module);return module
def read(name):return json.loads((ROOT/name).read_text())
models={
    'reference_robust_interpolation_EM_tuple':load_module(ROOT/'solution/inference.py','perturb_reference'),
    'independent_regression_LS_bitmask':load_module(ROOT/'authoring/provenance/independent_submission/inference.py','perturb_alternative'),
}

def panels_by_id(result):return {r['case_id']:r for r in result['records']}

def calibrated_batch(mod):
    assay=read('tests/sealed_batch_assay.json');base=panels_by_id(mod.reconcile_batches(assay)); altered=copy.deepcopy(assay)
    for row in altered['calibrators']:
        if row['batch_id']=='lane_B':
            c=row['counts'];shift=min(35,c['A2_unique']//3)
            c['A1_unique']+=shift;c['A2_unique']-=shift
    changed=panels_by_id(mod.reconcile_batches(altered))
    spike=[abs(base[c]['case']['a1_fraction']-changed[c]['case']['a1_fraction']) for c in base]
    changed=copy.deepcopy(assay)
    for row in changed['calibrators']:
        if row['known_a1_fraction']==.5 and row['replicate']==0:
            c=row['counts'];c['A1_unique']+=c['A2_unique'];c['A2_unique']=0
    bad=panels_by_id(mod.reconcile_batches(changed))
    corrupt=[abs(base[c]['case']['a1_fraction']-bad[c]['case']['a1_fraction']) for c in base]
    return {'lane_b_calibration_max_response':max(spike), 'lane_b_calibration_mean_response':statistics.mean(spike),
      'contaminated_replicate_max_deviation':max(corrupt),'contaminated_replicate_mean_deviation':statistics.mean(corrupt)}

def calibrated_phased(mod):
    assay=read('tests/sealed_phased_assay.json');base=panels_by_id(mod.reconcile_phased_rna(assay))
    alt=copy.deepcopy(assay)
    for row in alt['standards']:
        if row['batch_id'] in ('lane_A','lane_B') and row['standard_id'] not in ('pure_0','pure_1','pure_2','pure_3'):
            fr=row['known_product_fractions'];scaled={k:v*(1.35 if k.endswith('A1') else .65) for k,v in fr.items()};tot=sum(scaled.values())
            row['known_product_fractions']={k:v/tot for k,v in scaled.items()}
    corrected=panels_by_id(mod.reconcile_phased_rna(alt))
    mix=[abs(base[c]['case']['hap1_a1_fraction']-corrected[c]['case']['hap1_a1_fraction']) for c in base]
    alt=copy.deepcopy(assay)
    for row in alt['measurements']:
        if row['case_id']=='m03' and row['group']=='case':
            v=min(45,row['counts']['h1_A2']);row['counts']['h1_A2']-=v;row['counts']['h1_A1']+=v
    changed=panels_by_id(mod.reconcile_phased_rna(alt))
    hap1=changed['m03']['case']['hap1_a1_fraction']-base['m03']['case']['hap1_a1_fraction']
    hap0=abs(changed['m03']['case']['hap0_a1_fraction']-base['m03']['case']['hap0_a1_fraction'])
    alt=copy.deepcopy(assay)
    row=next(r for r in alt['standards'] if r['batch_id']=='lane_A' and r['standard_id']=='a1_rich' and r['replicate_id']=='r1')
    counts=row['counts'];total=sum(counts.values());counts.update({k:0 for k in counts});counts['unphased_A2']=total
    changed=panels_by_id(mod.reconcile_phased_rna(alt))
    contam=[abs(base[c]['case']['hap0_a1_fraction']-changed[c]['case']['hap0_a1_fraction']) for c in base]
    return {'mixed_calibration_max_response':max(mix),'mixed_calibration_mean_response':statistics.mean(mix),
            'reassigned_junction_hap1_signed_response':hap1,'reassigned_junction_hap0_abs_response':hap0,
            'contaminated_standard_max_hap0_deviation':max(contam),'contaminated_standard_mean_hap0_deviation':statistics.mean(contam)}

def calibrated_graph(mod):
    assay=read('tests/sealed_phase_graph_assay.json')
    base=panels_by_id(mod.reconcile_phase_graph(assay))
    first=sorted(assay['cases'],key=lambda row:row['case_id'])[0]['case_id']
    alt=copy.deepcopy(assay)
    for row in alt['calibration']:
        if row['batch_id']=='lane_A': row['known_relation']='opposite' if row['known_relation']=='same' else 'same'
    changed=panels_by_id(mod.reconcile_phase_graph(alt))
    cal_delta=sum(abs(a['p_opposite_anchor']-b['p_opposite_anchor']) for a,b in zip(base[first]['loci'],changed[first]['loci']))
    alt=copy.deepcopy(assay)
    for edge in alt['cases'][0]['edges']:
        edge['counts']['cis'],edge['counts']['trans']=edge['counts']['trans'],edge['counts']['cis']
    changed=panels_by_id(mod.reconcile_phase_graph(alt))
    cid=alt['cases'][0]['case_id']
    edge_diffs=[abs(a['p_opposite_anchor']-b['p_opposite_anchor']) for a,b in zip(base[cid]['loci'],changed[cid]['loci'])]
    return {'reversed_lane_a_calibration_sum_response':cal_delta,'swapped_edge_counts_max_site_response':max(edge_diffs),
            'swapped_edge_counts_mean_site_response':statistics.mean(edge_diffs)}

LIMITS={
    'lane_b_calibration_max_response': {'operator':'>','threshold':.004},
    'contaminated_replicate_max_deviation':{'operator':'<','threshold':.03},
    'mixed_calibration_max_response':{'operator':'>','threshold':.004},
    'reassigned_junction_hap1_signed_response':{'operator':'>','threshold':.035},
    'contaminated_standard_max_hap0_deviation':{'operator':'<','threshold':.05},
    'reversed_lane_a_calibration_sum_response':{'operator':'>','threshold':.10},
    'swapped_edge_counts_max_site_response':{'operator':'>','threshold':.03},
}
records={}
for name,mod in models.items():
    m={};m.update(calibrated_batch(mod));m.update(calibrated_phased(mod));m.update(calibrated_graph(mod))
    records[name]={'metrics':{k:round(v,9) for k,v in m.items()},'gate_pass':{
        k:(v>LIMITS[k]['threshold'] if LIMITS[k]['operator']=='>' else v<LIMITS[k]['threshold']) for k,v in m.items() if k in LIMITS}}
summary={k:{'min_observed':round(min(r['metrics'][k] for r in records.values()),9),
             'max_observed':round(max(r['metrics'][k] for r in records.values()),9),
             **limits,'both_methods_pass':all(r['gate_pass'][k] for r in records.values())} for k,limits in LIMITS.items()}
report={'purpose':'Method-separated calibration of seven verifier-enforced response/discrimination inequalities plus one non-gating hap0 diagnostic',
  'scope':'synthetic sealed fixtures only; not clinical data, not an empirical proof of universal robustness',
  'methods':{'reference_robust_interpolation_EM_tuple':'reference batch robust interpolation, phased EM and tuple-based exact phase posterior',
             'independent_regression_LS_bitmask':'separately implemented batch quadratic regression, phased constrained least squares and bitmask exact phase posterior'},
  'protocol':'Perturbation operations reproduce those in tests/test_workflow.py, but numerical outputs are computed without importing the verifier. Two distinct correct estimators process the same synthetic assay inputs.',
  'negative_controls':{'frozen_outputs':'A constant or cached output produces zero response to perturbations and fails all positive-change gates',
     'nonselective_junction':'Hap0 response to the hap1 perturbation is recorded as a diagnostic only; the verifier enforces a per-case hap1 increase, not a hap0 drift cap',
     'unrobust_contamination':'The upper-bound gates detect instability under a single calibration replicate corrupted without changing true sample observations'},
  'interpretation':'Lower thresholds require an observable response exceeding six-decimal rounding noise, but retain separation below measured correct-model responses. Upper thresholds reject large sensitivity to isolated damaged standards while remaining above independent correct-model observations. These are chosen for the fixed specified synthetic perturbations; they are not general biological calibration estimates.',
  'metrics':summary,'non_gating_diagnostics':{'reassigned_junction_hap0_abs_response':{'graded':False,'verifier_assertion':False,'description':'Measured for descriptive analysis only; no 0.04 hap0 drift bound is enforced by the actual verifier.', 'min_observed':round(min(r['metrics']['reassigned_junction_hap0_abs_response'] for r in records.values()),9),'max_observed':round(max(r['metrics']['reassigned_junction_hap0_abs_response'] for r in records.values()),9)}},'methods_results':records,
  'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'tests/sealed_batch_assay.json', ROOT/'tests/sealed_phased_assay.json',ROOT/'tests/sealed_phase_graph_assay.json']},
  'reproduction_command':'python3 authoring/provenance/validate_perturbation_bounds.py'}
dest=ROOT/'authoring/evidence/perturbation_bounds_validation_2026-11-10.json';dest.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'metrics':summary,'methods_results':records},indent=2))
if not all(v['both_methods_pass'] for v in summary.values()):raise SystemExit('A current perturbation bound did not pass independent-method calibration.')
