#!/usr/bin/env python3
"""Reproducible current-bundle checks; writes author-only non-agent evidence."""
import ast, hashlib, json, math, re, statistics, tomllib, copy, importlib.util
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def read(path):return json.loads((R/path).read_text())
meta=tomllib.loads((R/'task.toml').read_text())
assert meta['task']['name']=='afterquery/eda-protein-consequences'
assert meta['metadata']['category']=='Science'
inst=(R/'instruction.md').read_text()
suffix='You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.'
assert len(inst)<=10000 and inst.endswith('\n\n'+suffix+'\n')
assert not re.search('kep'+'ler', '\n'.join(
 (R/f).read_text(errors='ignore') for f in ('README.md','instruction.md','task.toml')),re.I)
tests=[n.name for n in ast.parse((R/'tests/test_workflow.py').read_text()).body
       if isinstance(n,ast.FunctionDef) and n.name.startswith('test_')]
assert len(tests)==72
batch=read('authoring/evidence/batch_assay_calibration_2026-11-09.json')
longitudinal=read('authoring/evidence/longitudinal_two_method_neutrality.json')
coding=read('authoring/evidence/coding_key_provenance_2026-11-09.json')
assert all(v['passed']==v['datasets']==22 for v in longitudinal['summary'].values())
assert all(v['passed']==v['datasets']==12 for v in batch['summary'].values())
assert sum(v['independent_resolved_checks'] for v in coding.values())==92
cases={};zero_baselines={};raw_baselines={}
for label,file in [('practice','environment/data/practice.json'),('sealed','tests/sealed_cases.json'),('hidden','tests/hidden_cases.json')]:
  d=read(file)['cases'];assert all(set(c)=={'case_id','variants'} for c in d)
  cases[label]={'cases':len(d),'case_fields':['case_id','variants'],'answer_label_keys':0}
for label in ('sealed','hidden'):
  truth=read(f'tests/{label}_longitudinal_truth.json')['records']
  donors=[abs(d['treatment_time_effect']) for r in truth for d in r['donors']]
  percase=[abs(r['mean_treatment_time_effect']) for r in truth]
  assert min(donors)>.15 and min(percase)>.15
  zero_baselines[label]={'donors':len(donors),'cases':len(percase),
    'min_absolute_donor_truth':round(min(donors),6),
    'min_absolute_case_truth':round(min(percase),6),
    'constant_zero_donor_MAE':round(statistics.mean(donors),6),
    'donors_with_truth_inside_fixed_[-0.15,0.15]':sum(v<=.15 for v in donors)}
  readings=read(f'tests/{label}_longitudinal_assay.json')['measurements']
  by_donor={}
  for row in readings:
    c=row['counts'];h0=c['h0_A1']/(c['h0_A1']+c['h0_A2']);h1=c['h1_A1']/(c['h1_A1']+c['h1_A2'])
    by_donor.setdefault((row['case_id'],row['donor_id']),{})[(row['condition'],row['time'])]=h1-h0
  est={}
  for (cid,donor),panel in by_donor.items():
    delta=panel['treated','followup']-panel['treated','baseline']-panel['untreated','followup']+panel['untreated','baseline']
    est.setdefault(cid,[]).append(delta)
  by_true={t['case_id']:t['mean_treatment_time_effect'] for t in truth}
  err=[abs(statistics.mean(v)-by_true[cid]) for cid,v in est.items()]
  raw_baselines[label]={'mean_abs_error':round(statistics.mean(err),6),
       'p90_lower_order':round(sorted(err)[int(.9*(len(err)-1))],6),
       'fails_published_case_gate':statistics.mean(err)>=.06 or sorted(err)[int(.9*(len(err)-1))]>=.105}
  assert raw_baselines[label]['fails_published_case_gate']

# Reproducible case-level naive-baseline and counterfactual audit. Every
# fitted method is scored against independently generated planted truth;
# response probes are run on altered input, not against reference predictions.
def audit_module(path,name):
  spec=importlib.util.spec_from_file_location(name,R/path)
  mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
  return mod
methods={
  'EM':audit_module('solution/inference.py','author_reference'),
  'independent_LS':audit_module('authoring/provenance/independent_submission/inference.py','author_alternative'),
}
casewise={}
for tag in ('sealed','hidden'):
  batch_truth={r['case_id']:r for r in read(f'tests/{tag}_batch_truth.json')['records']}
  phased_truth={r['case_id']:r for r in read(f'tests/{tag}_phased_truth.json')['records']}
  long_truth={r['case_id']:r for r in read(f'tests/{tag}_longitudinal_truth.json')['records']}
  batch_assay=read(f'tests/{tag}_batch_assay.json');phased_assay=read(f'tests/{tag}_phased_assay.json')
  long_assay=read(f'tests/{tag}_longitudinal_assay.json')
  batch_copy=[c for c,t in batch_truth.items() if abs(t['case_a1_fraction']-t['control_a1_fraction'])<.065]
  phased_copy=[c for c,t in phased_truth.items() if statistics.mean(abs(t['case'][h+'_a1_fraction']-t['control'][h+'_a1_fraction']) for h in ('hap0','hap1'))<.04]
  quarter=[c for c,t in long_truth.items() if min(abs(t['mean_treatment_time_effect']-.25),abs(t['mean_treatment_time_effect']+.25))<.105]
  lv=[t['mean_treatment_time_effect'] for t in long_truth.values()]
  median=statistics.median(lv)
  median_point_only=[c for c,t in long_truth.items() if abs(t['mean_treatment_time_effect']-median)<.105]
  assert not (batch_copy or phased_copy or quarter)
  casewise[tag]={'cases':len(batch_truth),'copy_control_batch_point_passes':len(batch_copy),
    'copy_control_phased_point_passes':len(phased_copy),
    'fixed_positive_negative_quarter_longitudinal_point_passes':len(quarter),
    'global_median_longitudinal_point_only_passes':len(median_point_only),
    'global_median_longitudinal_casewise_joint_passes':0,
    'note':'A global median may be numerically near some case means, but it cannot respond to the explicitly graded case-level treatment-junction perturbation; every case is checked for a response.'}
  batch_changed=copy.deepcopy(batch_assay)
  for row in batch_changed['measurements']:
    if row['group']=='case':
      cnt=row['counts'];moved=min(cnt['A2_unique'],round(.30*(cnt['A1_unique']+cnt['A2_unique'])))
      cnt['A2_unique']-=moved;cnt['A1_unique']+=moved
  phased_changed=copy.deepcopy(phased_assay)
  for row in phased_changed['measurements']:
    if row['group']=='case':
      cnt=row['counts'];moved=min(80,cnt['h1_A2']);cnt['h1_A2']-=moved;cnt['h1_A1']+=moved
  long_changed=copy.deepcopy(long_assay)
  for row in long_changed['measurements']:
    if row['condition']=='treated' and row['time']=='followup':
      cnt=row['counts']
      for a,b in (('h1_A1','h1_A2'),('unphased_A1','unphased_A2')):cnt[a],cnt[b]=cnt[b],cnt[a]
  for name,method in methods.items():
    base={r['case_id']:r for r in method.reconcile_batches(batch_assay)['records']}
    changed={r['case_id']:r for r in method.reconcile_batches(batch_changed)['records']}
    db=[changed[c]['case']['a1_fraction']-base[c]['case']['a1_fraction'] for c in base]
    base={r['case_id']:r for r in method.reconcile_phased_rna(phased_assay)['records']}
    changed={r['case_id']:r for r in method.reconcile_phased_rna(phased_changed)['records']}
    dp=[changed[c]['case']['hap1_a1_fraction']-base[c]['case']['hap1_a1_fraction'] for c in base]
    base={r['case_id']:r for r in method.reconcile_longitudinal(long_assay,phased_assay)['records']}
    changed={r['case_id']:r for r in method.reconcile_longitudinal(long_changed,phased_assay)['records']}
    dl=[abs(changed[c]['mean_treatment_time_effect']-base[c]['mean_treatment_time_effect']) for c in base]
    assert len(db)==len(dp)==len(dl)==len(batch_truth)
    assert min(db)>.02 and min(dp)>.03 and min(dl)>.03
    casewise[tag][name+'_response']={'minimum_batch':round(min(db),6),
      'minimum_phased':round(min(dp),6), 'minimum_longitudinal':round(min(dl),6),
      'all_cases_respond':True}

m=read('authoring/provenance/source_manifest.json')
for parent in ('batch_assay','phased_longitudinal_discrimination_revision'):
  for f,digest in m[parent]['sha256'].items():
    assert hashlib.sha256((R/f).read_bytes()).hexdigest()==digest, f
for file,digest in m['current_variant_input_sha256'].items():
  assert hashlib.sha256((R/file).read_bytes()).hexdigest()==digest
record={
 'scope':'Local only; no claim about external validation stages',
 'date':'2026-11-10',
 'problems_addressed':['no consequence labels in input','donor-level zero baseline defeated',
    'expected coding-key source disclosed and independently reconstructed',
    'batch accuracy/coverage/width and perturbation response calibrated',
    'README includes designated Difficulty, Reference solution and Verification sections'],
 'test_functions':len(tests),
 'reference_solver_full_pytest':'72 named tests; run new verifier for current outcomes',
 'independent_solver_full_pytest':'72 named tests; run independent verifier for current outcomes',
 'variant_input_schema':cases,'donor_zero_baselines':zero_baselines,
 'uncorrected_longitudinal_baselines':raw_baselines,
 'longitudinal_simulation_summary':longitudinal['summary'],
 'batch_simulation_summary':batch['summary'],
 'batch_lane_B_perturbation':batch['perturbation_check'],
 'coding_key_provenance':coding,
 'source_manifest_hashes':'all declared batch, phased, longitudinal and variant-input hashes recomputed and validated',
 'casewise_baseline_discrimination':casewise,
 'README_required_sections':all(('## '+s) in (R/'README.md').read_text() for s in ('Difficulty','Reference solution','Verification')),
 'instruction_characters':len(inst), 'instruction_suffix_exact':True
}
p=R/'authoring/evidence/current_submission_validation_2026-11-10.json'
p.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'input_case_schema':cases,'donor_baselines':zero_baselines,'raw_baselines':raw_baselines,
                  'batch_summary':batch['summary'],'longitudinal_summary':longitudinal['summary']},indent=2))
