#!/usr/bin/env python3
"""Author-only method-neutral audit: EM versus independently coded constrained LS.

Both methods are evaluated against identical independently generated ground truth,
not by matching one another or by comparing output to the reference output.
"""
from pathlib import Path
import importlib.util, json, hashlib, tempfile, contextlib, io, statistics
R=Path(__file__).resolve().parents[2]
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m
G=module(R/'authoring/provenance/generate_longitudinal_assay.py','g_long_audit')
solvers={
 'multinomial_EM':R/'solution/inference.py',
 'simplex_least_squares':R/'authoring/provenance/independent_submission/inference.py',
}
methods={k:module(v,'long_'+k) for k,v in solvers.items()}
seed_offsets=[0,100,251,419,731,1107,1763,2579,4207,6319,8191]
requirements={'case_mae':('.06','lt'),'case_p90':('.105','lt'),'donor_mae':('.115','lt'),
              'case_coverage':('.80','ge'),'donor_coverage':('.75','ge'),
              'case_mean_width':('.31','lt'),'donor_mean_width':('.48','lt')}

def score(out,truth):
    case_err=[];donor_err=[];case_cover=[];donor_cover=[];case_w=[];donor_w=[]
    for row in out['records']:
        t=truth[row['case_id']];est=row['mean_treatment_time_effect'];target=t['mean_treatment_time_effect']
        case_err.append(abs(est-target));case_cover.append(row['interval_95'][0]<=target<=row['interval_95'][1]);case_w.append(row['interval_95'][1]-row['interval_95'][0]);dtrue={d['donor_id']:d['treatment_time_effect'] for d in t['donors']}
        for donor in row['donors']:
            v=donor['treatment_time_effect'];x=dtrue[donor['donor_id']];donor_err.append(abs(v-x));donor_cover.append(donor['interval_95'][0]<=x<=donor['interval_95'][1]);donor_w.append(donor['interval_95'][1]-donor['interval_95'][0])
    return {'case_mae':statistics.mean(case_err),
            'case_p90':sorted(case_err)[int(.9*(len(case_err)-1))],
            'donor_mae':statistics.mean(donor_err),
            'case_coverage':statistics.mean(case_cover),
            'donor_coverage':statistics.mean(donor_cover),
            'case_mean_width':statistics.mean(case_w),
            'donor_mean_width':statistics.mean(donor_w)}
report=[]
with tempfile.TemporaryDirectory(prefix='independent_long_audit_') as temp:
    d=Path(temp);(d/'tests').mkdir()
    for tag in ('sealed','hidden'):
        (d/'tests'/f'{tag}_cases.json').write_bytes((R/'tests'/f'{tag}_cases.json').read_bytes())
    G.ROOT=d
    for seed in seed_offsets:
        G.SEED_OFFSET=seed
        for tag in ('sealed','hidden'):
            if seed==0:
                assay=json.loads((R/'tests'/f'{tag}_longitudinal_assay.json').read_text())
                truth={x['case_id']:x for x in json.loads((R/'tests'/f'{tag}_longitudinal_truth.json').read_text())['records']}
            else:
                with contextlib.redirect_stdout(io.StringIO()):G.generate(tag)
                assay=json.loads((d/'tests'/f'{tag}_longitudinal_assay.json').read_text())
                truth={x['case_id']:x for x in json.loads((d/'tests'/f'{tag}_longitudinal_truth.json').read_text())['records']}
            phased=json.loads((R/'tests'/f'{tag}_phased_assay.json').read_text())
            for method,impl in methods.items():
                scores=score(impl.reconcile_longitudinal(assay,phased),truth)
                passed=all((scores[k] < float(bound) if op=='lt' else scores[k] >= float(bound)) for k,(bound,op) in requirements.items())
                report.append({'estimator':method,'fixture':tag,'seed_offset':seed,'passes':passed,'metrics':{k:round(v,6) for k,v in scores.items()}})
summary={}
for method in methods:
    rec=[r for r in report if r['estimator']==method]
    summary[method]={'datasets':len(rec),'passed':sum(r['passes'] for r in rec),
      'max_case_mae':max(r['metrics']['case_mae'] for r in rec),
      'max_case_p90':max(r['metrics']['case_p90'] for r in rec),
      'max_donor_mae':max(r['metrics']['donor_mae'] for r in rec),
      'min_case_coverage':min(r['metrics']['case_coverage'] for r in rec),
      'min_donor_coverage':min(r['metrics']['donor_coverage'] for r in rec),
      'max_case_mean_width':max(r['metrics']['case_mean_width'] for r in rec),
      'max_donor_mean_width':max(r['metrics']['donor_mean_width'] for r in rec)}
doc={'purpose':'Assess method-neutral longitudinal pass bounds using two mathematically distinct correct estimators and independently seeded synthetic truth',
 'dataset_type':'Synthetic only; generated truth withheld from solver inputs',
 'solver_methods':{'multinomial_EM':'Likelihood-based 9-read-class emission deconvolution and capture correction',
                   'simplex_least_squares':'Active-set constrained Euclidean least squares on 4-product simplex followed by batch capture correction'},
 'calibration_scope':'Both estimate batch chemistry solely from provided pure and mixed independent standards; neither accesses test truth.',
 'bounds':{k:{'threshold':float(x),'comparison':y} for k,(x,y) in requirements.items()},
 'seed_offsets':seed_offsets,'summary':summary,'results':report,
 'methodology':'Each score is evaluated against independent latent donor and case truth; results are not matched to reference-estimator coefficients. Coverage/width gates are checked together with error gates. These are synthetic sensitivity results, not proof of calibration in real patients.'}
path=R/'authoring/evidence/longitudinal_two_method_neutrality.json';path.write_text(json.dumps(doc,indent=2,sort_keys=True)+'\n')
print(json.dumps(summary,indent=2))
assert all(row['passes'] for row in report),[(r['estimator'],r['fixture'],r['seed_offset']) for r in report if not r['passes']]
