#!/usr/bin/env python3
"""Independent-truth calibration evidence for the batch assay (synthetic only).

Scores a multinomial-likelihood reference and a separate constrained-LS
implementation on planted fractions, not on each other's predictions.
"""
import importlib.util, json, math, statistics, copy
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
G=load(R/'authoring/provenance/generate_batch_assay.py','batch_generator')
F={
  'reference':load(R/'solution/inference.py','batch_ref').reconcile_batches,
  'alternative':load(R/'authoring/provenance/independent_submission/inference.py','batch_alt').reconcile_batches,
}
records=[]
for tag in ('sealed','hidden'):
    case_ids=[v['case_id'] for v in json.loads((R/'tests'/f'{tag}_cases.json').read_text())['cases']]
    base_assay=json.loads((R/'tests'/f'{tag}_batch_assay.json').read_text())
    base_truth=json.loads((R/'tests'/f'{tag}_batch_truth.json').read_text())['records']
    for seed in (None,13021+501,13021+1013,13021+2027,19021+3337,19021+6061):
        if seed is None: assay,truth=base_assay,base_truth
        else: assay,truth=G.make(case_ids,seed,'EDA-batch-assay-v1')
        target={r['case_id']:r for r in truth}
        for name,method in F.items():
            out=method(assay)
            errors=[];covered=[];widths=[]
            for r in out['records']:
                for group in ('case','control'):
                    v=r[group];x=target[r['case_id']][group+'_a1_fraction']
                    assert v['quality_flag']=='estimated'
                    errors.append(abs(v['a1_fraction']-x))
                    covered.append(v['interval_95'][0]<=x<=v['interval_95'][1])
                    widths.append(v['interval_95'][1]-v['interval_95'][0])
            metric={
                'max_error':max(errors),'mean_error':statistics.mean(errors),
                'coverage':statistics.mean(covered),'max_width':max(widths),
                'mean_width':statistics.mean(widths)}
            # Same pre-existing published scientific bounds.
            ok=(metric['max_error']<.065 and metric['mean_error']<=.035
                and sum(covered)>=math.ceil(.85*len(errors))
                and metric['max_width']<.45 and metric['mean_width']<.32)
            records.append({'tag':tag,'seed':'fixture' if seed is None else seed,'method':name,
                            'passed':ok, 'metrics':{k:round(v,6) for k,v in metric.items()}})
summary={}
for method in F:
    sub=[r for r in records if r['method']==method]
    summary[method]={'datasets':len(sub),'passed':sum(r['passed'] for r in sub),
       'worst_max_error':max(r['metrics']['max_error'] for r in sub),
       'worst_mean_error':max(r['metrics']['mean_error'] for r in sub),
       'minimum_coverage':min(r['metrics']['coverage'] for r in sub),
       'largest_max_width':max(r['metrics']['max_width'] for r in sub),
       'largest_mean_width':max(r['metrics']['mean_width'] for r in sub)}
# Quantify existing response-test perturbations, measuring actual output change.
assay=json.loads((R/'tests/sealed_batch_assay.json').read_text())
base=F['reference'](assay)['records']; changed=copy.deepcopy(assay)
for row in changed['calibrators']:
    if row['batch_id']=='lane_B':
        c=row['counts'];shift=min(35,c['A2_unique']//3)
        c['A1_unique']+=shift;c['A2_unique']-=shift
updated=F['reference'](changed)['records']
diffs=[abs(a['case']['a1_fraction']-b['case']['a1_fraction']) for a,b in zip(base,updated)]
response={'lane_B_spikein_response_max':round(max(diffs),6), 'published_minimum':.004}
# Independently audit every dataset case against a deliberately uncalibrated
# null method. This does not grade the alternative estimator against the
# reference; it compares raw observed counts with planted fractions.
casewise={}
for label in ('sealed','hidden'):
    assay=json.loads((R/'tests'/f'{label}_batch_assay.json').read_text())
    truth=json.loads((R/'tests'/f'{label}_batch_truth.json').read_text())['records']
    counts={}
    for row in assay['measurements']:
        pair=counts.setdefault((row['case_id'],row['group']),[0,0])
        pair[0]+=row['counts']['A1_unique']
        pair[1]+=row['counts']['A2_unique']
    records_case=[]
    for row in truth:
        errs=[]
        for group in ('case','control'):
            a,b=counts[row['case_id'],group]
            errs.append(abs(a/(a+b)-row[group+'_a1_fraction']))
        records_case.append({'case_id':row['case_id'],
                             'max_abs_raw_fraction_error':round(max(errs),6),
                             'raw_pm_0_11_can_pass_both_panels':all(e<=.11 for e in errs)})
    casewise[label]={
        'case_count':len(records_case),
        'min_case_max_raw_error':round(min(r['max_abs_raw_fraction_error'] for r in records_case),6),
        'raw_pooling_solved_cases':sum(r['raw_pm_0_11_can_pass_both_panels'] for r in records_case),
        'per_case':records_case,
    }
    assert casewise[label]['raw_pooling_solved_cases']==0
report={'scope':'Synthetic planted truth only. No experimental samples or clinical claims.',
        'provenance':'generator=generate_batch_assay.make; truth independently drawn and not provided to agent',
        'reference':'standard-library nonlinear robust calibration',
        'alternative':'separately implemented batch estimator evaluated against planted truth',
        'existing_bounds':{'max_error_lt':.065,'mean_error_le':.035,'coverage_ge':.85,
                           'max_width_lt':.45,'mean_width_lt':.32},
        'summary':summary,'perturbation_check':response,'records':records,
        'casewise_raw_pooling_discrimination':{'method':'Uncalibrated A1_unique/(A1_unique+A2_unique) with +/-0.11 intervals',
                                                'per_dataset':casewise}}
p=R/'authoring/evidence/batch_assay_calibration_2026-11-09.json'
p.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'summary':summary,'response':response},indent=2))
if not all(r['passed'] for r in records):
    print('Failed seeds:',[(r['tag'],r['seed'],r['method'],r['metrics']) for r in records if not r['passed']])
    raise SystemExit(1)
