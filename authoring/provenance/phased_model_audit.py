"""Reproduce numerical evidence with independent reference and alternative solvers.

Only authoring-side files are read. Truth is not installed in agent or verifier.
"""
import json, math, sys, importlib.util
from pathlib import Path
R=Path(__file__).resolve().parents[2]

def import_solver(path,label):
    spec=importlib.util.spec_from_file_location(label,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod.reconcile_phased_rna

SOLVERS={
  'multinomial_EM':import_solver(R/'solution/inference.py','em_solver'),
  'projected_least_squares':import_solver(R/'authoring/provenance/independent_submission/inference.py','ls_solver')
}

def eval_one(fn, tag):
    assay=json.loads((R/'tests'/f'{tag}_phased_assay.json').read_text())
    truth={r['case_id']:r for r in json.loads((R/'tests'/f'{tag}_phased_truth.json').read_text())['records']}
    output=fn(assay)
    errs=[];coverage=0;width=[];maxwidth=0
    per_hap={h:{'errors':[],'covered':0,'widths':[]} for h in ('hap0','hap1')}
    for row in output['records']:
        t=truth[row['case_id']]
        for group in ['case','control']:
            for hap in ['hap0','hap1']:
                val=row[group][f'{hap}_a1_fraction'];ref=t[group][f'{hap}_a1_fraction']
                lo,hi=row[group][f'{hap}_interval_95']
                errs.append(abs(val-ref));coverage+=int(lo<=ref<=hi)
                width.append(hi-lo)
                per_hap[hap]['errors'].append(abs(val-ref));per_hap[hap]['covered']+=int(lo<=ref<=hi)
                per_hap[hap]['widths'].append(hi-lo)
    return {'num_fraction_estimates':len(errs), 'mean_absolute_error':round(sum(errs)/len(errs),6),
      'max_absolute_error':round(max(errs),6), 'interval_95_covered':coverage,'interval_95_total':len(errs),
      'mean_interval_width':round(sum(width)/len(width),6), 'max_interval_width':round(max(width),6),
      'per_haplotype':{h:{'MAE':round(sum(v['errors'])/len(v['errors']),6),
        'maximum_error':round(max(v['errors']),6),
        'errors_ge_0_08':sum(e>=.08 for e in v['errors']),
        'interval_coverage':round(v['covered']/len(v['errors']),6),
        'mean_width':round(sum(v['widths'])/len(v['widths']),6),
        'max_width':round(max(v['widths']),6),
        'meets_all_gates':sum(v['errors'])/len(v['errors'])<=.04 and max(v['errors'])<.09 and sum(e>=.08 for e in v['errors'])<=1 and v['covered']>=math.ceil(.85*len(v['errors'])) and sum(v['widths'])/len(v['widths'])<.28 and max(v['widths'])<.45} for h,v in per_hap.items()}}

results={name:{tag:eval_one(fn,tag) for tag in ('sealed','hidden')} for name,fn in SOLVERS.items()}
for tag in ('sealed','hidden'):
    assay=json.loads((R/'tests'/f'{tag}_phased_assay.json').read_text())
    truth={r['case_id']:r for r in json.loads((R/'tests'/f'{tag}_phased_truth.json').read_text())['records']}
    grouped={}
    for r in assay['measurements']:
        out=grouped.setdefault((r['case_id'],r['group']),{'hap0':[0,0],'hap1':[0,0]})
        for hap,c in [('hap0','h0'),('hap1','h1')]:
            out[hap][0]+=r['counts'][c+'_A1'];out[hap][1]+=r['counts'][c+'_A2']
    naive=[];fixed=[]
    for (cid,group),row in grouped.items():
        for h in ('hap0','hap1'):
            target=truth[cid][group][h+'_a1_fraction']
            x=row[h];naive.append(abs(x[0]/sum(x)-target));fixed.append(abs(.5-target))
    results.setdefault('baseline',{})[tag]={
      'uncorrected_haplotype_junction_MAE':round(sum(naive)/len(naive),6),
      'constant_fraction_0_5_MAE':round(sum(fixed)/len(fixed),6)}
results['provenance']={'data_type':'seeded synthetic RNA',
 'generator':'authoring/provenance/generate_phased_assay.py',
 'reference_method':'multinomial EM, standard-derived capture ratios',
 'alternate_method':'projected least squares, standard-derived capture ratios',
 'note':'No clinical samples, EDA patient reads, or experimentally measured allele-specific splicing are claimed.'}
path=R/'authoring/evidence/phased_multivariate_audit.json'
path.write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
print(json.dumps(results,indent=2))
