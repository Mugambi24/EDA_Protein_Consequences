#!/usr/bin/env python3
"""Reproduce planted-phase calibration and check a naive local-link baseline."""
import sys,json,math,hashlib,statistics
from pathlib import Path
R=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(R/'authoring/provenance/independent_submission'))
from inference import reconcile_phase_graph

def check(kind):
    assay=json.loads((R/f'tests/{kind}_phase_graph_assay.json').read_text())
    truth={row['case_id']:row['phase_bits'] for row in json.loads((R/f'tests/{kind}_phase_graph_truth.json').read_text())['records']}
    computed=reconcile_phase_graph(assay)
    abs_error=[];brier=[];nulls=[];confident=[];truth_bits=[]
    min_confident_opposite=1000000
    min_case_opposite_count=1000000
    trivial_zero_cases=[]
    example_case=None
    for rec in computed['records']:
        bits=truth[rec['case_id']]
        confident_opposite=sum(s['p_opposite_anchor']>.5 for s in rec['loci'][1:])
        min_confident_opposite=min(min_confident_opposite,confident_opposite)
        min_case_opposite_count=min(min_case_opposite_count,sum(bits))
        if all(s['p_opposite_anchor']<=.05 for s in rec['loci']) and rec['map_posterior']>=.975:
            trivial_zero_cases.append(rec['case_id'])
        if rec['case_id']=='h09':
            example_case={'case_id':'h09', 'planted_opposite_bits':sum(bits), 'posterior_sites_preferring_opposite_phase':confident_opposite,
                          'map_posterior':rec['map_posterior']}
        for s,t in zip(rec['loci'],bits):
            p=s['p_opposite_anchor']
            abs_error.append(abs(p-t));brier.append((p-t)**2)
            nulls.append(s['phase_bit'] is None)
            if s['phase_bit'] is not None:confident.append(s['phase_bit']==t)
            truth_bits.append(t)
    return {'dataset':kind,'cases':len(computed['records']),'loci':len(truth_bits),
       'fraction_truth_phase_1':sum(truth_bits)/len(truth_bits),
       'fixed_all_zero_error_fraction':sum(truth_bits)/len(truth_bits),
       'calibrated_mean_absolute_error_against_planted_phase':round(statistics.mean(abs_error),6),
       'calibrated_phase_brier_score':round(statistics.mean(brier),6),
       'confident_call_accuracy':round(statistics.mean(confident),6),
       'uncertain_site_count':sum(nulls),
       'partial_case_count':sum(rec['quality_flag']=='partial' for rec in computed['records']),
       'minimum_planted_opposite_bits_per_case':min_case_opposite_count,
       'minimum_sites_preferring_opposite_phase_per_case':min_confident_opposite,
       'casewise_near_certain_all_zero_shortcuts':trivial_zero_cases,
       'named_hidden_case_check':example_case,
       'map_posterior_range':[min(row['map_posterior'] for row in computed['records']), max(row['map_posterior'] for row in computed['records'])]}

def main():
    results=[check(k) for k in ('sealed','hidden')]
    assert all(r['calibrated_phase_brier_score']<.08 for r in results)
    assert all(r['fixed_all_zero_error_fraction']>.28 for r in results)
    assert all(r['uncertain_site_count']>0 for r in results)
    assert all(not r['casewise_near_certain_all_zero_shortcuts'] for r in results)
    assert all(r['minimum_sites_preferring_opposite_phase_per_case']>=1 for r in results)
    record={'model':'calibrated signed-cycle phasing','status':'fully synthetic',
         'method':'independent bitmask enumeration of all anchored configurations with batch-specific control-fitted binomial link likelihoods',
         'independent_generator':'authoring/provenance/generate_phase_graph.py',
         'results':results,'untrusted_data_warning':'Planted truth is verifier-only, never agent environment.',
         'sha256':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'tests/sealed_phase_graph_assay.json',R/'tests/hidden_phase_graph_assay.json',R/'tests/sealed_phase_graph_truth.json',R/'tests/hidden_phase_graph_truth.json']}}
    dest=R/'authoring/evidence/phase_graph_validation_2026-11-09.json';dest.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps(results,indent=2))
if __name__=='__main__':main()
