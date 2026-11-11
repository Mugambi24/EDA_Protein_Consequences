"""Aggregate independent junction-mapping observations after molecule-level QC.

No clinical pathogenicity interpretation is made. All read evidence in the task is
simulated, although protein identity is tied to an accessioned human sequence.
"""
from collections import defaultdict
from math import sqrt

ALLOWED = {'A1_unique', 'A2_unique', 'shared'}
GROUPS = ('case','control')

def summarize(data):
    observations = data['observations']
    by_case=defaultdict(list)
    for row in observations:
        if row.get('group') not in GROUPS or row.get('signal') not in ALLOWED:
            raise ValueError('Invalid group or junction support class')
        if type(row.get('mapping_quality')) is not int or type(row.get('base_quality')) is not int:
            raise ValueError('Invalid read quality data')
        if not isinstance(row.get('umi'),str) or not row['umi']:
            raise ValueError('Missing UMI')
        by_case[row['case_id']].append(row)
    reports=[]
    for case_id in sorted(by_case):
        groups={}
        for group in GROUPS:
            reads=[r for r in by_case[case_id] if r['group']==group]
            retained={}
            for r in reads:
                if r['mapping_quality'] < 30 or r['base_quality'] < 20:
                    continue
                key=r['umi']
                rank=(r['mapping_quality'],r['base_quality'],r['signal'])
                if key not in retained or rank>retained[key][0]:
                    retained[key]=(rank,r['signal'])
            counts={signal:sum(1 for _,sig in retained.values() if sig==signal) for signal in sorted(ALLOWED)}
            a=counts['A1_unique']+1.5
            b=counts['A2_unique']+1.5
            mean=a/(a+b)
            variance=a*b/((a+b)**2*(a+b+1))
            low=max(0,mean-1.96*sqrt(variance))
            high=min(1,mean+1.96*sqrt(variance))
            enough=(counts['A1_unique']+counts['A2_unique'])>=12
            groups[group]={'counts':counts,'deduplicated':len(retained),'usable_unique':counts['A1_unique']+counts['A2_unique'],
                           'a1_fraction':round(mean,6) if enough else None,
                           'interval_95':[round(low,6),round(high,6)] if enough else None}
        ca=groups['case']['a1_fraction'];co=groups['control']['a1_fraction']
        if ca is None or co is None:
            label='insufficient_evidence'; delta=None
        else:
            delta=round(ca-co,6)
            # The delta is a descriptive sample statistic, not a clinical significance test.
            label='A1_enriched' if delta>=0.2 else ('A2_enriched' if delta<=-0.2 else 'inconclusive')
        reports.append({'case_id':case_id,'case':groups['case'],'control':groups['control'],
                        'delta_a1_fraction':delta,'evidence_class':label})
    return reports
