"""Independent, iteration-oriented audit of molecule evidence.

Separate from the reference solver: sorts group rows and uses itertools.groupby,
manual category counting, and the beta moment identities.
"""
from itertools import groupby
from math import sqrt
import json
from pathlib import Path


def independently_verify(data):
    observations=data['observations']
    output=[]
    for cid,case_rows_iter in groupby(sorted(observations,key=lambda r:r['case_id']),key=lambda r:r['case_id']):
        case_rows=list(case_rows_iter)
        group_reports={}
        for group in ['case','control']:
            rows=[row for row in case_rows if row['group']==group and row['mapping_quality']>=30 and row['base_quality']>=20]
            rows=sorted(rows,key=lambda r:(r['umi'],-r['mapping_quality'],-r['base_quality'],r['signal']),reverse=False)
            unique=[]
            for key,parts in groupby(rows,key=lambda r:r['umi']):
                best=max(list(parts),key=lambda r:(r['mapping_quality'],r['base_quality'],r['signal']))
                unique.append(best)
            count1=sum(x['signal']=='A1_unique' for x in unique)
            count2=sum(x['signal']=='A2_unique' for x in unique)
            shared=sum(x['signal']=='shared' for x in unique)
            n=count1+count2
            alpha=count1+1.5; beta=count2+1.5
            total=alpha+beta
            frac=alpha/total
            sd=sqrt((frac*(1-frac))/(total+1))
            enough=n>=12
            group_reports[group]={
                'counts':{'A1_unique':count1,'A2_unique':count2,'shared':shared},
                'deduplicated':len(unique),'usable_unique':n,
                'a1_fraction':round(frac,6) if enough else None,
                'interval_95':[round(max(0,frac-1.96*sd),6),round(min(1,frac+1.96*sd),6)] if enough else None}
        a=group_reports['case']['a1_fraction']
        b=group_reports['control']['a1_fraction']
        delta=None if a is None or b is None else round(a-b,6)
        label=('insufficient_evidence' if delta is None else 'A1_enriched' if delta>=0.2 else 'A2_enriched' if delta<=-0.2 else 'inconclusive')
        output.append({'case_id':cid,'case':group_reports['case'],'control':group_reports['control'],'delta_a1_fraction':delta,'evidence_class':label})
    return {'records':output}

if __name__=='__main__':
    ROOT=Path(__file__).resolve().parents[2]
    for source,target in [('environment/data/practice_observations.json','authoring/evidence/practice_evidence_expected.json'),
                          ('tests/sealed_observations.json','tests/sealed_evidence_expected.json'),
                          ('tests/hidden_observations.json','tests/hidden_evidence_expected.json')]:
        (ROOT/target).write_text(json.dumps(independently_verify(json.loads((ROOT/source).read_text())),indent=2,sort_keys=True)+'\n')
