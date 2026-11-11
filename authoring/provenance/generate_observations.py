"""Independent seeded simulator for synthetic short-read junction-mapping evidence.

The classes are upstream aligner categories, not direct observations of binding
or transcriptional regulation in any individual. It is not clinical data.
"""
import json,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def synthesize(case_file, seed):
    rng=random.Random(seed)
    obj=json.loads(Path(case_file).read_text())
    observations=[]
    for i,case in enumerate(obj['cases']):
        cid=case['case_id']
        # Intentional underpowered panels belong only to the public practice data.
        # Sealed and held-out instances must require actual mixture estimation;
        # low-coverage null semantics are exercised separately in a dedicated test.
        low_coverage=(Path(case_file).name == 'practice.json' and i%11==0)
        for group in ('case','control'):
            n=(12 if low_coverage else (82+(i*7)%29))
            frac=(0.2+((i*13)%61)/100) if group=='case' else (0.33+((i*7)%29)/100)
            for j in range(n):
                p=rng.random()
                signal=('A1_unique' if p<frac*0.68 else 'A2_unique' if p<0.68 else 'shared')
                mq=rng.choice([10,18,28,30,35,42,50,60,60,60])
                bq=rng.choice([12,18,20,26,32,35,38,40])
                umi=f'{cid}_{group}_{j:04d}'
                row={'case_id':cid,'group':group,'umi':umi,'signal':signal,
                     'mapping_quality':mq,'base_quality':bq}
                observations.append(row)
                if j%9==0:  # PCR duplicates: sometimes more, sometimes less reliable
                    dup=dict(row)
                    dup['mapping_quality']=min(60,mq+3) if j%2==0 else max(0,mq-7)
                    observations.append(dup)
    rng.shuffle(observations)
    return {'observations':observations}

if __name__=='__main__':
    table=[('environment/data/practice.json','environment/data/practice_observations.json',100),
           ('tests/sealed_cases.json','tests/sealed_observations.json',200),
           ('tests/hidden_cases.json','tests/hidden_observations.json',300)]
    for a,b,s in table:
        (ROOT/b).write_text(json.dumps(synthesize(ROOT/a,s),separators=(',',':'))+'\n')
