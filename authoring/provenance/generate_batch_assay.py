"""Seeded synthetic assay for research-grade EDA-A1/A2 batch reconciliation.

These counts are not experimental observations, and the generator and truths must
remain outside the agent image. The agent sees only calibration and observations.
"""
import json, random, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
SIGNALS=('A1_unique','A2_unique','shared')
FRACTIONS=(0.0,0.2,0.5,0.8,1.0)

def counts_for(rng,n,f,batch):
    # Paired chemistry-specific response functions are monotone yet biased
    # relative to the truth.  A raw A1_unique ratio cannot be used as f.
    # The numerical responses are *not* supplied to the learner; they must
    # be calibrated from the public replicated standards.
    strength={'lane_A':0.6716, 'lane_B':0.9052, 'lane_C':0.803}[batch]
    relative=math.exp(strength)
    frac=f*relative/(1-f+f*relative) if 0<f<1 else f
    frac=max(.015,min(.985,frac))
    shared={'lane_A':.21,'lane_B':.33,'lane_C':.44}[batch]
    ps=max(.05,min(.7,shared+.045*(f-.5)))
    p1=(1-ps)*frac
    p2=(1-ps)*(1-frac)
    vals=[0,0,0]
    for _ in range(n):
        q=rng.random()
        vals[0 if q<p1 else 1 if q<p1+p2 else 2]+=1
    return dict(zip(SIGNALS,vals))

def make(case_ids,seed,aid):
    rng=random.Random(seed)
    calibrators=[]
    for batch in ('lane_A','lane_B','lane_C'):
        for f in FRACTIONS:
            for rep in range(3):
                counts=counts_for(rng,750+rep*37,f,batch)
                # A corrupted middle-standard aliquot: robust calibration must
                # not trust every raw replicate equally; the clean replicates
                # remain sufficient for recovery.
                if f==0.5 and rep==0:
                    counts['A1_unique'],counts['A2_unique']=counts['A2_unique'],counts['A1_unique']
                calibrators.append({'batch_id':batch,'standard_id':f'{batch}_{int(f*100):03}_{rep}',
                                    'known_a1_fraction':f,'replicate':rep,'counts':counts})
    measures=[];truth=[]
    for i,cid in enumerate(case_ids):
        baseline=0.13+0.56*rng.random()
        shift=rng.choice([-.24,-.13,-.03,.04,.16,.25])
        fc=max(.06,min(.94,baseline+shift)); f0=baseline
        # Stratified positive-control design: every case pairs at least one
        # non-neutral mixture with its baseline. This is specified at the
        # experiment-generation stage, not selected from a solver's errors.
        # Near-balanced controls remain present and still test appropriate
        # response to actual RNA evidence changes.
        if max(abs(fc-.5),abs(f0-.5))<.23:
            fc=max(.78,f0+.19) if shift>=0 else min(.22,f0-.19)
            fc=max(.06,min(.94,fc))
        # Every individual case/control pair is a genuine treatment contrast.
        # This is a prospective synthetic-design condition, not oracle-based grading.
        if abs(fc-f0)<.13:
            direction = 1 if shift >= 0 else -1
            candidate=f0+direction*.155
            if not (.06 <= candidate <= .94):candidate=f0-direction*.155
            fc=candidate
        assert abs(fc-f0)>=.13
        for group,f in [('case',fc),('control',f0)]:
            for batch in ('lane_A','lane_B','lane_C'):
                for rep in range(2):
                    n=580+((i*13+rep*17)%5)*70
                    measures.append({'case_id':cid,'group':group,'batch_id':batch,
                                    'replicate':rep,'counts':counts_for(rng,n,f,batch)})
        truth.append({'case_id':cid,'case_a1_fraction':round(fc,6),
                      'control_a1_fraction':round(f0,6)})
    rng.shuffle(calibrators); rng.shuffle(measures)
    return {'assay_id':aid,'synthetic':True,'calibrators':calibrators,'measurements':measures},truth

def save():
    specs=[('practice','environment/data/practice.json',8021),('sealed','tests/sealed_cases.json',13021),('hidden','tests/hidden_cases.json',19021)]
    for label,casfile,seed in specs:
        ids=[r['case_id'] for r in json.loads((ROOT/casfile).read_text())['cases']]
        obj,truth=make(ids,seed,'EDA-batch-assay-v1')
        target=ROOT/('environment/data/batch_assay.json' if label=='practice' else f'tests/{label}_batch_assay.json')
        target.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
        if label!='practice':
            (ROOT/f'tests/{label}_batch_truth.json').write_text(json.dumps({'records':truth},indent=2,sort_keys=True)+'\n')
        print(label,len(ids),'measurements',len(obj['measurements']),'calibrators',len(obj['calibrators']))

if __name__=='__main__':save()
