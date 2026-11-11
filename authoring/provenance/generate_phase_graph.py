#!/usr/bin/env python3
"""Deterministic, synthetic noisy long-range haplotype-link graphs.

Truth bits never enter the learner-visible assay; calibration panels contain only
known cis/trans reference links. Read counts simulate DNA/RNA molecule linkage.
"""
import json,random,sys,math,hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
RATES={'lane_A':(.86,.19),'lane_B':(.39,.63),'lane_C':(.74,.43)}

def draw_counts(rng,n,p):
    k=sum(rng.random()<p for _ in range(n))
    return {'cis':k,'trans':n-k,'ambiguous':rng.randrange(0,max(3,n//4))}

def build(seed,ids):
    rng=random.Random(seed)
    controls=[]
    for batch,(cis,trans) in RATES.items():
        for relation,p in [('same',cis),('opposite',trans)]:
            for r in range(3):
                n=700+rng.randrange(0,150)
                controls.append({'batch_id':batch,'standard_id':f'{batch}_{relation}',
                     'replicate_id':f'r{r}','known_relation':relation,'counts':draw_counts(rng,n,p)})
    cases=[]; truth=[]
    for caseid in ids:
        size=rng.randrange(10,13)
        loci=[f'v{k:02d}' for k in range(size)]
        bits=[0]+[rng.randrange(2) for _ in range(size-1)]
        # A connected but noisy signed network, with loops that contradict naive
        # local majority assignment. Site labels contain no truth information.
        pairs={(j,rng.randrange(j)) for j in range(1,size)}
        while len(pairs)<2*size:
            j,k=rng.sample(range(size),2)
            pairs.add((max(j,k),min(j,k)))
        edges=[]
        for j,k in sorted(pairs):
            # Sparse low-depth sites require global inference across cycles.
            batch=rng.choice(list(RATES))
            n=rng.randrange(7,25)
            relation=bits[j]==bits[k]
            p=RATES[batch][0 if relation else 1]
            edges.append({'left':loci[j],'right':loci[k],'batch_id':batch,
                 'replicate_id':f'{caseid}_link_{j:02d}_{k:02d}', 'counts':draw_counts(rng,n,p)})
        # Avoid single-case zero-information shortcuts: a planted all-zero
        # configuration combined with overwhelming evidence for zero can be
        # accepted without solving a case. Use a separate deterministic stream
        # so the correction does not perturb other reproducible test panels.
        if not any(bits[1:]):
            digest=hashlib.sha256(f'{seed}:{caseid}:balanced-phase-v1'.encode()).hexdigest()
            repair_rng=random.Random(int(digest,16))
            while True:
                candidate=[0]+[repair_rng.randrange(2) for _ in range(size-1)]
                if size//3<=sum(candidate)<=size-1-size//3:
                    bits=candidate
                    break
            for edge in edges:
                j=loci.index(edge['left']);k=loci.index(edge['right'])
                depth=edge['counts']['cis']+edge['counts']['trans']
                phase_same=(bits[j]==bits[k])
                probability=RATES[edge['batch_id']][0 if phase_same else 1]
                new=draw_counts(repair_rng,depth,probability)
                # Keep non-informative molecule depth unchanged.
                edge['counts']['cis']=new['cis']
                edge['counts']['trans']=new['trans']
        rng.shuffle(edges)
        cases.append({'case_id':caseid,'loci':loci,'anchor':loci[0],'edges':edges})
        truth.append({'case_id':caseid,'phase_bits':bits})
    rng.shuffle(controls);rng.shuffle(cases)
    return {'assay_id':f'eda_linkgraph_{seed}','synthetic':True,'calibration':controls,'cases':cases}, {'assay_id':f'eda_linkgraph_{seed}','records':truth}

def main():
    sets=[('practice',2026110911,'environment/data/practice.json','environment/data/phase_graph_assay.json',None),
          ('sealed',2026110913,'tests/sealed_cases.json','tests/sealed_phase_graph_assay.json','tests/sealed_phase_graph_truth.json'),
          ('hidden',2026110917,'tests/hidden_cases.json','tests/hidden_phase_graph_assay.json','tests/hidden_phase_graph_truth.json')]
    for _,seed,casefile,dest,truthdest in sets:
        ids=[r['case_id'] for r in json.loads((ROOT/casefile).read_text())['cases']]
        assay,truth=build(seed,ids)
        (ROOT/dest).write_text(json.dumps(assay,indent=2,sort_keys=True)+'\n')
        if truthdest:(ROOT/truthdest).write_text(json.dumps(truth,indent=2,sort_keys=True)+'\n')
        print(dest,len(ids))
if __name__=='__main__':main()
