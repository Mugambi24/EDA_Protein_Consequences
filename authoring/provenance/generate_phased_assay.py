"""Deterministic *synthetic* phased EDA isoform evidence simulation.

Model: four latent phased transcripts (hap0-A1, hap0-A2, hap1-A1,
hap1-A2) have lane-specific capture efficiencies and 9-category
allele/junction assignment errors. Standards include pure and mixed
composition, and one deliberately mislabeled/contaminated replicate.
Truth is generated separately and must never be copied into the agent image.
"""
import json,random,math,sys
from pathlib import Path
R=Path(__file__).resolve().parents[2]
CATS=['h0_A1','h0_A2','h1_A1','h1_A2','h0_shared','h1_shared','unphased_A1','unphased_A2','unassigned']
PRODS=['hap0_A1','hap0_A2','hap1_A1','hap1_A2']
LANES=['lane_A','lane_B','lane_C']

def mc(rng, probs, n):
    # reproducible multinomial without external dependencies
    cs=[];s=0
    for p in probs:
        s+=p;cs.append(s)
    out=[0]*len(probs)
    for _ in range(n):
        v=rng.random();
        for k,p in enumerate(cs):
            if v<=p:out[k]+=1;break
        else:out[-1]+=1
    return dict(zip(CATS,out))

def make_emission(rng,lane):
    matrix=[]
    for prod in range(4):
        hap=prod//2;splice=prod%2
        h_correct=[.72,.78,.68][lane]+rng.uniform(-.025,.025)
        h_wrong=[.028,.035,.055][lane]+rng.uniform(-.006,.006)
        s_correct=[.80,.74,.83][lane]+rng.uniform(-.02,.02)
        s_wrong=[.07,.10,.055][lane]+rng.uniform(-.01,.01)
        hh=[h_correct,h_wrong,1-h_correct-h_wrong] if hap==0 else [h_wrong,h_correct,1-h_correct-h_wrong]
        ss=[s_correct,s_wrong,1-s_correct-s_wrong] if splice==0 else [s_wrong,s_correct,1-s_correct-s_wrong]
        probs=[hh[0]*ss[0],hh[0]*ss[1],hh[1]*ss[0],hh[1]*ss[1],hh[0]*ss[2],hh[1]*ss[2],hh[2]*ss[0],hh[2]*ss[1],hh[2]*ss[2]]
        matrix.append(probs)
    return matrix

def run(tag):
    seed={'sealed':910241,'hidden':732897,'practice':445012}[tag]
    rng=random.Random(seed)
    if tag=='practice':
        ids=[r['case_id'] for r in json.loads((R/'environment/data/practice.json').read_text())['cases']]
    else:
        original=json.loads((R/'tests'/f'{tag}_cases.json').read_text())['cases']
        ids=[r['case_id'] for r in original]
    arrays=[]; efficiencies=[]
    for lane in range(3):
        arrays.append(make_emission(rng,lane))
        # Opposing allele/isoform capture distortions persist after pooling:
        # uncorrected haplotype-junction ratios are systematically biased.
        # Each batch also has distinct random product efficiencies, recoverable
        # only from its mixed calibration standards.
        capture_bias=([1.40,.68,.68,1.48],
                      [1.52,.61,.61,1.44],
                      [1.34,.62,.57,1.59])
        efficiencies.append([bias*rng.uniform(.94,1.06)
                             for bias in capture_bias[lane]])
    def draw(lane,theta,n):
        w=efficiencies[lane];corr=[a*b for a,b in zip(theta,w)];den=sum(corr)
        q=[v/den for v in corr]
        y=[sum(arrays[lane][j][s]*q[j] for j in range(4)) for s in range(9)]
        return mc(rng,y,n)
    standards=[]
    mixes=[([1 if j==k else 0 for j in range(4)],f'pure_{k}') for k in range(4)]
    # Four distinct three-positive mixed recipes plus two four-positive
    # recipes. The zero entries are biological absences, not invalid reads.
    mixes += [([.25,.25,.25,.25],'balanced'),
              ([.55,0,.40,.05],'a1_rich'),
              ([.05,.50,0,.45],'a2_rich'),
              ([.44,.44,.12,0],'hap0_rich'),
              ([0,.13,.43,.44],'hap1_rich'),
              ([.35,.10,.10,.45],'crossed')]
    for lane in range(3):
        for v,(mix,standard_id) in enumerate(mixes):
            for rep in range(4):
                mix2=mix
                if lane==1 and standard_id=='hap0_rich' and rep==2:
                    mix2=[.02,.94,.02,.02] # contaminated physical aliquot
                standards.append({'batch_id':LANES[lane],'standard_id':standard_id,'replicate_id':f'r{rep+1}',
                                  'known_product_fractions':dict(zip(PRODS,mix)),
                                  'counts':draw(lane,mix2,2300+110*rep)})
    truth=[];measure=[]
    for idx,cid in enumerate(ids):
        h0c=rng.uniform(.23,.78);h1c=rng.uniform(.2,.77)
        mod=rng.uniform(-.045,.045)
        change=(1 if idx%3!=1 else -1)*rng.uniform(.16,.30) if idx%5!=0 else (1 if idx%2 else -1)*rng.uniform(.155,.22)
        h0a=max(.075,min(.925,h0c+mod))
        h1a=max(.075,min(.925,h1c+change))
        # Prospective per-case discriminating strata: no four-target case may
        # be solved by a fixed, central 0.275-wide interval with no analysis.
        # This construction uses latent truth only during offline fixture
        # generation, never inside the verifier or runtime agent image.
        if max(abs(value-.5) for value in (h0a,h0c,h1a,h1c)) < .17:
            h0c=rng.uniform(.19,.26) if idx%2 else rng.uniform(.74,.81)
            h0a=max(.075,min(.925,h0c+mod))
        case_h0=rng.uniform(.36,.63);ctrl_h0=rng.uniform(.36,.63)
        proportions={
            'control':[ctrl_h0*h0c,ctrl_h0*(1-h0c),(1-ctrl_h0)*h1c,(1-ctrl_h0)*(1-h1c)],
            'case':[case_h0*h0a,case_h0*(1-h0a),(1-case_h0)*h1a,(1-case_h0)*(1-h1a)]}
        truth.append({'case_id':cid,'case':{'hap0_a1_fraction':h0a,'hap1_a1_fraction':h1a},
                        'control':{'hap0_a1_fraction':h0c,'hap1_a1_fraction':h1c},
                        'variant_linked_shift':(h1a-h0a)-(h1c-h0c)})
        for group in ('case','control'):
            for lane in range(3):
                n=rng.randint(1400,2100)
                if idx%9==2 and lane==2:n=110  # third batch may be underpowered
                measure.append({'case_id':cid,'group':group,'batch_id':LANES[lane],
                                'replicate_id':'rna1','counts':draw(lane,proportions[group],n)})
    rng.shuffle(standards);rng.shuffle(measure)
    output={'assay_id':'EDA-phased-RNA-v1','synthetic':True,'signal_categories':CATS,'product_classes':PRODS,
            'standards':standards,'measurements':measure}
    if tag!='practice':
        (R/'tests'/f'{tag}_phased_assay.json').write_text(json.dumps(output,sort_keys=True,indent=2)+'\n')
        (R/'tests'/f'{tag}_phased_truth.json').write_text(json.dumps({'records':truth},sort_keys=True,indent=2)+'\n')
    else:
        (R/'environment/data/phased_assay.json').write_text(json.dumps(output,sort_keys=True,indent=2)+'\n')
    print(tag,len(ids),len(standards),len(measure),len(json.dumps(output)))

if __name__=='__main__':
    for name in ['sealed','hidden','practice']:run(name)
