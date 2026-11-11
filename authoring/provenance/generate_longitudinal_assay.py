"""Deterministic synthetic paired-donor EDA splice intervention assay.

The assay reuses independently measured batch-dependent emissions and capture
from the phased RNA standards, but no latent truth is present in the input.
Each donor has four matched aliquots (untreated/treated x baseline/followup).
Donor heterogeneity and deliberately unbalanced batch assignment mean pooled
read ratios can have the wrong sign despite the paired intervention effect.
"""
import importlib.util, json, random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('phased_gen',Path(__file__).with_name('generate_phased_assay.py'))
ph=importlib.util.module_from_spec(spec);spec.loader.exec_module(ph)
CATS=ph.CATS; PRODS=ph.PRODS; LANES=ph.LANES
SEED_OFFSET=0  # Author-only sensitivity parameter; shipped fixture uses zero.

def clamp(x): return max(.02,min(.98,x))

def generate(tag):
    standards_seed={'practice':445012,'sealed':910241,'hidden':732897}[tag]
    rng_cal=random.Random(standards_seed)
    emissions=[];captures=[]
    capture_bias=([1.40,.68,.68,1.48],[1.52,.61,.61,1.44],[1.34,.62,.57,1.59])
    for lane in range(3):
        emissions.append(ph.make_emission(rng_cal,lane))
        captures.append([x*rng_cal.uniform(.94,1.06) for x in capture_bias[lane]])
    rng=random.Random({'practice':714003,'sealed':937331,'hidden':339001}[tag]+SEED_OFFSET)
    if tag=='practice': data=json.loads((ROOT/'environment/data/practice.json').read_text())
    else: data=json.loads((ROOT/'tests'/f'{tag}_cases.json').read_text())
    ids=sorted(row['case_id'] for row in data['cases'])
    def draw(lane,mix,n):
        scaled=[a*b for a,b in zip(mix,captures[lane])]; total=sum(scaled)
        q=[x/total for x in scaled]
        observed=[sum(q[j]*emissions[lane][j][k] for j in range(4)) for k in range(9)]
        return ph.mc(rng,observed,n)
    measure=[]; truth=[]
    for idx,cid in enumerate(ids):
        donor_truth=[]
        base_effect= rng.choice((-1,1))*rng.uniform(.39,.455)
        for di in range(6):
            donor=f'D{di+1:02d}'
            h0= rng.uniform(.25,.70);h1=rng.uniform(.39,.47) if base_effect>0 else rng.uniform(.53,.61)
            trend=rng.uniform(-.065,.065)
            base_condition=rng.uniform(-.055,.055)
            background=rng.uniform(-.025,.025)
            donor_effect=base_effect+rng.uniform(-.019,.019)
            prop={};contrast={}
            for treatment in ('untreated','treated'):
                for time in ('baseline','followup'):
                    follow=time=='followup'; treated=treatment=='treated'
                    a0=clamp(h0+int(follow)*trend+int(treated)*base_condition+int(follow and treated)*background)
                    a1=clamp(h1+int(follow)*trend*.6+int(treated)*base_condition*.65+int(follow and treated)*(background+donor_effect))
                    hap0=rng.uniform(.35,.64)
                    prop[(treatment,time)]=[hap0*a0,hap0*(1-a0),(1-hap0)*a1,(1-hap0)*(1-a1)]
                    contrast[(treatment,time)]=a1-a0
            effect=(contrast['treated','followup']-contrast['treated','baseline'])-(contrast['untreated','followup']-contrast['untreated','baseline'])
            donor_truth.append({'donor_id':donor,'treatment_time_effect':effect})
            for t in ('untreated','treated'):
                for tm in ('baseline','followup'):
                    # Systematic treatment/time lane imbalance requires explicit
                    # chemistry correction, not pooling or donor-only differences.
                    base={('untreated','baseline'):1,('untreated','followup'):2,('treated','baseline'):0,('treated','followup'):1}[(t,tm)]
                    lane=(base+(1 if di%4==0 else 0)+(1 if di%5==0 else 0))%3
                    n=rng.randint(1700,2500)
                    measure.append({'case_id':cid,'donor_id':donor,'condition':t,'time':tm,'batch_id':LANES[lane],
                                   'replicate_id':'library1','counts':draw(lane,prop[(t,tm)],n)})
        mean=sum(x['treatment_time_effect'] for x in donor_truth)/len(donor_truth)
        truth.append({'case_id':cid,'mean_treatment_time_effect':mean,'donors':donor_truth})
    rng.shuffle(measure)
    assay={'assay_id':'EDA-paired-perturbation-v1','synthetic':True,
           'calibration_source':'phased_assay.standards','measurements':measure}
    if tag=='practice':
        out=ROOT/'environment/data/longitudinal_assay.json'
    else:
        out=ROOT/'tests'/f'{tag}_longitudinal_assay.json'
        (ROOT/'tests'/f'{tag}_longitudinal_truth.json').write_text(json.dumps({'records':truth},indent=2,sort_keys=True)+'\n')
    out.write_text(json.dumps(assay,indent=2,sort_keys=True)+'\n')
    print(tag,len(ids),len(measure),out.stat().st_size)
if __name__=='__main__':
    for label in ('practice','sealed','hidden'):generate(label)
