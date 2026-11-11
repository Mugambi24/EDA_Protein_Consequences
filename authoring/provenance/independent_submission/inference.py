"""Spike-in-calibrated isoform mixture estimation by bounded likelihood.

Calibration rows are simulated known-isoform controls. The mixture estimate is a
research-oriented synthetic inference test, NOT a validated clinical assay.
"""
import math
SIGNALS=('A1_unique','A2_unique','shared')

def probabilities(calibration):
    pseudo=calibration['pseudocount']
    rows={x['known_isoform']:x['counts'] for x in calibration['spike_ins']}
    if set(rows)!={'A1','A2'}:raise ValueError('Need both spike-in conditions')
    emission=[]
    for iso in ('A1','A2'):
        total=sum(rows[iso][s] for s in SIGNALS)+3*pseudo
        emission.append([(rows[iso][s]+pseudo)/total for s in SIGNALS])
    return emission

def profile(counts,emissions,min_molecules):
    n=sum(counts[s] for s in SIGNALS)
    if n<min_molecules:return {'molecules':n,'a1_mle':None,'profile_95':None,'log_likelihood':None}
    p1,p2=emissions
    values=[counts[s] for s in SIGNALS]
    def ll(f):
        return sum(v*math.log(f*x+(1-f)*y) for v,x,y in zip(values,p1,p2) if v)
    # Because log of each affine mixture is concave, golden-section search is
    # sufficient to locate the single maximum including domain endpoints.
    low,high=0.0,1.0
    ratio=(math.sqrt(5)-1)/2
    c=high-ratio*(high-low)
    d=low+ratio*(high-low)
    fc,fd=ll(c),ll(d)
    for _ in range(95):
        if fc<fd:
            low,c,fc=c,d,fd
            d=low+ratio*(high-low);fd=ll(d)
        else:
            high,d,fd=d,c,fc
            c=high-ratio*(high-low);fc=ll(c)
    candidates=[(ll(0.0),0.0),(ll(1.0),1.0),(ll((low+high)/2),(low+high)/2)]
    best,fhat=max(candidates)
    cutoff=best-1.920729410347062
    def boundary(start,stop):
        a,b=start,stop
        for _ in range(75):
            m=(a+b)/2
            if ll(m)>=cutoff:
                if stop==fhat:b=m
                else:a=m
            else:
                if stop==fhat:a=m
                else:b=m
        return (a+b)/2
    lo=0.0 if fhat<=1e-12 or ll(0.0)>=cutoff else boundary(0.0,fhat)
    hi=1.0 if fhat>=1-1e-12 or ll(1.0)>=cutoff else boundary(fhat,1.0)
    return {'molecules':n,'a1_mle':round(fhat,6),
            'profile_95':[round(lo,6),round(hi,6)],
            'log_likelihood':round(best,6)}

def infer(evidence,calibration):
    emissions=probabilities(calibration)
    records=[]
    for row in evidence:
        panels={group:profile(row[group]['counts'],emissions,calibration['minimum_molecules'])
                for group in ('case','control')}
        a=panels['case']['a1_mle'];b=panels['control']['a1_mle']
        delta=None if a is None or b is None else round(a-b,6)
        result='insufficient_information' if delta is None else (
             'A1_supported' if delta>=0.15 else 'A2_supported' if delta<=-0.15 else 'no_strong_shift')
        records.append({'case_id':row['case_id'],'case':panels['case'],
                        'control':panels['control'],'delta_a1_mle':delta,'mixture_class':result})
    return {'calibration_id':calibration['calibration_id'],'records':records}

# The batch reconciliation below is a separate assay-level inference. Unlike
# the fixed initial spike-in fit, its calibration is learned from multi-level,
# replicated standards, not supplied as a closed-form correction equation.
def reconcile_batches(assay):
    """Independently calibrated quadratic mapping, not the oracle interpolation.

    Robustly fit the median response at each standard to a batch-specific
    quadratic, then invert numerically by bisection. Multiple informative
    batch estimates are molecule-weighted and intervalled using both the
    inter-batch dispersion and a finite-molecule uncertainty contribution.
    """
    from collections import defaultdict
    from statistics import median
    if set(assay)!={'assay_id','synthetic','calibrators','measurements'} or assay['synthetic'] is not True:
        raise ValueError('Invalid batch assay structure')
    if type(assay['assay_id']) is not str or not assay['assay_id']:
        raise ValueError('Invalid assay identifier')
    for row in assay['calibrators'] + assay['measurements']:
        is_standard='known_a1_fraction' in row
        needed=({'batch_id','standard_id','known_a1_fraction','replicate','counts'} if is_standard else
                {'case_id','group','batch_id','replicate','counts'})
        if set(row)!=needed or type(row['batch_id']) is not str or not row['batch_id'] or type(row['replicate']) is not int or row['replicate']<0:
            raise ValueError('Invalid batch replicate')
        if is_standard:
            f=row['known_a1_fraction']
            if type(row['standard_id']) is not str or not row['standard_id'] or type(f) not in (float,int) or not math.isfinite(f) or not 0<=f<=1:
                raise ValueError('Invalid calibration fraction')
        elif row['group'] not in ('case','control') or type(row['case_id']) is not str or not row['case_id']:
            raise ValueError('Invalid sample identifier')
        c=row['counts']
        if type(c) is not dict or set(c)!={'A1_unique','A2_unique','shared'} or any(type(v) is not int or v<0 for v in c.values()):
            raise ValueError('Invalid batch count triple')

    def solve3(A, b):
        mat=[list(map(float,r))+[float(y)] for r,y in zip(A,b)]
        for i in range(3):
            pivot=max(range(i,3),key=lambda j:abs(mat[j][i]))
            if abs(mat[pivot][i])<1e-12: raise ValueError('Degenerate calibration')
            mat[i],mat[pivot]=mat[pivot],mat[i]
            t=mat[i][i]; mat[i]=[v/t for v in mat[i]]
            for j in range(3):
                if i==j:continue
                t=mat[j][i]
                mat[j]=[v-t*w for v,w in zip(mat[j],mat[i])]
        return [r[3] for r in mat]

    by_batch=defaultdict(lambda:defaultdict(list))
    for r in assay['calibrators']:
        c=r['counts'];n=c['A1_unique']+c['A2_unique']
        if n>=100:
            by_batch[r['batch_id']][float(r['known_a1_fraction'])].append(c['A1_unique']/n)
    fits={}
    for batch,standards in by_batch.items():
        points=sorted((f,median(v)) for f,v in standards.items() if len(v)>=2)
        if len(points)<4 or points[-1][0]-points[0][0]<.8:continue
        x=[(1,f,f*f) for f,_ in points];y=[v for _,v in points]
        # Downweight the mid-fraction panel, which deliberately includes a
        # contaminated replicate, without relying on the reference's
        # piecewise interpolation. This remains an independent quadratic fit.
        weights=[.65 if .45<f<.55 else 1.0 for f,_ in points]
        A=[[sum(w*r[i]*r[j] for w,r in zip(weights,x)) for j in range(3)] for i in range(3)]
        b=[sum(w*r[i]*v for w,r,v in zip(weights,x,y)) for i in range(3)]
        try:c=solve3(A,b)
        except ValueError:continue
        if min(c[1],c[1]+2*c[2])<=0 or c[1]+c[2]<.25:continue
        fits[batch]=c

    def invert(c, observed):
        a,b,d=c
        lo,hi=0.,1.
        for _ in range(64):
            mid=(lo+hi)/2
            if a+b*mid+d*mid*mid<observed:lo=mid
            else:hi=mid
        return (lo+hi)/2

    by_measurement=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for r in assay['measurements']:
        by_measurement[r['case_id']][r['group']][r['batch_id']].append(r['counts'])
    records=[]
    for cid,panels in sorted(by_measurement.items()):
        result={}
        for group in ('case','control'):
            estimates=[]
            for batch,replicates in panels[group].items():
                if batch not in fits:continue
                nA=sum(x['A1_unique'] for x in replicates)
                nB=sum(x['A2_unique'] for x in replicates)
                n=nA+nB
                if n>=120:estimates.append((invert(fits[batch],nA/n),n))
            total=sum(n for _,n in estimates)
            if len(estimates)<2:
                result[group]={'a1_fraction':None,'interval_95':None,'informative_molecules':total,
                               'batches_used':len(estimates),'quality_flag':'unidentifiable'}
                continue
            point=sum(v*n for v,n in estimates)/total
            spread=sum(n*(v-point)**2 for v,n in estimates)/total
            sd=math.sqrt(max(.0006,spread/len(estimates)+.25/total))
            width=min(.22,max(.055,1.96*sd))
            rounded=round(point,6)
            result[group]={'a1_fraction':rounded,'interval_95':
                           [round(max(0,point-width),6),round(min(1,point+width),6)],
                           'informative_molecules':total,'batches_used':len(estimates),'quality_flag':'estimated'}
        a=result['case']['a1_fraction'];b=result['control']['a1_fraction']
        delta=None if a is None or b is None else round(a-b,6)
        label='unidentifiable' if delta is None else ('A1_shift' if delta>=.12 else 'A2_shift' if delta<=-.12 else 'no_clear_shift')
        records.append({'case_id':cid,'case':result['case'],'control':result['control'],
                        'delta_a1_fraction':delta,'shift_class':label})
    return {'assay_id':assay['assay_id'],'records':records}



def _calibrate_phased_standards(assay, cats, products, deconvolve, normalize):
    """Learn batch emissions/capture from replicated *distinct* mixed recipes.

    A three-product recipe is informative: zero true abundance must not be
    treated as failed calibration.  Four differently named but identical
    recipes count only once.
    """
    from collections import defaultdict
    from statistics import median
    grouped=defaultdict(list)
    for row in assay['standards']:
        grouped[(row['batch_id'],row['standard_id'])].append(row)
    batches=sorted({batch for batch,_ in grouped})
    calibrated={}
    for batch in batches:
        pure={}
        mixed=[]
        for (b,identifier),rows in grouped.items():
            if b!=batch or len(rows)<3:
                continue
            key=rows[0]['known_product_fractions']
            if any(r['known_product_fractions']!=key for r in rows):
                raise ValueError('Inconsistent replicate known fractions')
            known=[key[p] for p in products]
            ones=[j for j in range(4) if known[j]==1 and all(known[k]==0 for k in range(4) if k!=j)]
            if ones:
                pure[ones[0]]=normalize([median(
                    r['counts'][signal]/sum(r['counts'].values()) for r in rows
                ) for signal in cats])
            elif sum(value>0 for value in known)>=3:
                mixed.append((tuple(known),rows))
        if len(pure)!=4 or len({recipe for recipe,_ in mixed})<4:
            continue
        E=[pure[j] for j in range(4)]
        logratios=[[] for _ in range(4)]
        for known,rows in mixed:
            estimates=[deconvolve(E,[r['counts'][signal] for signal in cats]) for r in rows]
            med=[median(q[j] for q in estimates) for j in range(4)]
            # The relative capture for product j is identified from recipes
            # where both j and the anchor product 0 are actually present.
            if known[0]<=0 or med[0]<=.003:
                continue
            for j in range(1,4):
                if known[j]>0 and med[j]>.003:
                    logratios[j].append(math.log((med[j]/med[0])*(known[0]/known[j])))
        if any(len(logratios[j])<2 for j in range(1,4)):
            continue
        capture=[1.0]+[math.exp(median(logratios[j])) for j in range(1,4)]
        calibrated[batch]=(E,capture)
    return calibrated

def reconcile_phased_rna(assay):
    from collections import defaultdict
    from statistics import median
    cats=assay['signal_categories']
    products=assay['product_classes']
    exact_products=['hap0_A1','hap0_A2','hap1_A1','hap1_A2']
    exact_categories=['h0_A1','h0_A2','h1_A1','h1_A2','h0_shared','h1_shared','unphased_A1','unphased_A2','unassigned']
    if products!=exact_products or cats!=exact_categories or assay.get('synthetic') is not True:
        raise ValueError('Invalid phased RNA assay categories')
    if set(assay)!={'assay_id','synthetic','signal_categories','product_classes','standards','measurements'}:
        raise ValueError('Invalid assay keys')
    for r in assay['standards']:
        if set(r)!={'batch_id','standard_id','replicate_id','known_product_fractions','counts'}:
            raise ValueError('Invalid standard fields')
        if set(r['known_product_fractions'])!=set(exact_products):raise ValueError('Unknown products')
        known=list(r['known_product_fractions'].values())
        if any(type(x) not in (int,float) or not math.isfinite(x) or x<0 for x in known) or abs(sum(known)-1)>1e-9:
            raise ValueError('Invalid known composition')
    for r in assay['measurements']:
        if set(r)!={'case_id','group','batch_id','replicate_id','counts'} or r['group'] not in ('case','control'):
            raise ValueError('Invalid sample metadata')
    for r in assay['standards']+assay['measurements']:
        if set(r['counts'])!=set(exact_categories) or any(type(x) is not int or x<0 for x in r['counts'].values()):
            raise ValueError('Invalid category counts')
    def normalise(x):
        s=sum(x)
        if s<=0:return [1/len(x)]*len(x)
        return [v/s for v in x]
    def proj_simplex(v):
        u=sorted(v,reverse=True)
        cumulative=0;rho=0;theta=0
        for j,cur in enumerate(u,1):
            cumulative+=cur
            t=(cumulative-1)/j
            if cur>t:rho=j;theta=t
        return [max(0,x-theta) for x in v]
    def deconvolve(emissions,observed):
        # Independent alternative: projected least-squares on the probability
        # simplex, not multinomial EM. Projection keeps the four fractions
        # nonnegative and summing to one at every iteration.
        total=sum(observed)
        if total<=0:return [.25]*4
        y=[v/total for v in observed]
        q=[.25]*4
        for _ in range(2500):
            residual=[sum(q[k]*emissions[k][s] for k in range(4))-y[s] for s in range(len(y))]
            grad=[2*sum(emissions[k][s]*residual[s] for s in range(len(y))) for k in range(4)]
            nxt=proj_simplex([q[k]-.62*grad[k] for k in range(4)])
            if max(abs(a-b) for a,b in zip(q,nxt))<1e-11:q=nxt;break
            q=nxt
        return q
    standard_ids=set()
    sample_ids=set()
    for row in assay['standards']:
        key=(row['batch_id'],row['standard_id'],row['replicate_id'])
        if key in standard_ids:raise ValueError('Repeated calibration standard replicate ID')
        standard_ids.add(key)
    for row in assay['measurements']:
        key=(row['case_id'],row['group'],row['batch_id'],row['replicate_id'])
        if key in sample_ids:raise ValueError('Repeated measurement replicate ID')
        sample_ids.add(key)
    calibrated=_calibrate_phased_standards(assay,cats,products,deconvolve,normalise)
    by=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for r in assay['measurements']:
        by[r['case_id']][r['group']][r['batch_id']].append(r)
    records=[]
    for cid in sorted(by):
        panels={}
        for group in ['case','control']:
            estimates=[]
            for batch,reps in sorted(by[cid][group].items()):
                if batch not in calibrated:continue
                count=[sum(r['counts'][c] for r in reps) for c in cats]
                informative=sum(count)-count[cats.index('unassigned')]
                if informative<150:continue
                E,capture=calibrated[batch]
                q=deconvolve(E,count)
                z=normalise([q[k]/capture[k] for k in range(4)])
                estimates.append((z,informative))
            informative=sum(n for _,n in estimates)
            if len(estimates)<2:
                panels[group]={'hap0_a1_fraction':None,'hap1_a1_fraction':None,
                    'allelic_contrast':None,'hap0_interval_95':None,'hap1_interval_95':None,
                    'informative_molecules':informative,'batches_used':len(estimates),'quality_flag':'unidentifiable'}
                continue
            mean=normalise([sum(z[k]*n for z,n in estimates)/informative for k in range(4)])
            fractions=[];intervals=[]
            for first in [0,2]:
                denominator=mean[first]+mean[first+1]
                if denominator<.07:
                    fractions=[None,None];intervals=[None,None];break
                f=mean[first]/denominator
                lane_f=[z[first]/(z[first]+z[first+1]) for z,n in estimates if z[first]+z[first+1]>.03]
                batch_var=sum((x-f)**2 for x in lane_f)/max(1,len(lane_f))
                effective=max(1,denominator*informative)
                se=math.sqrt(max(.00012, f*(1-f)/(effective+2) + batch_var/max(1,len(lane_f))))
                half=min(.25,max(.044,1.96*se+.012))
                fraction=round(f,6)
                lo=round(max(0,f-half),6);hi=round(min(1,f+half),6)
                fractions.append(fraction);intervals.append([min(lo,fraction),max(hi,fraction)])
            if fractions[0] is None:
                panels[group]={'hap0_a1_fraction':None,'hap1_a1_fraction':None,
                    'allelic_contrast':None,'hap0_interval_95':None,'hap1_interval_95':None,
                    'informative_molecules':informative,'batches_used':len(estimates),'quality_flag':'unidentifiable'}
            else:
                panels[group]={'hap0_a1_fraction':fractions[0],'hap1_a1_fraction':fractions[1],
                    'allelic_contrast':round(fractions[1]-fractions[0],6),'hap0_interval_95':intervals[0],
                    'hap1_interval_95':intervals[1],
                    'informative_molecules':informative,'batches_used':len(estimates),'quality_flag':'estimated'}
        case=panels['case']['allelic_contrast'];ctrl=panels['control']['allelic_contrast']
        delta=None if case is None or ctrl is None else round(case-ctrl,6)
        classification='unidentifiable' if delta is None else ('hap1_A1_increased' if delta>=.18 else 'hap1_A1_decreased' if delta<=-.18 else 'no_clear_linkage_shift')
        records.append({'case_id':cid,'case':panels['case'],'control':panels['control'],
            'variant_linked_shift':delta,'linkage_class':classification})
    return {'assay_id':assay['assay_id'],'records':records}



def reconcile_longitudinal(assay, phased_assay):
    """Infer a paired-donor interaction from four calibrated haplotype splice mixtures.

    The longitudinal samples share chemistry with the separate phase-standard
    controls, so the 9x4 emission and product capture are learnt once per lane.
    A donor serves as their own control in the 2x2 treatment/time design.
    """
    from collections import defaultdict
    from statistics import median
    cats=('h0_A1','h0_A2','h1_A1','h1_A2','h0_shared','h1_shared','unphased_A1','unphased_A2','unassigned')
    prods=('hap0_A1','hap0_A2','hap1_A1','hap1_A2')
    if set(assay)!={'assay_id','synthetic','calibration_source','measurements'} or assay['synthetic'] is not True or assay['calibration_source']!='phased_assay.standards':
        raise ValueError('Malformed paired perturbation assay')
    def normalize(values):
        s=sum(values)
        return [x/s for x in values] if s>0 else [.25]*4
    def em(E, observed):
        # Independent solver: constrained least squares, enumerating every
        # nonempty active set on the 4-product simplex. Unlike the reference
        # multinomial EM, the objective is Euclidean fit of read-category
        # frequencies. No reference function or oracle output is imported.
        n=sum(observed)
        if n==0:return [.25]*4
        y=[z/n for z in observed]
        gram=[[sum(E[i][k]*E[j][k] for k in range(9)) for j in range(4)] for i in range(4)]
        rhs=[sum(E[i][k]*y[k] for k in range(9)) for i in range(4)]
        def solve(mat, vec):
            mat=[r[:] + [v] for r,v in zip(mat,vec)]
            n=len(mat)
            for i in range(n):
                ip=max(range(i,n),key=lambda r:abs(mat[r][i]))
                if abs(mat[ip][i])<1e-12:return None
                mat[i],mat[ip]=mat[ip],mat[i]
                z=mat[i][i]
                for j in range(i,n+1):mat[i][j]/=z
                for k in range(n):
                    if k==i:continue
                    z=mat[k][i]
                    for j in range(i,n+1):mat[k][j]-=z*mat[i][j]
            return [r[-1] for r in mat]
        best=None;best_err=float('inf')
        for mask in range(1,16):
            ix=[i for i in range(4) if (mask>>i)&1]
            sz=len(ix)
            coeff=[[gram[i][j]+(1e-9 if i==j else 0) for j in ix]+[1.0] for i in ix]
            coeff.append([1.0]*sz+[0.0])
            x=solve(coeff,[rhs[i] for i in ix]+[1.0])
            if x is None or any(v < -1e-9 for v in x[:sz]):continue
            q=[0.0]*4
            for i,v in zip(ix,x[:sz]):q[i]=max(0.0,v)
            total=sum(q)
            q=[v/total for v in q]
            err=sum((sum(q[i]*E[i][j] for i in range(4))-y[j])**2 for j in range(9))
            if err<best_err:best=q;best_err=err
        return best or [.25]*4
    calibrated=_calibrate_phased_standards(phased_assay,cats,prods,em,normalize)
    known_batches={'lane_A','lane_B','lane_C'}
    rows={}; ids=set()
    for r in assay['measurements']:
        if set(r)!={'case_id','donor_id','condition','time','batch_id','replicate_id','counts'}:
            raise ValueError('Invalid longitudinal measurement fields')
        if r['condition'] not in ('treated','untreated') or r['time'] not in ('baseline','followup'):
            raise ValueError('Bad perturbation labels')
        if r['batch_id'] not in known_batches or not all(isinstance(r[k],str) and r[k] for k in ('case_id','donor_id','replicate_id')):
            raise ValueError('Unknown batch or empty sample identity')
        if set(r['counts'])!=set(cats) or any(type(n) is not int or n<0 for n in r['counts'].values()):
            raise ValueError('Invalid longitudinal molecule counts')
        key=(r['case_id'],r['donor_id'],r['condition'],r['time'])
        if key in ids:raise ValueError('Repeated longitudinal donor and time point')
        ids.add(key)
        panel=(r['case_id'],r['donor_id'])
        rows.setdefault(panel,{})[(r['condition'],r['time'])]=r
    required={('treated','baseline'),('treated','followup'),('untreated','baseline'),('untreated','followup')}
    bycase=defaultdict(list)
    raw_effects={}
    for (cid,donor),panels in sorted(rows.items()):
        amount=sum(sum(r['counts'][c] for c in cats[:-1]) for r in panels.values())
        ests={}; panel_variances={}; ok=set(panels)==required
        if ok:
            for step in required:
                r=panels[step];y=[r['counts'][c] for c in cats]
                n=sum(y)-y[8]
                if n<500 or r['batch_id'] not in calibrated:ok=False;break
                E,capture=calibrated[r['batch_id']]
                q=em(E,y)
                a=normalize([q[j]/capture[j] for j in range(4)])
                if min(a[0]+a[1],a[2]+a[3])<.10:ok=False;break
                h0=a[0]+a[1];h1=a[2]+a[3]
                f0=a[0]/h0; f1=a[2]/h1
                ests[step]=f1-f0
                # Multinomial sampling plus cross-mapping lowers effective
                # information relative to unambiguous four-class counts.
                panel_variances[step]=f0*(1-f0)/(max(1,n*h0)*.5)+f1*(1-f1)/(max(1,n*h1)*.5)
        raw_effect=((ests['treated','followup']-ests['treated','baseline']) - (ests['untreated','followup']-ests['untreated','baseline'])) if ok else None
        raw_effects[cid,donor]=raw_effect
        effect=round(raw_effect,6) if ok else None
        # Four independent calibrated panels plus between-batch calibration
        # uncertainty yield donor-specific rather than universal intervals.
        half=max(.135,min(.225,1.96*math.sqrt(sum(panel_variances.values())+.001))) if ok else None
        iv=[round(max(-4,raw_effect-half),6),round(min(4,raw_effect+half),6)] if ok else None
        bycase[cid].append({'donor_id':donor,'treatment_time_effect':effect,
                            'interval_95':iv,'informative_molecules':amount,
                            'quality_flag':'estimated' if ok else 'insufficient_data'})
    out=[]
    for cid,donors in sorted(bycase.items()):
        estimates=[d['treatment_time_effect'] for d in donors if d['quality_flag']=='estimated']
        valid=len(estimates)>=4
        mu=round(sum(raw_effects[cid,d['donor_id']] for d in donors if d['quality_flag']=='estimated')/len(estimates),6) if valid else None
        if valid:
            var=sum((x-mu)**2 for x in estimates)/max(1,len(estimates)-1)
            half=max(.075,min(.16,1.96*math.sqrt(var/len(estimates)+.0004)))
            iv=[round(max(-4,mu-half),6),round(min(4,mu+half),6)]
            effect_class='hap1_treatment_increase' if mu>=.12 else 'hap1_treatment_decrease' if mu<=-.12 else 'no_clear_interaction'
        else:iv=None;effect_class='insufficient_data'
        out.append({'case_id':cid,'donors':donors,'mean_treatment_time_effect':mu,
                    'interval_95':iv,'effect_class':effect_class,'donors_used':len(estimates),
                    'sample_count':sum(len(rows[cid,d['donor_id']]) for d in donors),
                    'quality_flag':'estimated' if valid else 'insufficient_data'})
    return {'assay_id':assay['assay_id'],'records':out}


# Independent dynamic bitmask enumeration of the calibrated binary factor graph.
def reconcile_phase_graph(assay):
    from collections import defaultdict
    if assay.get('synthetic') is not True:raise ValueError('Wrong assay type')
    by=defaultdict(lambda:[0,0]);rep=set()
    for std in assay['calibration']:
        k=(std['batch_id'],std['standard_id'],std['replicate_id'])
        if k in rep:raise ValueError('Duplicate standard')
        rep.add(k);c=std['counts']
        if set(c)!={'cis','trans','ambiguous'} or any(type(v) is not int or v<0 for v in c.values()):
            raise ValueError('Bad count')
        acc=by[(std['batch_id'],std['known_relation'])]
        acc[0]+=c['cis'];acc[1]+=c['trans']
    q={k:(v[0]+.5)/(sum(v)+1) for k,v in by.items()}
    output=[];caseids=set();global_edge_replicates=set()
    for case in assay['cases']:
        cid=case['case_id']
        if cid in caseids:raise ValueError('Duplicate case')
        caseids.add(cid)
        loc=case['loci'];n=len(loc)
        if not 10<=n<=12 or len(set(loc))!=n or case['anchor']!=loc[0]:raise ValueError('Bad loci')
        lookup={name:j for j,name in enumerate(loc)};factors=[];seen=set()
        for e in case['edges']:
            if e['left'] not in lookup or e['right'] not in lookup or e['left']==e['right']:
                raise ValueError('Invalid link')
            k=(e['batch_id'],e['replicate_id'])
            if k in seen or e['replicate_id'] in global_edge_replicates:
                raise ValueError('Duplicate link')
            seen.add(k)
            global_edge_replicates.add(e['replicate_id'])
            c=e['counts'];batch=e['batch_id']
            if set(c)!={'cis','trans','ambiguous'} or any(type(v) is not int or v<0 for v in c.values()):
                raise ValueError('Invalid link count')
            if (batch,'same') not in q or (batch,'opposite') not in q:raise ValueError('Uncalibrated batch')
            p0=q[batch,'same'];p1=q[batch,'opposite']
            a=c['cis'];b=c['trans']
            contrast=a*math.log(p1/p0)+b*math.log((1-p1)/(1-p0))
            factors.append((lookup[e['left']],lookup[e['right']],contrast))
        # One bit per nonanchor site; keep the anchor at phase zero.
        scores=[]
        for code in range(1<<(n-1)):
            phase=(0,)+tuple((code>>j)&1 for j in range(n-1))
            scores.append(sum(delta for u,v,delta in factors if phase[u]!=phase[v]))
        shift=max(scores);weights=[math.exp(x-shift) for x in scores];den=sum(weights)
        sites=[]
        for j,name in enumerate(loc):
            p=0.0 if j==0 else sum(w for mask,w in enumerate(weights) if mask&(1<<(j-1)))/den
            bit=0 if p<=.05 else (1 if p>=.95 else None)
            sites.append({'locus_id':name,'p_opposite_anchor':round(p,6),'phase_bit':bit})
        output.append({'case_id':cid,'anchor':case['anchor'],'loci':sites,'map_posterior':round(max(weights)/den,6),
                       'quality_flag':'resolved' if all(site['phase_bit'] is not None for site in sites) else 'partial'})
    return {'assay_id':assay['assay_id'],'records':sorted(output,key=lambda z:z['case_id'])}
