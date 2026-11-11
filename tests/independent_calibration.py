"""Independent spike-in likelihood audit via likelihood-derivative bisection.

This does not import the reference implementation. Values are used for numerical
cross-validation and not copied into the agent's visible environment.
"""
import json, math
CATS=('A1_unique','A2_unique','shared')

def independent_panel(counts,cal):
    n=sum(counts[k] for k in CATS)
    if n<cal['minimum_molecules']:
        return {'molecules':n,'a1_mle':None,'profile_95':None,'log_likelihood':None}
    spike={x['known_isoform']:x['counts'] for x in cal['spike_ins']}
    alpha=cal['pseudocount']
    probs=[]
    for iso in ('A1','A2'):
        total=sum(spike[iso].values())+3*alpha
        probs.append([(spike[iso][c]+alpha)/total for c in CATS])
    x,y=probs
    o=[counts[k] for k in CATS]
    def f(t): return sum(v*math.log(t*a+(1-t)*b) for v,a,b in zip(o,x,y) if v)
    def gradient(t): return sum(v*(a-b)/(t*a+(1-t)*b) for v,a,b in zip(o,x,y) if v)
    if gradient(0)<=0: estimate=0.
    elif gradient(1)>=0: estimate=1.
    else:
        l,r=0.,1.
        for _ in range(90):
            m=(l+r)/2
            if gradient(m)>0:l=m
            else:r=m
        estimate=(l+r)/2
    cutoff=f(estimate)-1.920729410347062
    def endpoint(l,r,lower):
        for _ in range(90):
            m=(l+r)/2
            if (f(m)>=cutoff)==lower: r=m
            else: l=m
        return (l+r)/2
    low=0. if estimate<1e-12 or f(0.)>=cutoff else endpoint(0,estimate,True)
    high=1. if estimate>1-1e-12 or f(1.)>=cutoff else endpoint(estimate,1.,False)
    return {'molecules':n,'a1_mle':round(estimate,6),'profile_95':[round(low,6),round(high,6)],'log_likelihood':round(f(estimate),6)}

def calculate(rows,cal):
    output=[]
    for r in rows:
        a=independent_panel(r['case']['counts'],cal)
        b=independent_panel(r['control']['counts'],cal)
        delta=None if a['a1_mle'] is None or b['a1_mle'] is None else round(a['a1_mle']-b['a1_mle'],6)
        cat=('insufficient_information' if delta is None else 'A1_supported' if delta>=0.15 else 'A2_supported' if delta<=-0.15 else 'no_strong_shift')
        output.append({'case_id':r['case_id'],'case':a,'control':b,'delta_a1_mle':delta,'mixture_class':cat})
    return {'calibration_id':cal['calibration_id'],'records':output}

