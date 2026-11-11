#!/usr/bin/env python3
"""Entrypoint for the EDA isoform reconstruction and junction-evidence workflow."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coding import compute
from evidence import summarize
from audit import check_reference
from inference import infer, reconcile_batches, reconcile_phased_rna, reconcile_longitudinal, reconcile_phase_graph

def run():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference',required=True)
    ap.add_argument('--variants',required=True)
    ap.add_argument('--observations',required=True)
    ap.add_argument('--fasta',required=True)
    ap.add_argument('--output-dir',required=True)
    ap.add_argument('--calibration',required=True)
    ap.add_argument('--curated',required=True)
    ap.add_argument('--batch-assay',required=True)
    ap.add_argument('--phased-assay',required=True)
    ap.add_argument('--longitudinal-assay',required=True)
    ap.add_argument('--phase-graph',required=True)
    args=ap.parse_args()
    ref=json.loads(Path(args.reference).read_text())
    cases=json.loads(Path(args.variants).read_text())
    obs=json.loads(Path(args.observations).read_text())
    calibration=json.loads(Path(args.calibration).read_text())
    assay=json.loads(Path(args.batch_assay).read_text())
    phased=json.loads(Path(args.phased_assay).read_text())
    longitudinal=json.loads(Path(args.longitudinal_assay).read_text())
    graph=json.loads(Path(args.phase_graph).read_text())
    dna_ids=[c['case_id'] for c in cases['cases']]
    if any(type(x) is not str or not x for x in dna_ids) or len(dna_ids)!=len(set(dna_ids)):
        raise ValueError('Duplicate or invalid DNA case identifiers')
    result=compute(ref,cases)
    evidence=summarize(obs)
    audit=check_reference(ref,args.fasta,args.curated)
    seen={row['case_id'] for row in result['records']}
    if set(row['case_id'] for row in evidence)!=seen:
        raise ValueError('RNA evidence and DNA variant case IDs do not match')
    if {x['case_id'] for x in assay['measurements']}!=seen:
        raise ValueError('Batch assay and DNA variant case IDs do not match')
    if {x['case_id'] for x in phased['measurements']}!=seen:
        raise ValueError('Phased allele/RNA data and DNA variant case IDs do not match')
    if {x['case_id'] for x in longitudinal['measurements']}!=seen:
        raise ValueError('Paired donor assay case IDs do not match DNA cases')
    if {x['case_id'] for x in graph['cases']}!=seen:
        raise ValueError('Phase graph case IDs do not match DNA cases')
    out=Path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    for name,content in [('coding.json',result),('evidence.json',{'records':evidence}),('audit.json',audit),('likelihood.json',infer(evidence,calibration)),('batch_inference.json',reconcile_batches(assay)),('phased_rna.json',reconcile_phased_rna(phased)),('longitudinal.json',reconcile_longitudinal(longitudinal,phased)),('phase_graph.json',reconcile_phase_graph(graph))]:
        (out/name).write_text(json.dumps(content,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':run()
