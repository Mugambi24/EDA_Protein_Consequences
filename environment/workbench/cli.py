#!/usr/bin/env python3
"""EDA multi-stage computational workbench, currently incomplete.

The task is to repair and integrate the modules while extending the CLI to the contract in instruction.md.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coding import compute
from evidence import summarize
from audit import check_reference
from inference import infer

def main():
    parser=argparse.ArgumentParser()
    for flag in ('reference','variants','observations','fasta','output-dir'):
        parser.add_argument('--'+flag,required=True)
    parser.add_argument('--calibration',required=True)
    parser.add_argument('--curated',required=True)
    parser.add_argument('--batch-assay',required=True)
    parser.add_argument('--phased-assay',required=True)
    parser.add_argument('--longitudinal-assay',required=True)
    parser.add_argument('--phase-graph',required=True)
    args=vars(parser.parse_args())
    reference=json.loads(Path(args['reference']).read_text())
    cases=json.loads(Path(args['variants']).read_text())
    observations=json.loads(Path(args['observations']).read_text())
    calibration=json.loads(Path(args['calibration']).read_text())
    summary=summarize(observations)
    outputs={
       'coding.json':compute(reference,cases),
       'evidence.json':{'records':summary},
       'likelihood.json':infer(summary,calibration),
       'audit.json':check_reference(reference,args['fasta'],args['curated']),
    }
    cases_seen={c['case_id'] for c in cases['cases']}
    if set(c['case_id'] for c in outputs['evidence.json']['records'])!=cases_seen:
        raise ValueError('Case identifiers do not match between variant and RNA evidence')
    root=Path(args['output_dir']);root.mkdir(parents=True,exist_ok=True)
    for name,table in outputs.items():
        (root/name).write_text(json.dumps(table,sort_keys=True,indent=2)+'\n')
if __name__=='__main__':main()
