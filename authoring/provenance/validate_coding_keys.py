#!/usr/bin/env python3
"""Provenance for expected coding keys: reference-oracle generation with an independent genomic reconstruction.

The independent check reconstructs mutated exon DNA on the genomic plus strand,
then reverse complements each exon and translates until the first stop. It
neither imports nor calls the reference CDS-edit projection implementation.
"""
import hashlib, importlib.util, json
from pathlib import Path
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parents[2]
def sha(s): return hashlib.sha256(s.encode('ascii')).hexdigest()
def load(path):return json.loads(path.read_text())
spec=importlib.util.spec_from_file_location('source_oracle',ROOT/'solution/coding.py')
oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
ref=load(ROOT/'environment/data/reference.json');genome=ref['genome']; report={}
for tag in ('sealed','hidden'):
    cases=load(ROOT/'tests'/f'{tag}_cases.json');expected=load(ROOT/'tests'/f'{tag}_expected.json')
    actual=oracle.compute(ref,cases)
    assert actual==expected, f'{tag} expected coding JSON differs from the named producer'
    by_case={c['case_id']:c for c in cases['cases']};validated=0;unresolved=0
    for record in expected['records']:
        case=by_case[record['case_id']]
        tx=next(x for x in ref['transcripts'] if x['transcript_id']==record['transcript_id'])
        if record['status']!='resolved':
            assert record['product'] is None
            unresolved+=1;continue
        exons=[]
        for a,b in tx['exons']:
            dna=genome[a:b]
            edits=[]
            for v in case['variants']:
                start=v['start'];length=len(v['ref'])
                if (length==0 and a<start<b) or (length>0 and a<=start and start+length<=b):
                    edits.append(v)
            for edit in sorted(edits,key=lambda e:e['start'],reverse=True):
                x=edit['start']-a;dna=dna[:x]+edit['alt']+dna[x+len(edit['ref']):]
            # Negative strand: genomic edit first, then independently assemble reverse complement.
            exons.append(str(Seq(dna).reverse_complement()) if ref['strand']=='-' else dna)
        altered=''.join(exons)
        original=''.join(str(Seq(genome[a:b]).reverse_complement()) if ref['strand']=='-' else genome[a:b]
                         for a,b in tx['exons'])
        def tr(seq):
            from Bio.Data.CodonTable import unambiguous_dna_by_id
            code=unambiguous_dna_by_id[1]
            protein=[]
            for idx in range(0,len(seq)-2,3):
                cod=seq[idx:idx+3]
                if cod in code.stop_codons: break
                protein.append(code.forward_table[cod])
            return ''.join(protein)
        product=record['product'];checks={
            'reference_cds_sha256':sha(original), 'altered_cds_sha256':sha(altered),
            'reference_protein_sha256':sha(tr(original)), 'altered_protein_sha256':sha(tr(altered)),
            'reference_cds_length':len(original),'altered_cds_length':len(altered),
            'reference_protein_length':len(tr(original)),'altered_protein_length':len(tr(altered))}
        for field,value in checks.items():
            assert product[field]==value, (tag,record['case_id'],record['transcript_id'],field)
        validated+=1
    report[tag]={'expected_file':f'tests/{tag}_expected.json',
                 'sha256':hashlib.sha256((ROOT/'tests'/f'{tag}_expected.json').read_bytes()).hexdigest(),
                 'source':'solution/coding.py:compute() (reference-produced, then persisted)',
                 'independent_resolved_checks':validated,'unresolved_records_oracle':unresolved,
                 'checks':'Genomic-plus exon edit/reverse complement; coding and peptide SHA-256/length separately reconstructed with Biopython codon table'}
path=ROOT/'authoring/evidence/coding_key_provenance_2026-11-09.json'
path.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2))
