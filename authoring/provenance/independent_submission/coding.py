#!/usr/bin/env python3
"""Reference implementation for a *synthetic* isoform-aware CDS editing task.

This is never copied into the agent image.
"""
import argparse
import hashlib
import json
from pathlib import Path

RC = str.maketrans('ACGT', 'TGCA')
def revcomp(s): return s.translate(RC)[::-1]

# Standard genetic code, unambiguous DNA only.
CODON_TABLE = '''
TTT F TTC F TTA L TTG L TCT S TCC S TCA S TCG S TAT Y TAC Y TAA * TAG * TGT C TGC C TGA * TGG W
CTT L CTC L CTA L CTG L CCT P CCC P CCA P CCG P CAT H CAC H CAA Q CAG Q CGT R CGC R CGA R CGG R
ATT I ATC I ATA I ATG M ACT T ACC T ACA T ACG T AAT N AAC N AAA K AAG K AGT S AGC S AGA R AGG R
GTT V GTC V GTA V GTG V GCT A GCC A GCA A GCG A GAT D GAC D GAA E GAG E GGT G GGC G GGA G GGG G
'''.split()
GENETIC = dict(zip(CODON_TABLE[::2], CODON_TABLE[1::2]))
assert len(GENETIC) == 64

def h(s): return hashlib.sha256(s.encode('ascii')).hexdigest()

def translated(seq):
    protein = []
    stop_index = None
    for i in range(0, len(seq)-2, 3):
        aa = GENETIC[seq[i:i+3]]
        if aa == '*':
            stop_index = len(protein)+1
            break
        protein.append(aa)
    return ''.join(protein), stop_index

def native_stop_retained(original, altered, edits):
    """Follow the original termination site's ancestry through CDS edits.

    A later stop sliding into its vacated codon index is *not* the native stop.
    In contrast, an upstream in-frame indel can shift an intact stop, and a
    same-length stop-to-stop substitution can preserve termination at its site.
    All edit coordinates refer to the unmodified transcript CDS.
    """
    _, ref_stop = translated(original)
    if ref_stop is None:
        return False
    start = (ref_stop-1)*3
    shift = 0
    for pos, ref, alt in sorted(edits):
        end = pos+len(ref)
        # Inserting within a triplet interrupts its provenance; an
        # unequal-length replacement overlapping the stop removes its site.
        if (not ref and start < pos < start+3) or (
                ref and pos < start+3 and end > start and len(ref)!=len(alt)):
            return False
        if pos <= start:
            shift += len(alt)-len(ref)
    target = start+shift
    if target < 0 or target % 3 or target+3 > len(altered):
        return False
    # Preserve ancestry even if a NEW upstream stop is encountered first.
    # Early termination and destruction of the ancestral stop are different events.
    return GENETIC[altered[target:target+3]] == '*'


def consequence(original, altered, had_exonic_edit, native_stop_intact=None):
    ref_pep, ref_stop = translated(original)
    alt_pep, alt_stop = translated(altered)
    premature = alt_stop is not None and alt_stop < len(altered)//3
    # An untraced direct call can still handle unchanged same-length stops.
    if native_stop_intact is None:
        native_stop_intact = bool(ref_stop is not None and alt_stop==ref_stop
                                  and len(original)==len(altered)
                                  and altered[(ref_stop-1)*3:ref_stop*3] in ('TAA','TAG','TGA'))
    if not had_exonic_edit or original == altered: klass = 'no_change'
    elif altered[:3] != 'ATG': klass = 'start_lost'
    elif (len(altered)-len(original)) % 3: klass = 'frameshift'
    elif premature and (ref_stop is None or alt_stop < ref_stop): klass = 'stop_gained'
    elif ref_stop is not None and not native_stop_intact: klass = 'stop_lost'
    elif len(altered)>len(original): klass = 'inframe_insertion'
    elif len(altered)<len(original): klass = 'inframe_deletion'
    elif alt_pep == ref_pep: klass = 'synonymous'
    else: klass = 'protein_altering'
    first = next((i+1 for i,(a,b) in enumerate(zip(ref_pep,alt_pep)) if a!=b), None)
    if first is None and len(ref_pep)!=len(alt_pep): first=min(len(ref_pep),len(alt_pep))+1
    return {
      'consequence':klass,'reference_cds_sha256':h(original), 'altered_cds_sha256':h(altered),
      'reference_protein_sha256':h(ref_pep), 'altered_protein_sha256':h(alt_pep),
      'reference_cds_length':len(original),'altered_cds_length':len(altered),
      'reference_protein_length':len(ref_pep),'altered_protein_length':len(alt_pep),
      'first_changed_residue':first,'premature_stop':premature,'stop_observed':alt_stop is not None,
    }

def analyze_variant_set(case, transcript, genome, strand):
    exons = transcript['exons']
    ref_cds = ''.join(revcomp(genome[s:e]) if strand=='-' else genome[s:e] for s,e in exons)
    variants = case['variants']
    invalid = set()
    projected = []
    # Detect invalid coordinates and genomic reference errors independent of transcript inclusion.
    for v in variants:
        start=v.get('start'); ref=v.get('ref'); alt=v.get('alt')
        if (type(start) is not int or start<0 or not isinstance(ref,str) or not isinstance(alt,str)
             or len(ref)==len(alt)==0 or not all(b in 'ACGT' for b in ref+alt)
             or start+len(ref)>len(genome)):
            invalid.add('invalid_edit');continue
        if genome[start:start+len(ref)]!=ref:invalid.add('reference_mismatch')
    valid=[v for v in variants if type(v.get('start')) is int and isinstance(v.get('ref'),str)
           and isinstance(v.get('alt'),str) and v['start']>=0]
    # Reject identical or overlapping non-empty genomic variants within the same phased set.
    for i, a in enumerate(valid):
        for b in valid[i+1:]:
            a1,a2=a['start'],a['start']+len(a['ref'])
            b1,b2=b['start'],b['start']+len(b['ref'])
            if (a1==b1 or (a1<b2 and b1<a2)):
                invalid.add('overlapping_edits')
    for v in valid:
        start=v['start']; finish=start+len(v['ref']); fl=[]
        # A single deletion/replacement that spans multiple exons is ambiguous
        # because the edit also traverses intronic sequence.
        if len(v['ref']) and sum(start < e and finish > s for s, e in exons) > 1:
            invalid.add('multiple_exon_overlap')
        tx_offset=0
        for exon_start, exon_end in exons:
            if len(v['ref']) == 0:
                if exon_start < start < exon_end:
                    tx_pos = tx_offset + (exon_end-start if strand=='-' else start-exon_start)
                    fl.append((tx_pos, '', revcomp(v['alt']) if strand=='-' else v['alt']))
                elif start in (exon_start,exon_end):
                    invalid.add('exon_boundary_ambiguity')
            else:
                overlap=start<exon_end and finish>exon_start
                if overlap:
                    if start < exon_start or finish > exon_end:
                        invalid.add('partial_exon_overlap')
                    else:
                        tx_pos = tx_offset+(exon_end-finish if strand=='-' else start-exon_start)
                        fl.append((tx_pos,revcomp(v['ref']) if strand=='-' else v['ref'],revcomp(v['alt']) if strand=='-' else v['alt']))
            tx_offset += exon_end-exon_start
        if len(fl)>1:invalid.add('multiple_exon_overlap')
        projected+=fl
    if invalid:
        return {'case_id':case['case_id'],'transcript_id':transcript['transcript_id'],
            'status':'unresolved','reason_codes':sorted(invalid),'product':None}
    # Optional independent transcript sequence anchor. The exonic reference is
    # projected from the genome; a discordant supplied CDS anchor is a genuine
    # projection discrepancy, but only otherwise valid exonic edits carry it.
    cds_anchor=transcript.get('cds_sequence')
    if projected and cds_anchor is not None and cds_anchor != ref_cds:
        return {'case_id':case['case_id'],'transcript_id':transcript['transcript_id'],
            'status':'unresolved','reason_codes':['projection_mismatch'],'product':None}
    changes=[]
    for pos,ref,alt in projected:
        if ref_cds[pos:pos+len(ref)] != ref:
            return {'case_id':case['case_id'],'transcript_id':transcript['transcript_id'],
                'status':'unresolved','reason_codes':['projection_mismatch'],'product':None}
        changes.append((pos,ref,alt))
    new_cds=ref_cds
    for pos,ref,alt in sorted(changes,reverse=True):
        new_cds=new_cds[:pos]+alt+new_cds[pos+len(ref):]
    return {'case_id':case['case_id'],'transcript_id':transcript['transcript_id'],
            'status':'resolved','reason_codes':[],
            'product':consequence(ref_cds,new_cds,bool(changes),native_stop_retained(ref_cds,new_cds,changes))}

def compute(reference, data):
    genome=reference['genome']; transcripts=reference['transcripts']
    records=[]
    for case in data['cases']:
        for tx in transcripts:
            records.append(analyze_variant_set(case, tx, genome, reference['strand']))
    return {'reference_id':reference['reference_id'],
            'records':sorted(records,key=lambda d:(d['case_id'],d['transcript_id']))}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--reference',required=True)
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    output=compute(json.loads(Path(a.reference).read_text()),json.loads(Path(a.input).read_text()))
    Path(a.output).write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
