#!/usr/bin/env python3
"""Deterministic EDA-A1/A2-*inspired* synthetic locus. NOT NCBI sequences."""
from pathlib import Path
import json, random, hashlib
R = Path(__file__).resolve().parents[2]
random.seed(13052026)
CODONS = ['GCT','GCC','GCA','GCG','CTG','TTG','GTT','GTC','CCT','CCG','TCT','AGC','ACT','ACC','AAC','AAT','GAA','GAG','GAT','GAC','AAA','AAG','AGA','AGG','TAC','TAT','TGG','CAT','CAC','ATC','ATT','ATG','GGT','GGC','CAG','CAA','TTC','TTT','CGC','CGT']
COMP=str.maketrans('ACGT','TGCA')
rc=lambda s:s.translate(COMP)[::-1]
# Protein-anchored construct. The protein is NP_001390.1 (391 aa), not randomly generated.
# Nucleotide codons and splice coordinate intervals remain artificial; this is NOT NM_001399.5.
# UniProt Q92838-3 (EDA-A2) has Val307-Glu308 deleted relative to EDA-A1.
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Seq import Seq
protein = ''.join(s.strip() for s in (R/'environment/data/NP_001390.1.faa').read_text().splitlines() if not s.startswith('>'))
assert len(protein)==391 and protein[306:308]=='VE' and protein[0]=='M'
codons_by_aa={}
for codon,aa in unambiguous_dna_by_id[1].forward_table.items():
    codons_by_aa.setdefault(aa,[]).append(codon)
A1=''.join(random.choice(codons_by_aa[aa]) for aa in protein)+'TAA'
A2=A1[:306*3]+A1[308*3:]
assert len(A1)==1176 and len(A2)==1170
assert str(Seq(A1[:-3]).translate())==protein
assert str(Seq(A2[:-3]).translate())==protein[:306]+protein[308:]
# Six artificial coding exons; A1-specific VE triplets are last six bases of exon 5.
lengths=[180,225,213,171,135,252]
coords=[]
starts=[10000,8200,6550,4890,3200,1420]
for s,L in zip(starts,lengths): coords.append([s,s+L])
contig=list(''.join(random.choice('ACGT') for _ in range(12000)))
for (s,e),cds_segment in zip(coords,[A1[sum(lengths[:i]):sum(lengths[:i+1])] for i in range(6)]):
    contig[s:e]=rc(cds_segment)
contig=''.join(contig)
ref={
 'reference_id':'EDA_SIM_NEGSTRAND_2026_1','synthetic':True,'contig_name':'chrSimX','strand':'-',
 'genome':contig,
 'transcripts':[
    {'transcript_id':'EDA_SIM_A1.1','protein_model':'A1','exons':coords,'cds_start':0,'cds_end':1176},
    {'transcript_id':'EDA_SIM_A2.1','protein_model':'A2','exons':[ *coords[:4], [coords[4][0]+6,coords[4][1]],coords[5]],'cds_start':0,'cds_end':1170}
 ],'note':'EDA-A1 and EDA-A2 polypeptide identities are anchored to NCBI NP_001390.1 and UniProt Q92838 A2 deletion; DNA codons, chromosome coordinates, exon structure and cases are simulated and are not clinical data.'
}
# verify assembled transcript sequences
for t, seq in zip(ref['transcripts'],[A1,A2]):
    assembled=''.join(rc(contig[a:b]) for a,b in t['exons'])
    assert assembled==seq and assembled.startswith('ATG') and assembled.endswith('TAA')
(R/'environment/data/reference.json').write_text(json.dumps(ref,separators=(',',':'))+'\n')

def variant(pos,L,alt=None):
    v={'start':int(pos),'ref':contig[pos:pos+L]}
    if alt is None: alt=rc('A'*(L+1))
    v['alt']=alt
    return v

def snv(pos):
    n=contig[pos]
    return variant(pos,1,{'A':'C','C':'G','G':'T','T':'A'}[n])

def insertion(pos,nt):
    return {'start':pos,'ref':'','alt':rc(nt)}

def deletion(pos,n):return {'start':pos,'ref':contig[pos:pos+n],'alt':''}

def coding_offset_to_genomic(cpos, iso='A1'):
    t=ref['transcripts'][0 if iso=='A1' else 1]; offset=0
    for a,b in t['exons']:
        if offset <= cpos < offset+(b-a):return b-1-(cpos-offset)
        offset+=b-a
    raise ValueError(cpos)

def phased_case(case_id, variants, label):
    # Human-readable consequence labels belong to authoring, never to
    # the graded variant inputs: they would disclose intended outputs.
    return {'case_id':case_id,'variants':variants}

base=[
 phased_case('m01',[snv(coding_offset_to_genomic(39))],'one_snv'),
 phased_case('m02',[snv(coding_offset_to_genomic(61))],'early_snv'),
 phased_case('m03',[insertion(coding_offset_to_genomic(300),'G')],'frameshift_insertion'),
 phased_case('m04',[insertion(coding_offset_to_genomic(245),'AAG')],'inframe_insertion'),
 phased_case('m05',[deletion(coding_offset_to_genomic(470)-2,3)],'inframe_deletion'),
 phased_case('m06',[snv(coding_offset_to_genomic(0))],'initiation_change'),
 phased_case('m07',[snv(coding_offset_to_genomic(105)),snv(coding_offset_to_genomic(106))],'cis_same_codon'),
 phased_case('m08',[snv(coding_offset_to_genomic(919))],'a1_specific'),
 phased_case('m09',[deletion(coords[4][0]+3,6)],'isoform_boundary'),
 phased_case('m10',[deletion(coords[3][0]+22,4),snv(coords[3][0]+24)],'overlapping_variants'),
 phased_case('m11',[variant(coords[0][0]+4,1,'N')],'non_acgt'),
 phased_case('m12',[variant(11200,1,'G')],'intronic_variant'),
 phased_case('m13',[{'start':coords[0][0]+10,'ref':'G','alt':'A'}],'reference_mismatch_possible'),
 phased_case('m14',[deletion(coords[2][0]-2,7)],'spanning_exon_boundary'),
 phased_case('m15',[snv(coding_offset_to_genomic(532)),snv(coding_offset_to_genomic(877))],'independent_exons'),
 phased_case('m16',[deletion(coding_offset_to_genomic(600)-2,2)],'frameshift_deletion'),
 phased_case('m17',[insertion(coding_offset_to_genomic(677),'C'),snv(coding_offset_to_genomic(302))],'compound_frameshift'),
 phased_case('m18',[insertion(coords[4][0],'GTT')],'boundary_insertion'),
 phased_case('m19',[snv(coding_offset_to_genomic(920)),snv(coding_offset_to_genomic(940)),snv(coding_offset_to_genomic(960))],'three_edits'),
 phased_case('m20',[deletion(coding_offset_to_genomic(1169),3)],'terminal_edit')
]
# Target scientifically distinct protein consequences, not accidental cryptic stops.
base[4]=phased_case('m05',[deletion(coding_offset_to_genomic(470),3)],'codon_aligned_deletion')
base[19]=phased_case('m20',[snv(coding_offset_to_genomic(1175))],'stop_lost')
# Synonymous SNV in a synonymous four-fold codon, selected by translation.
from Bio.Seq import Seq
for c in range(90,600,3):
    orig=A1[c:c+3]
    if c<=918<=c+3: continue
    for x in 'ACGT':
        altcodon=orig[:2]+x
        if x!=orig[2] and str(Seq(orig).translate())==str(Seq(altcodon).translate()):
            p=c+2;g=coding_offset_to_genomic(p)
            snv_ref=contig[g]
            snv_alt={'A':'C','C':'G','G':'T','T':'A'}[snv_ref]
            desired_plus={'A':'T','C':'G','G':'C','T':'A'}[x]
            base.append(phased_case('m21',[variant(g,1,desired_plus)],'synonymous_codon'))
            break
    else: continue
    break
else: raise AssertionError('No synonym selected')
base.append(phased_case('m22',[insertion(11300,'GT'),insertion(coding_offset_to_genomic(1051),'ACG')],'intronic_plus_coding_insertion'))
base.append(phased_case('m23',[insertion(coding_offset_to_genomic(219),'GCA'),snv(coding_offset_to_genomic(747))],'phased_compound'))
# Ensure valid non-exonic ref variant; and one deliberate wrong ref.
base[11]['variants']=[snv(11200),insertion(coding_offset_to_genomic(1012),'GG')]
base[12]['variants'][0]={'start':coords[0][0]+10,'ref':{'A':'C','C':'G','G':'T','T':'A'}[contig[coords[0][0]+10]],'alt':'A'}
# Post-review independent compound cases: eliminate accidental no-variant
# identity matches in A2 for m08 and m18 while preserving the biologically
# meaningful isoform-specific boundary situations in the first component.
base[7]['variants'].append(dict(base[0]['variants'][0]))
# Case labels are deliberately not exported to agent-visible variant JSON.
base[17]['variants'].append(dict(base[1]['variants'][0]))
# Case labels are deliberately not exported to agent-visible variant JSON.
# early public practice, private sealed, plus unseen challenge cases
(R/'environment/data/practice.json').write_text(json.dumps({'cases':base[:6]},indent=2)+'\n')
(R/'tests/sealed_cases.json').write_text(json.dumps({'cases':base},separators=(',',':'))+'\n')
# Synthetic seeded additional cases covering diverse locus positions.
rng=random.Random(17831)
extra=[]
for i in range(28):
    p=rng.randrange(0,1170)
    start=coding_offset_to_genomic(p)
    kind=i%7
    if kind==0: v=[snv(start)]
    elif kind==1: v=[insertion(start,'G')]
    elif kind==2: v=[insertion(start,'ACG')]
    elif kind==3: v=[deletion(start,2)]
    elif kind==4: v=[snv(start),snv(coding_offset_to_genomic((p+42)%1169))]
    elif kind==5: v=[insertion(start,'CTAGTA')]
    else: v=[deletion(start,1)]
    extra.append(phased_case('h%02d'%i,v,'seeded_hidden'))
(R/'tests/hidden_cases.json').write_text(json.dumps({'cases':extra},separators=(',',':'))+'\n')
(R/'authoring/provenance/source_manifest.json').write_text(json.dumps({
 'is_synthetic':True,'transcript_labels':['EDA_SIM_A1.1','EDA_SIM_A2.1'],
 'biological_analogy':'Authentic NP_001390.1 A1 polypeptide and the A2 loss of Val307-Glu308, backtranslated into synthetic synonymous codons and artificial negative-strand exons' ,
 'human_accessions_for_context_only':['NM_001399.5','NM_001005609.2'],
 'fixture_script':'authoring/provenance/generate_fixture.py',
 'seed_primary':13052026,'seed_hidden':17831,
 'reference_sha256':hashlib.sha256((R/'environment/data/reference.json').read_bytes()).hexdigest(),
 'not_patient_data':True, 'protein_reference':'NP_001390.1', 'protein_provenance':'user-provided NCBI FASTA checked against UniProt Q92838 canonical 391 aa; residue 307-308 VE; nucleotide model synthetic'
},indent=2)+'\n')
print('Generated',len(base), 'sealed cases +',len(extra),'hidden cases')
