"""Sequence-level checks anchored to a real accessioned human protein sequence."""
from hashlib import sha256
from pathlib import Path

REV=str.maketrans('ACGT','TGCA')

def check_reference(reference, fasta_path, curated_path):
    contents=Path(fasta_path).read_text().splitlines()
    if not contents or not contents[0].startswith('>NP_001390.1'):
        raise ValueError('Wrong accession for protein anchor')
    authentic=''.join(row.strip() for row in contents[1:])
    from coding import translated
    products={}
    for tr in reference['transcripts']:
        cds=''.join(reference['genome'][s:e].translate(REV)[::-1] for s,e in tr['exons'])
        products[tr['protein_model']]=translated(cds)[0]
    expected_A2=authentic[:306]+authentic[308:]
    if products.get('A1')!=authentic or products.get('A2')!=expected_A2:
        raise ValueError('Protein anchor failed')
    import json
    panel=json.loads(Path(curated_path).read_text())
    curated=public_variant_concordance(panel,authentic)
    return {'curated_concordance':curated,'accession':'NP_001390.1','a1_length':len(authentic),'a2_length':len(expected_A2),
            'a1_sha256':sha256(authentic.encode()).hexdigest(),
            'a2_sha256':sha256(expected_A2.encode()).hexdigest(),
            'missing_a2_residues':'Val307-Glu308',
            'coordinates_are_synthetic':True}


# Publicly reported EDA protein substitutions are audited against the same
# versioned protein references used by the simulated transcript workbench.
# This is variant consequence reconciliation, not a diagnostic classifier.
THREE_TO_ONE = dict(zip(
    ('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val').split(),
    'ARNDCQEGHILKMFPSTWYV'))


def public_variant_concordance(panel, authentic):
    import re
    if panel['protein_anchor'] != 'NP_001390.1' or panel['transcript_anchor'] != 'NM_001399.5':
        raise ValueError('Public variant panel has an unsupported reference version')
    a2 = authentic[:306] + authentic[308:]
    rows = []
    for v in panel['variants']:
        m = re.fullmatch(r'p\.([A-Z][a-z]{2})([1-9][0-9]*)([A-Z][a-z]{2})', v['protein_change'])
        if not m or m.group(1) not in THREE_TO_ONE or m.group(3) not in THREE_TO_ONE:
            raise ValueError('Unsupported protein-change notation')
        before, number, after = THREE_TO_ONE[m.group(1)], int(m.group(2)), THREE_TO_ONE[m.group(3)]
        observations = v['classification_observations']
        assertions = sorted({o['classification'] for o in observations})
        accessions = sorted({o['record'] for o in observations})
        if not assertions or not accessions:
            raise ValueError('Missing evidence provenance')
        mapped={}
        for iso, peptide in [('A1',authentic),('A2',a2)]:
            if iso == 'A2' and number in (307,308):
                mapped[iso]={'status':'isoform_absent','position':None,'reference_residue':None,
                             'alternate_residue':None,'altered_protein_sha256':None}
                continue
            pos = number if iso=='A1' or number<=306 else number-2
            if pos < 1 or pos > len(peptide) or peptide[pos-1] != before:
                mapped[iso]={'status':'reference_mismatch','position':pos,'reference_residue':peptide[pos-1] if 1<=pos<=len(peptide) else None,
                             'alternate_residue':after,'altered_protein_sha256':None}
                continue
            altered=peptide[:pos-1]+after+peptide[pos:]
            mapped[iso]={'status':'mapped','position':pos,'reference_residue':before,
                         'alternate_residue':after,'altered_protein_sha256':sha256(altered.encode()).hexdigest()}
        rows.append({'variant_id':v['variant_id'],'protein_change':v['protein_change'],
                     'source_accessions':accessions,'reported_assertions':assertions,
                     'assertion_concordance':'discordant' if len(assertions)>1 else 'single_class',
                     'isoforms':mapped})
    if len({r['variant_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate public variant ID')
    # Curated, assembly-specific evidence is kept separate from the artificial
    # negative-strand sequence simulator. No chromosome alleles are validated
    # without a vendored authentic genomic sequence.
    context=panel['genome_context']
    chrom=context['chromosome_accession']
    begin,end=context['gene_interval_1based_closed']
    if not (isinstance(begin,int) and isinstance(end,int) and begin<end):
        raise ValueError('Invalid published gene boundaries')
    if context['gene_orientation'] not in ('+','-'):
        raise ValueError('Invalid published strand')
    if context['mane_select_transcript']!='NM_001399.5' or context['a2_transcript']!='NM_001005609.2':
        raise ValueError('Unsupported transcript accessions')
    crosswalk=[]
    known={row['variant_id']:row for row in rows}
    for item in panel['genomic_crosswalk_records']:
        match=re.fullmatch(r'([^:]+):g\.([1-9][0-9]*)([ACGT])>([ACGT])',item['genomic_hgvs'])
        if not match or match.group(1)!=chrom:
            raise ValueError('Genomic description must use the declared chromosome reference')
        pos=int(match.group(2)); inside=begin<=pos<=end
        if item['variant_id'] not in known or item['protein_change'] != known[item['variant_id']]['protein_change']:
            raise ValueError('Genomic description does not match curated protein variant')
        if not all(item['transcript_hgvs'].get(iso,'').startswith(
            context['mane_select_transcript' if iso=='A1' else 'a2_transcript']+':c.')
                   for iso in ('A1','A2')):
            raise ValueError('Transcript crosswalk lacks an accession match')
        if len(item['transcript_hgvs'])!=2:
            raise ValueError('Incorrect isoform set')
        # Validate transcript HGVS against the protein substitution, rather
        # than accepting only the accession prefix. No genuine chromosome
        # sequence is bundled: this validates a *possible* codon-level mapping,
        # not the actual alleles of the patient's chromosome.
        from coding import GENETIC
        transcript_positions={}
        protein_match=re.fullmatch(r'p\.([A-Z][a-z]{2})([1-9][0-9]*)([A-Z][a-z]{2})',known[item['variant_id']]['protein_change'])
        if protein_match is None:
            raise ValueError('Invalid protein substitution')
        for iso in ('A1','A2'):
            c_hgvs=item['transcript_hgvs'][iso]
            cm=re.fullmatch(r'([^:]+):c\.([1-9][0-9]*)([ACGT])>([ACGT])',c_hgvs)
            if cm is None:
                raise ValueError('Unusable transcript HGVS substitution')
            cdna_pos=int(cm.group(2));ref_base=cm.group(3);alt_base=cm.group(4)
            expected_iso=known[item['variant_id']]['isoforms'][iso]
            protein_pos=expected_iso['position']
            if expected_iso['status']!='mapped' or protein_pos is None or (cdna_pos-1)//3+1!=protein_pos:
                raise ValueError('Transcript coding coordinate does not map to the protein residue')
            # For the positive-strand public locus, genomic and cDNA SNV
            # alleles must match. A negative-strand mapping uses complements.
            gref,galt=match.group(3),match.group(4)
            if context['gene_orientation']=='-':
                complement=str.maketrans('ACGT','TGCA')
                gref,galt=gref.translate(complement),galt.translate(complement)
            if (ref_base,alt_base)!=(gref,galt):
                raise ValueError('Transcript and genomic allele changes disagree')
            ref_aa,alt_aa=expected_iso['reference_residue'],expected_iso['alternate_residue']
            codon_offset=(cdna_pos-1)%3
            if not any(GENETIC[codon]==ref_aa and codon[codon_offset]==ref_base and
                       GENETIC.get(codon[:codon_offset]+alt_base+codon[codon_offset+1:])==alt_aa
                       for codon in GENETIC):
                raise ValueError('Transcript nucleotide substitution cannot yield reported protein change')
        matched=known[item['variant_id']]['isoforms']
        crosswalk.append({
            'variant_id':item['variant_id'], 'genomic_hgvs':item['genomic_hgvs'],
            'in_annotated_gene':inside,
            'transcript_hgvs':dict(sorted(item['transcript_hgvs'].items())),
            'isoform_status':{iso:matched[iso]['status'] for iso in ('A1','A2')},
            'clinical_record':item['clinical_record'],
            'genomic_allele_sequence_validated':False
        })
    crosswalk=sorted(crosswalk,key=lambda x:x['variant_id'])
    if len({x['variant_id'] for x in crosswalk})!=len(crosswalk):
        raise ValueError('Duplicate genomic variant IDs')
    genome_report={
        'assembly':context['assembly'], 'chromosome_accession':chrom,
        'gene_orientation':context['gene_orientation'],
        'gene_interval_1based_closed':[begin,end],
        'simulated_workbench_orientation':'-',
        'actual_chromosome_sequence_supplied':False,
        'records':crosswalk,
    }
    return {'panel_id':panel['panel_id'],'genomic_crosswalk':genome_report,'source_type':panel['source_type'],
            'clinical_interpretation_inferred':False,
            'variants':sorted(rows,key=lambda r:r['variant_id'])}
