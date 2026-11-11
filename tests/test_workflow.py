"""Sealed scientific validation. The verifier never imports agent-authored modules."""
import json,os,signal,stat,shutil,subprocess,sys,tempfile,hashlib,ast,math,copy,time
from decimal import Decimal
from pathlib import Path
import pytest
from process_hygiene import clean_agent_processes

ROOT=Path(__file__).resolve().parent
SOURCE=Path(os.getenv('EDA_WORKBENCH_ARTIFACT','/app/submission'))
CODE_NAMES=['cli.py','coding.py','evidence.py','audit.py','inference.py']

class WorkflowNonzeroExit(AssertionError):
    """The submitted program explicitly rejected an invalid input.

    Unlike generic AssertionError, this excludes timeout, missing output,
    verifier assertion, or process-cleanup failures from rejection controls.
    """

def load_exact(p):
    def pairs_unique(pairs):
        d={}
        for key,val in pairs:
            if key in d: raise ValueError('Duplicate JSON property')
            d[key]=val
        return d
    return json.loads(Path(p).read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)),object_pairs_hook=pairs_unique)

def run_agent(cases,observations,calibration_override=None,curated_override=None,reference_override=None,fasta_override=None,batch_assay_override=None,phased_assay_override=None,longitudinal_assay_override=None,phase_graph_override=None):
    assert SOURCE.is_dir() and not SOURCE.is_symlink()
    assert all((SOURCE/name).is_file() and not (SOURCE/name).is_symlink() for name in CODE_NAMES)
    assert all(0<(SOURCE/name).stat().st_size<100_000 for name in CODE_NAMES)
    assert sum((SOURCE/name).stat().st_size for name in CODE_NAMES)<=400_000
    with tempfile.TemporaryDirectory(prefix='eda_validate_') as tmp:
        scratch=Path(tmp)
        os.chown(scratch,65534,65534)
        os.chmod(scratch,0o700)
        code=scratch/'workbench';code.mkdir();code.chmod(0o755)
        for name in CODE_NAMES:
            dst=code/name
            shutil.copyfile(SOURCE/name,dst)
            dst.chmod(0o644)
        ref=scratch/'reference.json';ref.write_text(json.dumps(reference_override) if reference_override is not None else (ROOT/'assets/reference.json').read_text());ref.chmod(0o644)
        fasta=scratch/'NP_001390.1.faa';fasta.write_text(fasta_override if fasta_override is not None else (ROOT/'assets/NP_001390.1.faa').read_text());fasta.chmod(0o644)
        curated=scratch/'public_eda_variants.json';curated.write_text(json.dumps(curated_override) if curated_override is not None else (ROOT/'assets/public_eda_variants.json').read_text());curated.chmod(0o644)
        cal=scratch/'spikein_calibration.json';cal.write_text(json.dumps(calibration_override) if calibration_override is not None else (ROOT/'assets/spikein_calibration.json').read_text());cal.chmod(0o644)
        data=scratch/'cases.json';data.write_text(json.dumps(cases));data.chmod(0o644)
        obs=scratch/'molecules.json';obs.write_text(json.dumps(observations));obs.chmod(0o644)
        ids={case['case_id'] for case in cases['cases']}
        src='hidden_batch_assay.json' if any(str(cid).startswith('h') for cid in ids) else 'sealed_batch_assay.json'
        assay_data=copy.deepcopy(batch_assay_override) if batch_assay_override is not None else load_exact(ROOT/src)
        known={r['case_id'] for r in assay_data['measurements']}
        missing=ids-known
        # Newly constructed boundary probes must still exercise the batch stage;
        # use replicated public-fixture rows rather than weakening the case-ID contract.
        if missing and batch_assay_override is None:
            first=sorted(known)[0]
            exemplar=[r for r in assay_data['measurements'] if r['case_id']==first]
            for cid in sorted(missing):
                assay_data['measurements'].extend({**r,'case_id':cid} for r in exemplar)
        if batch_assay_override is None:
            assay_data['measurements']=[r for r in assay_data['measurements'] if r['case_id'] in ids]
        assay=scratch/'batch_assay.json';assay.write_text(json.dumps(assay_data));assay.chmod(0o644)
        src_phase='hidden_phased_assay.json' if any(str(cid).startswith('h') for cid in ids) else 'sealed_phased_assay.json'
        phased_data=copy.deepcopy(phased_assay_override) if phased_assay_override is not None else load_exact(ROOT/src_phase)
        covered={r['case_id'] for r in phased_data['measurements']}
        missing=ids-covered
        if missing and phased_assay_override is None:
            exemplar=sorted(covered)[0]
            original=[r for r in phased_data['measurements'] if r['case_id']==exemplar]
            for cid in sorted(missing):
                phased_data['measurements'].extend({**r,'case_id':cid} for r in original)
        if phased_assay_override is None:
            phased_data['measurements']=[r for r in phased_data['measurements'] if r['case_id'] in ids]
        phased=scratch/'phased_assay.json';phased.write_text(json.dumps(phased_data));phased.chmod(0o644)
        long_src='hidden_longitudinal_assay.json' if any(str(cid).startswith('h') for cid in ids) else 'sealed_longitudinal_assay.json'
        long_data=copy.deepcopy(longitudinal_assay_override) if longitudinal_assay_override is not None else load_exact(ROOT/long_src)
        coverage={r['case_id'] for r in long_data['measurements']}
        missing=ids-coverage
        if missing and longitudinal_assay_override is None:
            original=[r for r in long_data['measurements'] if r['case_id']==sorted(coverage)[0]]
            for cid in sorted(missing):
                long_data['measurements'].extend({**r,'case_id':cid} for r in original)
        if longitudinal_assay_override is None:
            long_data['measurements']=[r for r in long_data['measurements'] if r['case_id'] in ids]
        long_path=scratch/'longitudinal_assay.json';long_path.write_text(json.dumps(long_data));long_path.chmod(0o644)

        graph_src='hidden_phase_graph_assay.json' if any(str(cid).startswith('h') for cid in ids) else 'sealed_phase_graph_assay.json'
        graph_data=copy.deepcopy(phase_graph_override) if phase_graph_override is not None else load_exact(ROOT/graph_src)
        known={r['case_id'] for r in graph_data['cases']}
        missing=ids-known
        if missing and phase_graph_override is None:
            source_row=next(r for r in graph_data['cases'] if r['case_id']==sorted(known)[0])
            for cid in sorted(missing):
                added=copy.deepcopy(source_row)
                added['case_id']=cid
                for edge in added['edges']:
                    edge['replicate_id']=f'{cid}_' + edge['replicate_id']
                graph_data['cases'].append(added)
        if phase_graph_override is None:
            graph_data['cases']=[row for row in graph_data['cases'] if row['case_id'] in ids]
        graph_path=scratch/'phase_graph_assay.json';graph_path.write_text(json.dumps(graph_data));graph_path.chmod(0o644)

        out=scratch/'results';out.mkdir();os.chown(out,65534,65534);out.chmod(0o700)
        stdout=scratch/'stdout.txt';stderr=scratch/'stderr.txt'
        def demote():
            os.setgroups([]);os.setgid(65534);os.setuid(65534)
        cmd=[sys.executable,'-I','-S','-B',str(code/'cli.py'),'--reference',str(ref),
            '--variants',str(data),'--observations',str(obs),'--fasta',str(fasta),'--calibration',str(cal),'--curated',str(curated),'--batch-assay',str(assay),'--phased-assay',str(phased),'--longitudinal-assay',str(long_path),'--phase-graph',str(graph_path),'--output-dir',str(out)]
        invocation_started=time.monotonic()
        with stdout.open('wb') as so,stderr.open('wb') as se:
            proc=subprocess.Popen(cmd,cwd=tmp,
                env={'PATH':'/usr/local/bin:/usr/bin:/bin','HOME':tmp,'PYTHONDONTWRITEBYTECODE':'1'},
                stdout=so,stderr=se,preexec_fn=demote,start_new_session=True)
            try:
                try:proc.wait(timeout=35)
                except subprocess.TimeoutExpired:raise AssertionError('Workflow exceeded 35 seconds')
            finally:
                try:os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                proc.wait()
                # Count actual solver launches and their elapsed time across all
                # isolated verifier groups. This trusted file is not agent-writable.
                run_dir=os.environ.get('EDA_VERIFIER_RUN_DIR')
                if run_dir:
                    with (Path(run_dir)/'invocations_seconds.log').open('a') as accounting:
                        accounting.write(f'{time.monotonic()-invocation_started:.8f}\n')
                # Process groups alone miss detached grandchildren (setsid).
                # Root-owned verification tests must leave no live untrusted
                # descendants before any output files are examined.
                assert clean_agent_processes(), 'Untrusted child processes survived cleanup'
        if proc.returncode != 0:
            raise WorkflowNonzeroExit('Workflow exited nonzero: '+stderr.read_text(errors='replace')[:800])
        result={}
        for name in ['coding.json','evidence.json','audit.json','likelihood.json','batch_inference.json','phased_rna.json','longitudinal.json','phase_graph.json']:
            path=out/name
            assert path.is_file() and not path.is_symlink(),f'Missing regular output {name}'
            assert 0<path.stat().st_size<500_000,f'Oversized or empty {name}'
            result[name]=load_exact(path)
        assert set(out.iterdir())=={out/n for n in result},'Unexpected output files'
        return result

@pytest.fixture(scope='session')
def inputs():
    return (load_exact(ROOT/'sealed_cases.json'),load_exact(ROOT/'sealed_observations.json'),
            load_exact(ROOT/'hidden_cases.json'),load_exact(ROOT/'hidden_observations.json'))

@pytest.fixture(scope='session')
def outputs(inputs):
    # Reuse baseline calls across isolated pytest groups within ONE verifier
    # attempt. The runner creates a fresh root-only cache every attempt;
    # agent UID 65534 cannot alter it. All scientific checks remain intact.
    run_dir=os.environ.get('EDA_VERIFIER_RUN_DIR')
    cache=Path(run_dir)/'baseline_outputs.json' if run_dir else None
    if cache is not None and cache.exists():
        stored=load_exact(cache)
        assert type(stored) is list and len(stored)==2
        return tuple(stored)
    baseline=(run_agent(inputs[0],inputs[1]),run_agent(inputs[2],inputs[3]))
    if cache is not None:
        temporary=cache.with_suffix('.tmp')
        temporary.write_text(json.dumps(baseline,separators=(',',':'),allow_nan=False))
        temporary.chmod(0o600)
        temporary.replace(cache)
    return baseline

def test_authentic_protein_audit(outputs):
    protein=''.join(line.strip() for line in (ROOT/'assets/NP_001390.1.faa').read_text().splitlines() if not line.startswith('>'))
    assert protein[306:308]=='VE' and len(protein)==391
    for output in outputs:
        audit=output['audit.json']
        assert set(audit)=={'accession','a1_length','a2_length','a1_sha256','a2_sha256','missing_a2_residues','coordinates_are_synthetic','curated_concordance'}
        assert audit['accession']=='NP_001390.1'
        assert audit['a1_length']==391 and audit['a2_length']==389
        assert audit['a1_sha256']==hashlib.sha256(protein.encode()).hexdigest()
        assert audit['a2_sha256']==hashlib.sha256((protein[:306]+protein[308:]).encode()).hexdigest()
        # Both documented JSON encodings must describe the same ordered loss.
        annotation=audit['missing_a2_residues']
        assert (type(annotation) is str and annotation=='Val307-Glu308') or (
            type(annotation) is list and annotation==['Val307', 'Glu308']
        )
        assert audit['coordinates_are_synthetic'] is True

def permissible_evidence_deltas(row):
    """The two documented six-decimal rounding orders for beta evidence.

    Rounding each panel before subtraction can differ by one unit in the
    last decimal from subtracting the full-precision beta estimates first.
    These are alternative deterministic representations, not a score margin.
    """
    a=row['case']['a1_fraction'];b=row['control']['a1_fraction']
    if a is None or b is None:
        return {None}
    shown=round(a-b,6)
    raw=[]
    for label in ('case','control'):
        c=row[label]['counts']
        assert row[label]['usable_unique']==c['A1_unique']+c['A2_unique']
        raw.append((c['A1_unique']+1.5)/(row[label]['usable_unique']+3))
    return {shown,round(raw[0]-raw[1],6)}


def compare_independent_evidence(got,expected):
    assert set(got)=={'records'} and set(expected)=={'records'}
    assert len(got['records'])==len(expected['records'])
    for row,correct in zip(got['records'],expected['records']):
        # Retain exact count/coverage/interval/class verification. Only the
        # last-decimal rounding order for the delta has two valid forms.
        assert {k:v for k,v in row.items() if k!='delta_a1_fraction'} == {
            k:v for k,v in correct.items() if k!='delta_a1_fraction'
        }
        value=row['delta_a1_fraction']
        if value is not None:
            assert_six_decimal_number(value)
            assert -1<=value<=1
        assert value in permissible_evidence_deltas(row), (row['case_id'],value)


def test_sealed_coding_evidence(outputs):
    assert outputs[0]['coding.json']==load_exact(ROOT/'sealed_expected.json')
    compare_independent_evidence(outputs[0]['evidence.json'],load_exact(ROOT/'sealed_evidence_expected.json'))

def test_holdout_coding_evidence(outputs):
    assert outputs[1]['coding.json']==load_exact(ROOT/'hidden_expected.json')
    compare_independent_evidence(outputs[1]['evidence.json'],load_exact(ROOT/'hidden_evidence_expected.json'))
    # Demonstrate that the independent unrounded-first convention passes
    # without making count, interval, classification or schema checks lax.
    modified=copy.deepcopy(outputs[1]['evidence.json'])
    changed=0
    for row in modified['records']:
        possibilities=permissible_evidence_deltas(row)
        value=next((v for v in possibilities if v!=row['delta_a1_fraction']),None)
        if value is not None:
            changed+=1
            row['delta_a1_fraction']=value
    assert changed>0,'Fixtures must exercise intermediate-rounding distinction'
    compare_independent_evidence(modified,load_exact(ROOT/'hidden_evidence_expected.json'))
    corrupted=copy.deepcopy(modified)
    target=next(r for r in corrupted['records'] if r['delta_a1_fraction'] is not None)
    target['delta_a1_fraction']=round(target['delta_a1_fraction']+.000003,6)
    with pytest.raises(AssertionError):
        compare_independent_evidence(corrupted,load_exact(ROOT/'hidden_evidence_expected.json'))

def test_output_schemas_and_order(outputs):
    for output in outputs:
        coding=output['coding.json'];obs=output['evidence.json'];audit=output['audit.json']
        assert set(coding)=={'records','reference_id'}
        assert set(obs)=={'records'} and isinstance(obs['records'],list)
        assert [r['case_id'] for r in obs['records']]==sorted(set(r['case_id'] for r in obs['records']))
        ids=[(r['case_id'],r['transcript_id']) for r in coding['records']]
        assert ids==sorted(ids) and len(set(ids))==len(ids)
        assert len(coding['records'])==len(obs['records'])*2
        for item in coding['records']:
            assert set(item)=={'case_id','transcript_id','status','reason_codes','product'}
            assert type(item['case_id']) is str and type(item['transcript_id']) is str
            assert item['status'] in ('resolved','unresolved')
            assert type(item['reason_codes']) is list
            if item['status']=='unresolved':
                assert item['product'] is None and item['reason_codes']==sorted(set(item['reason_codes']))
                assert item['reason_codes']
            else:
                assert item['reason_codes']==[]
                product=item['product']
                assert type(product) is dict
                assert set(product)=={'consequence','reference_cds_sha256','altered_cds_sha256',
                                     'reference_protein_sha256','altered_protein_sha256',
                                     'reference_cds_length','altered_cds_length',
                                     'reference_protein_length','altered_protein_length',
                                     'first_changed_residue','premature_stop','stop_observed'}
                for field in ('reference_cds_length','altered_cds_length','reference_protein_length',
                              'altered_protein_length'):
                    assert type(product[field]) is int and product[field]>=0, (item['case_id'],field)
                assert product['first_changed_residue'] is None or (type(product['first_changed_residue']) is int and product['first_changed_residue']>0)
                assert type(product['premature_stop']) is bool and type(product['stop_observed']) is bool
                for field in ('reference_cds_sha256','altered_cds_sha256','reference_protein_sha256','altered_protein_sha256'):
                    assert type(product[field]) is str and len(product[field])==64
                    assert all(ch in '0123456789abcdef' for ch in product[field])
        for row in obs['records']:
            assert set(row)=={'case_id','case','control','delta_a1_fraction','evidence_class'}
            assert row['evidence_class'] in {'insufficient_evidence','A1_enriched','A2_enriched','inconclusive'}
            for k in ('case','control'):
                panel=row[k]
                assert set(panel)=={'counts','deduplicated','usable_unique','a1_fraction','interval_95'}
                assert set(panel['counts'])=={'A1_unique','A2_unique','shared'}
                assert all(type(v) is int and v>=0 for v in panel['counts'].values())
                assert panel['usable_unique']==panel['counts']['A1_unique']+panel['counts']['A2_unique']
                assert panel['deduplicated']==sum(panel['counts'].values())
                if panel['usable_unique']<12:
                    assert panel['a1_fraction'] is None and panel['interval_95'] is None
                else:
                    assert type(panel['a1_fraction']) is float
                    assert 0<=panel['a1_fraction']<=1
                    assert type(panel['interval_95']) is list and len(panel['interval_95'])==2
                    assert 0<=panel['interval_95'][0]<=panel['a1_fraction']<=panel['interval_95'][1]<=1

def test_known_negative_strand_consequences(outputs):
    rows={(x['case_id'],x['transcript_id']):x for x in outputs[0]['coding.json']['records']}
    # The 3-nt insertion in m04 shifts the coordinate of an intact stop.
    # That is not termination loss, so both transcripts must be in-frame insertions.
    for tx in ('EDA_SIM_A1.1','EDA_SIM_A2.1'):
        product=rows[('m04',tx)]['product']
        assert product['consequence']=='inframe_insertion'
        assert product['stop_observed'] is True and product['premature_stop'] is False
        assert product['altered_cds_length']==product['reference_cds_length']+3
        assert product['altered_protein_length']==product['reference_protein_length']+1
    # The zero-based insertion start is the cut BEFORE the genomic base.
    # Independently reconstruct m04 on the negative strand to pin the
    # before-base convention (using after-base would yield a different hash).
    ref=load_exact(ROOT/'assets/reference.json')
    genome=ref['genome']
    mutation=next(c for c in load_exact(ROOT/'sealed_cases.json')['cases'] if c['case_id']=='m04')['variants'][0]
    start=mutation['start']
    rc=lambda seq:seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
    for tx in ref['transcripts']:
        exons=tx['exons']
        original=''.join(rc(genome[a:b]) for a,b in exons)
        interior=[(sum(e-s for s,e in exons[:i])+b-start) for i,(a,b) in enumerate(exons) if a<start<b]
        assert len(interior)==1
        pos=interior[0]
        projected=original[:pos]+rc(mutation['alt'])+original[pos:]
        digest=hashlib.sha256(projected.encode('ascii')).hexdigest()
        observed=rows[('m04',tx['transcript_id'])]['product']['altered_cds_sha256']
        assert observed==digest
    assert rows[('m08','EDA_SIM_A1.1')]['product']['consequence']!='no_change'
    assert rows[('m08','EDA_SIM_A2.1')]['product']['consequence']!='no_change'
    assert rows[('m09','EDA_SIM_A2.1')]['status']=='unresolved'
    for name in ['EDA_SIM_A1.1','EDA_SIM_A2.1']:
        assert rows[('m10',name)]['reason_codes']==['overlapping_edits']
        assert rows[('m13',name)]['reason_codes']==['reference_mismatch']
        assert rows[('m17',name)]['product']['consequence']=='frameshift'

def test_molecule_coverage_and_uncertainty(outputs):
    for output in outputs:
        records=output['evidence.json']['records']
        # Every sealed/held-out molecule panel is informative. Underpowered
        # panels are verified separately below as purposeful controls.
        assert all(r['evidence_class']!='insufficient_evidence' for r in records)
        assert all(r[group]['usable_unique']>=12 for r in records for group in ('case','control'))
        assert any(r['case']['counts']['shared']>0 for r in records)
        assert any(r['case']['deduplicated']>r['case']['usable_unique'] for r in records)

def test_case_order_independence(inputs):
    cases=dict(inputs[0]);cases['cases']=list(reversed(cases['cases']))
    assert run_agent(cases,inputs[1])['coding.json']==load_exact(ROOT/'sealed_expected.json')

def test_molecule_order_independence(inputs):
    data={'observations':list(reversed(inputs[1]['observations']))}
    got=run_agent(inputs[0],data)['evidence.json']
    expected=load_exact(ROOT/'sealed_evidence_expected.json')
    compare_independent_evidence(got,expected)
    alternate=copy.deepcopy(got)
    for row in alternate['records']:
        options=permissible_evidence_deltas(row)
        different=options-{row['delta_a1_fraction']}
        if different: row['delta_a1_fraction']=next(iter(different))
    compare_independent_evidence(alternate,expected)

def test_single_case_reconciliation(inputs):
    cid=inputs[0]['cases'][8]['case_id']
    cases={'cases':[c for c in inputs[0]['cases'] if c['case_id']==cid]}
    obs={'observations':[o for o in inputs[1]['observations'] if o['case_id']==cid]}
    actual=run_agent(cases,obs)
    exp_coding=load_exact(ROOT/'sealed_expected.json')
    assert actual['coding.json']['records']==[r for r in exp_coding['records'] if r['case_id']==cid]
    exp_ev=load_exact(ROOT/'sealed_evidence_expected.json')
    compare_independent_evidence(actual['evidence.json'],
                                 {'records':[r for r in exp_ev['records'] if r['case_id']==cid]})


def test_no_unexplained_identity_baseline(outputs,inputs):
    """Every heldout case must depend on its supplied variant evidence."""
    for fixture in (inputs[0], inputs[2], load_exact(ROOT/'assets/practice.json')):
        assert all(set(case)=={'case_id','variants'} for case in fixture['cases']), 'Input must not disclose answer labels'
    actual={(r['case_id'],r['transcript_id']):r for r in outputs[0]['coding.json']['records']}
    for case in inputs[0]['cases']:
        assert case['variants'], 'Every sealed challenge contains at least one phased edit'
        assert set(case)=={'case_id','variants'}, 'No consequence labels may leak in graded inputs'
        # Test the previously reported zero-information cases independently.
        if case['case_id'] in ('m08','m18'):
            for tx in ('EDA_SIM_A1.1','EDA_SIM_A2.1'):
                r=actual[(case['case_id'],tx)]
                assert r['status']=='unresolved' or r['product']['consequence']!='no_change'



def assert_six_decimal_number(value):
    """Assert JSON numeric type, finiteness and at most six significant decimal places.

    The tolerance for comparing *two independently computed* estimates is
    separate from the output precision contract.  E.g., an estimate with a
    nonzero seventh decimal place may be numerically close but is invalid.
    Decimal(str(...)) validates the published JSON numeric value rather than
    the binary floating point storage representation.
    """
    assert type(value) in (int, float), 'Likelihood scalar must be a JSON number, not bool/string'
    assert math.isfinite(value), 'Likelihood scalar must be finite'
    value_decimal = Decimal(str(value))
    parts = value_decimal.as_tuple()
    trailing_zeros = 0
    for digit in reversed(parts.digits):
        if digit != 0:
            break
        trailing_zeros += 1
    fractional_places = max(0, -parts.exponent - trailing_zeros)
    assert fractional_places <= 6, f'Likelihood value exceeds six-decimal precision: {value}'


def assert_nullable_six_decimal(value):
    if value is not None:
        assert_six_decimal_number(value)


def assert_phased_no_call(panel):
    assert set(panel)=={'hap0_a1_fraction','hap1_a1_fraction','allelic_contrast',
                       'hap0_interval_95','hap1_interval_95','informative_molecules',
                       'batches_used','quality_flag'}
    assert panel['quality_flag']=='unidentifiable'
    assert type(panel['informative_molecules']) is int and panel['informative_molecules']>=0
    assert type(panel['batches_used']) is int and 0<=panel['batches_used']<=3
    for field in ('hap0_a1_fraction','hap1_a1_fraction','allelic_contrast',
                  'hap0_interval_95','hap1_interval_95'):
        assert panel[field] is None, f'Incomplete phased no-call: {field}'


def compare_likelihood(observed,expected):
    assert set(observed)=={'calibration_id','records'}
    assert observed['calibration_id']==expected['calibration_id']
    assert type(observed['records']) is list
    assert len(observed['records'])==len(expected['records'])
    assert [r['case_id'] for r in observed['records']]==sorted(r['case_id'] for r in expected['records'])
    for got,want in zip(observed['records'],expected['records']):
        assert set(got)=={'case_id','case','control','delta_a1_mle','mixture_class'}
        assert got['case_id']==want['case_id'] and got['mixture_class']==want['mixture_class']
        assert_nullable_six_decimal(got['delta_a1_mle'])
        if want['delta_a1_mle'] is None:
            assert got['delta_a1_mle'] is None
        else:
            assert got['delta_a1_mle']==pytest.approx(want['delta_a1_mle'],abs=3e-6,rel=0)
        # A displayed delta must be consistent with the already-rounded panel
        # estimates; agreement with a fixture alone cannot replace this rule.
        estimates=(got['case']['a1_mle'],got['control']['a1_mle'])
        if any(v is None for v in estimates):
            assert got['delta_a1_mle'] is None
        else:
            assert -1<=got['delta_a1_mle']<=1
            assert got['delta_a1_mle']==round(estimates[0]-estimates[1],6)
        for group in ('case','control'):
            a,b=got[group],want[group]
            assert set(a)=={'molecules','a1_mle','profile_95','log_likelihood'}
            assert type(a['molecules']) is int and a['molecules']>=0, 'molecules must be a nonnegative integer'
            assert a['molecules']==b['molecules']
            for key in ('a1_mle','log_likelihood'):
                assert_nullable_six_decimal(a[key])
                if b[key] is None:
                    assert a[key] is None
                else:
                    assert a[key]==pytest.approx(b[key],abs=3e-6,rel=0)
                    if key=='a1_mle':
                        assert 0<=a[key]<=1, 'Likelihood estimate must lie in [0,1]'
            if b['profile_95'] is None:
                assert a['profile_95'] is None
            else:
                assert type(a['profile_95']) is list and len(a['profile_95'])==2
                for endpoint in a['profile_95']:
                    assert_six_decimal_number(endpoint)
                assert a['profile_95']==pytest.approx(b['profile_95'],abs=3e-6,rel=0)
                assert 0<=a['profile_95'][0]<=a['a1_mle']<=a['profile_95'][1]<=1

def test_calibrated_mixture_against_independent_expected(outputs):
    compare_likelihood(outputs[0]['likelihood.json'],load_exact(ROOT/'sealed_likelihood_expected.json'))
    compare_likelihood(outputs[1]['likelihood.json'],load_exact(ROOT/'hidden_likelihood_expected.json'))


def test_likelihood_molecule_counts_require_integer_json_type(outputs):
    """Regression: 43.0 == 43 must not allow a floating molecule count."""
    expected=load_exact(ROOT/'sealed_likelihood_expected.json')
    valid=json.loads(json.dumps(outputs[0]['likelihood.json']))
    compare_likelihood(valid,expected)
    for bad_count in (float(valid['records'][0]['case']['molecules']), True,
                      str(valid['records'][0]['case']['molecules'])):
        corrupt=json.loads(json.dumps(valid))
        corrupt['records'][0]['case']['molecules']=bad_count
        with pytest.raises(AssertionError):
            compare_likelihood(corrupt,expected)


def test_likelihood_numeric_output_must_be_six_decimal_rounded(outputs):
    """Regression: values within the numerical tolerance still require valid precision."""
    expected=load_exact(ROOT/'sealed_likelihood_expected.json')
    valid=json.loads(json.dumps(outputs[0]['likelihood.json']))
    compare_likelihood(valid,expected)
    mutations=[
        ('delta_a1_mle',None),
        ('a1_mle','case'),
        ('log_likelihood','control'),
        ('profile_95','case'),
    ]
    for key,group in mutations:
        corrupt=json.loads(json.dumps(valid))
        row=corrupt['records'][0]
        if key=='profile_95':
            value=row[group][key][0]
            assert value is not None
            row[group][key][0]=round(value,6)+1e-7
        elif group is None:
            value=row[key]
            assert value is not None
            row[key]=round(value,6)+1e-7
        else:
            value=row[group][key]
            assert value is not None
            row[group][key]=round(value,6)+1e-7
        with pytest.raises(AssertionError):
            compare_likelihood(corrupt,expected)
    # Keep expected and observed *identically* corrupted: these mutation
    # controls only fail if our standalone schema/identity checks work.
    baseline=json.loads(json.dumps(valid))
    for edit in ('estimate_above_one','negative_interval_endpoint','delta_identity'):
        corrupt=json.loads(json.dumps(baseline))
        row=corrupt['records'][0]
        if edit=='estimate_above_one':
            row['case']['a1_mle']=1.01
            row['case']['profile_95']=[1.01,1.01]
            row['delta_a1_mle']=round(1.01-row['control']['a1_mle'],6)
        elif edit=='negative_interval_endpoint':
            row['case']['profile_95'][0]=-0.01
        else:
            row['delta_a1_mle']=round(row['delta_a1_mle']+0.012345,6)
        with pytest.raises(AssertionError):
            compare_likelihood(corrupt,corrupt)


def test_each_graded_likelihood_instance_requires_estimation(outputs):
    """Every sealed instance distinguishes true inference from a null-fit shortcut.

    Missing-information decisions are independently exercised by the dedicated
    low-coverage control below, not counted as graded evidence-rich instances.
    """
    for output in outputs:
        rows=output['likelihood.json']['records']
        assert rows, 'Expected nonempty cases'
        assert all(r['case']['molecules']>=12 and r['control']['molecules']>=12 for r in rows)
        assert all(r['case']['a1_mle'] is not None and r['control']['a1_mle'] is not None for r in rows)
        assert all(r['case']['profile_95'] is not None and r['control']['profile_95'] is not None for r in rows)
        assert all(r['delta_a1_mle'] is not None for r in rows)
        assert all(r['mixture_class'] != 'insufficient_information' for r in rows)
        # Nonconstant estimates: the inference is evidence-dependent.
        estimates={(r['case']['a1_mle'],r['control']['a1_mle']) for r in rows}
        assert len(estimates)>=6

def test_calibration_change_is_reflected(inputs,outputs):
    calibr=load_exact(ROOT/'assets/spikein_calibration.json')
    changed=json.loads(json.dumps(calibr))
    changed['spike_ins'][0]['counts']={'A1_unique':650,'A2_unique':20,'shared':330}
    changed['spike_ins'][1]['counts']={'A1_unique':18,'A2_unique':620,'shared':362}
    sample_cases={'cases':inputs[0]['cases'][:3]}
    allowed={c['case_id'] for c in sample_cases['cases']}
    observations={'observations':[o for o in inputs[1]['observations'] if o['case_id'] in allowed]}
    output=run_agent(sample_cases,observations,changed)['likelihood.json']
    # The verification expected values are generated with a separate derivative-based method.
    import sys
    sys.path.insert(0,str(ROOT.parent/'authoring/provenance')) if (ROOT.parent/'authoring/provenance').exists() else None
    # In the verifier image this independent audit module is sealed in /tests.
    from independent_calibration import calculate
    evidence=load_exact(ROOT/'sealed_evidence_expected.json')['records']
    expected=calculate([e for e in evidence if e['case_id'] in allowed],changed)
    compare_likelihood(output,expected)
    original=[r for r in outputs[0]['likelihood.json']['records'] if r['case_id'] in allowed]
    assert any(x['case']['a1_mle']!=y['case']['a1_mle'] for x,y in zip(output['records'],original) if x['case']['a1_mle'] is not None)

def test_mixture_low_coverage_is_explicit(inputs):
    """Purposeful null-control, separate from the informative graded cases.

    The case panel has less than 12 retained molecules, while its paired
    control panel has enough evidence to estimate a mixture. Independently
    recompute both panels, so a constant-null solver is not sufficient.
    """
    case=inputs[0]['cases'][0]
    cid=case['case_id']
    sample=[o for o in inputs[1]['observations'] if o['case_id']==cid]
    eligible=lambda r: r['mapping_quality']>=30 and r['base_quality']>=20
    chosen_case=[r for r in sample if r['group']=='case' and eligible(r)]
    chosen_control=[r for r in sample if r['group']=='control' and eligible(r)]
    # Select distinct UMIs to guarantee seven retained case molecules.
    unique_case={r['umi']:r for r in chosen_case}
    selected=list(unique_case.values())[:7]
    assert len(selected)==7
    assert len(set(r['umi'] for r in chosen_control))>=20
    observations={'observations': selected + chosen_control}
    actual=run_agent({'cases':[case]},observations)
    from independent_calibration import calculate
    # The independent molecule calculator uses a separate grouping algorithm.
    import sys
    if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
    from independent_evidence import independently_verify
    independently_computed=independently_verify(observations)
    compare_independent_evidence(actual['evidence.json'],independently_computed)
    expected=calculate(independently_computed['records'],load_exact(ROOT/'assets/spikein_calibration.json'))
    compare_likelihood(actual['likelihood.json'],expected)
    row=actual['likelihood.json']['records'][0]
    assert 0<row['case']['molecules']<12
    assert row['case']['a1_mle'] is None and row['case']['profile_95'] is None
    assert row['control']['molecules']>=12 and row['control']['a1_mle'] is not None
    assert row['delta_a1_mle'] is None and row['mixture_class']=='insufficient_information'


def test_variants_order_independence(inputs):
    """Phased alleles may be presented in a different record order."""
    swapped=json.loads(json.dumps(inputs[0]))
    for row in swapped['cases']:
        row['variants']=list(reversed(row['variants']))
    result=run_agent(swapped,inputs[1])
    assert result['coding.json']==load_exact(ROOT/'sealed_expected.json')

def test_low_quality_duplicate_molecules_do_not_change_results(inputs,outputs):
    """Additional low-quality duplicate records must not alter molecule counts."""
    source=inputs[1]['observations']
    added=[]
    for row in source[::73]:
        r=dict(row)
        r['mapping_quality']=min(row['mapping_quality'],5)
        r['base_quality']=min(row['base_quality'],5)
        added.append(r)
    mutated={'observations':added+list(reversed(source))}
    actual=run_agent(inputs[0],mutated)
    compare_independent_evidence(actual['evidence.json'],
                                 load_exact(ROOT/'sealed_evidence_expected.json'))
    assert actual['likelihood.json']==outputs[0]['likelihood.json']

def test_ignoring_all_variants_fails_every_coding_record(inputs,outputs):
    """A zero-edit baseline cannot match any sealed coding outcome."""
    for cases,observations,oracle in [(inputs[0],inputs[1],outputs[0]),(inputs[2],inputs[3],outputs[1])]:
        scrubbed=json.loads(json.dumps(cases))
        for case in scrubbed['cases']:
            case['variants']=[]
        baseline=run_agent(scrubbed,observations)['coding.json']['records']
        truth=oracle['coding.json']['records']
        assert len(baseline)==len(truth)
        assert all(b!=t for b,t in zip(baseline,truth))


def test_null_and_boundary_cases_are_scientifically_checked(inputs):
    """Unchanged, missing-field, partial-boundary and multi-exon cases are tested.

    Expected reasons are defined from geometry and input validity, not from a
    private invocation of the reference solver.
    """
    ref=load_exact(ROOT/'assets/reference.json')
    genome=ref['genome']
    alt=lambda b: 'A' if b != 'A' else 'C'
    probes=[
        {'case_id':'qa_empty','variants':[]},
        {'case_id':'qa_intronic','variants':[{'start':9500,'ref':genome[9500],'alt':alt(genome[9500])}]},
        {'case_id':'qa_missing','variants':[{'start':10140,'ref':genome[10140]}]},
        {'case_id':'qa_out_of_range','variants':[{'start':len(genome)+1,'ref':'A','alt':'C'}]},
        {'case_id':'qa_boundary','variants':[{'start':10180,'ref':'','alt':'A'}]},
        {'case_id':'qa_partial','variants':[{'start':10178,'ref':genome[10178:10183],'alt':''}]},
        {'case_id':'qa_multiexon','variants':[{'start':8290,'ref':genome[8290:10020],'alt':''}]},
    ]
    # Duplicate a valid, simple observation for each probe to exercise the full
    # end-to-end application while keeping the separate RNA stage deterministic.
    example=[row for row in inputs[1]['observations'] if row['case_id']==inputs[0]['cases'][0]['case_id']]
    assert example
    observations={'observations':[{**row,'case_id':case['case_id']} for case in probes for row in example]}
    output=run_agent({'cases':probes},observations)
    expected={
      'qa_empty':('resolved',set()),'qa_intronic':('resolved',set()),
      'qa_missing':('unresolved',{'invalid_edit'}),
      'qa_out_of_range':('unresolved',{'invalid_edit'}),
      'qa_boundary':('unresolved',{'exon_boundary_ambiguity'}),
      'qa_partial':('unresolved',{'partial_exon_overlap'}),
      'qa_multiexon':('unresolved',{'multiple_exon_overlap','partial_exon_overlap'}),
    }
    rows=output['coding.json']['records']
    assert len(rows)==len(probes)*2
    for row in rows:
        status,reasons=expected[row['case_id']]
        assert row['status']==status
        assert set(row['reason_codes'])==reasons
        if status=='resolved':
            assert row['product']['consequence']=='no_change'
            assert row['product']['reference_cds_sha256']==row['product']['altered_cds_sha256']
        else:
            assert row['product'] is None


def test_deliverables_are_standard_library_only():
    """Agent deliverables may import only standard-library or supplied modules.

    The subprocess uses `-I -S` to exclude runtime site packages. The source
    inspection catches direct forbidden imports even when a branch is not used.
    """
    allowed = set(sys.stdlib_module_names) | {Path(name).stem for name in CODE_NAMES}
    for name in CODE_NAMES:
        path=SOURCE/name
        tree=ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):
                modules=[item.name.split('.')[0] for item in node.names]
            elif isinstance(node,ast.ImportFrom):
                assert node.level==0, f'Relative import not supported: {name}'
                modules=[(node.module or '').split('.')[0]]
            else:
                continue
            for module in modules:
                assert module in allowed, f'Non-standard import in {name}: {module}'


def test_full_order_invariance_all_outputs(inputs,outputs):
    """Coding, RNA evidence and likelihood must all survive input reorderings."""
    cases={'cases':list(reversed(inputs[0]['cases']))}
    observations={'observations':list(reversed(inputs[1]['observations']))}
    for case in cases['cases']:
        case['variants']=list(reversed(case['variants']))
    got=run_agent(cases,observations)
    for filename in ('coding.json','likelihood.json','audit.json'):
        assert got[filename]==outputs[0][filename],filename
    compare_independent_evidence(got['evidence.json'],
                                 load_exact(ROOT/'sealed_evidence_expected.json'))


def test_umi_ranking_with_conflicting_signal_and_qualities(inputs):
    """The highest mapping-quality / base-quality / signal rank must win per UMI."""
    from independent_evidence import independently_verify
    from independent_calibration import calculate
    chosen=inputs[0]['cases'][0]
    cid=chosen['case_id']
    rows=[]
    for group in ('case','control'):
        for i in range(21):
            main_signal=('A1_unique' if i<4 else 'A2_unique' if i<8 else ('A1_unique' if ((i + (group=='case')) % 2 == 0) else 'A2_unique'))
            base={'case_id':cid,'group':group,'umi':f'{group}_quality_{i:03d}',
                  'signal':main_signal,'mapping_quality':50,'base_quality':40}
            rows.append(base)
            if i<4:
                rows.append({**base,'signal':'shared','mapping_quality':35,'base_quality':95})
            elif i<8:
                rows.append({**base,'signal':'shared','mapping_quality':60,'base_quality':23})
            elif i<12:
                rows.append({**base,'signal':'shared'})
            if i%7==0:
                rows.append({**base,'signal':'shared','mapping_quality':5,'base_quality':99})
    sample={'observations':rows}
    got=run_agent({'cases':[chosen]},sample)
    ev=independently_verify(sample)
    compare_independent_evidence(got['evidence.json'],ev)
    assert all(r[group]['deduplicated']==21 and r[group]['usable_unique']==13
               for r in ev['records'] for group in ('case','control'))
    expected=calculate(ev['records'],load_exact(ROOT/'assets/spikein_calibration.json'))
    compare_likelihood(got['likelihood.json'],expected)
    reversed_got=run_agent({'cases':[chosen]},{'observations':list(reversed(rows))})
    compare_independent_evidence(reversed_got['evidence.json'],ev)
    assert reversed_got['likelihood.json']==got['likelihood.json']


def test_shared_molecules_are_not_unique_isoform_support(inputs):
    """19 retained UMIs but only 9 isoform-specific: evidence fractions must be null.

    The calibrated mixture does use shared junction molecules, so its likelihood
    is independently tested rather than incorrectly forced to null.
    """
    from independent_evidence import independently_verify
    from independent_calibration import calculate
    chosen=inputs[0]['cases'][0]
    cid=chosen['case_id']
    rows=[]
    for group in ('case','control'):
        for i in range(19):
            sig='shared' if i>=9 else ('A1_unique' if i%2==0 else 'A2_unique')
            rows.append({'case_id':cid,'group':group,'umi':f'{group}_test_{i:03d}',
                         'signal':sig,'mapping_quality':55,'base_quality':30})
    observations={'observations':rows}
    got=run_agent({'cases':[chosen]},observations)
    ev=independently_verify(observations)
    compare_independent_evidence(got['evidence.json'],ev)
    for group in ('case','control'):
        panel=got['evidence.json']['records'][0][group]
        assert panel['deduplicated']==19 and panel['usable_unique']==9
        assert panel['a1_fraction'] is None and panel['interval_95'] is None
    expected=calculate(ev['records'],load_exact(ROOT/'assets/spikein_calibration.json'))
    compare_likelihood(got['likelihood.json'],expected)
    assert got['likelihood.json']['records'][0]['case']['a1_mle'] is not None


def test_calibration_changes_pseudocount_coverage_and_mle_endpoints(inputs):
    """Altered calibration, profile-likelihood boundary optima and coverage gate."""
    from independent_evidence import independently_verify
    from independent_calibration import calculate
    chosen=inputs[0]['cases'][0]
    cid=chosen['case_id']
    rows=[]
    for group in ('case','control'):
        sig='A1_unique' if group=='case' else 'A2_unique'
        for i in range(55):
            rows.append({'case_id':cid,'group':group,'umi':f'{group}_endpoint_{i:03d}',
                         'signal':sig,'mapping_quality':61,'base_quality':39})
        for i in range(6):
            rows.append({'case_id':cid,'group':group,'umi':f'{group}_ambig_{i:03d}',
                         'signal':'shared','mapping_quality':61,'base_quality':39})
    observations={'observations':rows}
    ev=independently_verify(observations)
    calibration=load_exact(ROOT/'assets/spikein_calibration.json')
    changed=json.loads(json.dumps(calibration))
    changed['calibration_id']='EDA_TEST_ALTERNATE'
    changed['pseudocount']=2.0
    changed['minimum_molecules']=40
    changed['spike_ins'][0]['counts']={'A1_unique':700,'A2_unique':20,'shared':280}
    changed['spike_ins'][1]['counts']={'A1_unique':25,'A2_unique':700,'shared':275}
    got=run_agent({'cases':[chosen]},observations,changed)
    compare_likelihood(got['likelihood.json'],calculate(ev['records'],changed))
    record=got['likelihood.json']['records'][0]
    assert record['case']['a1_mle']==1.0 and record['control']['a1_mle']==0.0
    assert record['case']['profile_95'][1]==1.0 and record['control']['profile_95'][0]==0.0
    # A higher minimum causes the same 61 unique molecules to be insufficient.
    higher=json.loads(json.dumps(changed))
    higher['minimum_molecules']=70
    output=run_agent({'cases':[chosen]},observations,higher)
    compare_likelihood(output['likelihood.json'],calculate(ev['records'],higher))
    for group in ('case','control'):
        assert output['likelihood.json']['records'][0][group]['a1_mle'] is None


def test_invalid_allele_pair_and_coordinate_types(inputs):
    """Check prohibited both-empty, boolean/float positions, and wrong bases."""
    reference=load_exact(ROOT/'assets/reference.json')
    nucleotide=reference['genome'][10140]
    invalid=[
      {'case_id':'qc_empty_both','variants':[{'start':10140,'ref':'','alt':''}]},
      {'case_id':'qc_float_start','variants':[{'start':10140.5,'ref':nucleotide,'alt':'A'}]},
      {'case_id':'qc_boolean_start','variants':[{'start':True,'ref':nucleotide,'alt':'A'}]},
      {'case_id':'qc_invalid_base','variants':[{'start':10140,'ref':nucleotide,'alt':'N'}]},
      {'case_id':'qc_paired_insertions','variants':[{'start':10140,'ref':'','alt':'A'},
                                                    {'start':10140,'ref':'','alt':'C'}]},
      {'case_id':'qc_insertion_snv_same_start','variants':[{'start':10140,'ref':'','alt':'G'},
                                               {'start':10140,'ref':nucleotide,'alt':next(b for b in 'ACGT' if b!=nucleotide)}]},
    ]
    original_case=inputs[0]['cases'][0]['case_id']
    sample=[r for r in inputs[1]['observations'] if r['case_id']==original_case]
    observations={'observations':[{**r,'case_id':c['case_id']} for c in invalid for r in sample]}
    got=run_agent({'cases':invalid},observations)['coding.json']['records']
    assert len(got)==len(invalid)*2
    for row in got:
        assert row['status']=='unresolved'
        assert row['product'] is None
        expected_code = ('overlapping_edits' if row['case_id'] in
                         ('qc_paired_insertions','qc_insertion_snv_same_start')
                         else 'invalid_edit')
        assert row['reason_codes']==[expected_code]


def test_real_eda_protein_variant_concordance(outputs):
    # These protein changes and classifications are sourced from public EDA
    # curation records, not from the synthetic genome and RNA simulator.
    result=outputs[0]['audit.json']['curated_concordance']
    assert result['clinical_interpretation_inferred'] is False
    assert result['panel_id']=='EDA_PUBLIC_PROTEIN_HGVS_SNAPSHOT_2026_11_06'
    rows={r['variant_id']:r for r in result['variants']}
    assert len(rows)==6 and len(result['variants'])==len(rows)
    protein=''.join(x.strip() for x in (ROOT/'assets/NP_001390.1.faa').read_text().splitlines() if not x.startswith('>'))
    a2=protein[:306]+protein[308:]
    for ident,position,ref,alt in [
        ('EDA-NM_001399.5-c457C_T',153,'R','C'),
        ('EDA-NM_001399.5-c463C_T',155,'R','C'),
        ('EDA-A1-pArg156His',156,'R','H'),
        ('EDA-NM_001399.5-c1067C_A',356,'A','D'),
        ('EDA-NM_001399.5-c1093G_A',365,'V','M'),
        ('EDA-A1-pArg153His',153,'R','H'),
    ]:
        row=rows[ident]
        for iso,pep,index in [('A1',protein,position),('A2',a2,position if position<=306 else position-2)]:
            data=row['isoforms'][iso]
            assert data['status']=='mapped' and data['position']==index
            assert data['reference_residue']==ref and data['alternate_residue']==alt
            assert pep[index-1]==ref
            mutated=pep[:index-1]+alt+pep[index:]
            assert data['altered_protein_sha256']==hashlib.sha256(mutated.encode()).hexdigest()
    assert rows['EDA-A1-pArg153His']['assertion_concordance']=='discordant'
    assert rows['EDA-A1-pArg153His']['reported_assertions']==['Likely benign','Pathogenic']


def test_curated_protein_panel_reorder_and_isoform_deletion(inputs):
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    panel['variants']=list(reversed(panel['variants']))
    output=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']
    assert output==run_agent(inputs[0],inputs[1])['audit.json']['curated_concordance']
    panel['variants'].append({
        'variant_id':'artificial-Val307-mapping-control','protein_change':'p.Val307Ala',
        'classification_observations':[{'classification':'Unknown','record':'SYNTHETIC-CONTROL'}]
    })
    rows=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']['variants']
    v=next(r for r in rows if r['variant_id']=='artificial-Val307-mapping-control')
    assert v['isoforms']['A1']['status']=='mapped'
    absent=v['isoforms']['A2']
    assert absent['status']=='isoform_absent'
    assert all(absent[field] is None for field in
               ('position','reference_residue','alternate_residue','altered_protein_sha256'))


def test_curated_variant_mismatch_and_assertion_conflicts(inputs):
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    for entry in panel['variants']:
        if entry['variant_id']=='EDA-A1-pArg153His':
            entry['classification_observations']=entry['classification_observations'][:1]
    panel['variants'].append({'variant_id':'artificial-reference-mismatch-control',
        'protein_change':'p.Ala365Asp','classification_observations':[
            {'classification':'Uncertain','record':'SYNTHETIC-CONTROL'}]})
    rows=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']['variants']
    lookup={r['variant_id']:r for r in rows}
    assert lookup['EDA-A1-pArg153His']['assertion_concordance']=='single_class'
    assert lookup['EDA-A1-pArg153His']['reported_assertions']==['Pathogenic']
    mismatch=lookup['artificial-reference-mismatch-control']
    assert all(mismatch['isoforms'][x]['status']=='reference_mismatch' for x in ('A1','A2'))
    assert all(mismatch['isoforms'][x]['altered_protein_sha256'] is None for x in ('A1','A2'))


def test_authentic_human_eda_grch38_crosswalk(outputs):
    # Independently documented GRCh38 EDA locus and ClinVar RCV006905138.1.
    cross=outputs[0]['audit.json']['curated_concordance']['genomic_crosswalk']
    assert cross['assembly']=='GRCh38.p14'
    assert cross['chromosome_accession']=='NC_000023.11'
    assert cross['gene_orientation']=='+'
    assert cross['simulated_workbench_orientation']=='-'
    assert cross['actual_chromosome_sequence_supplied'] is False
    assert cross['gene_interval_1based_closed']==[69616113,70039472]
    assert len(cross['records'])==1
    v=cross['records'][0]
    assert v['genomic_hgvs']=='NC_000023.11:g.69957087C>T'
    assert v['in_annotated_gene'] is True
    assert v['transcript_hgvs']=={'A1':'NM_001399.5:c.457C>T','A2':'NM_001005609.2:c.457C>T'}
    assert v['isoform_status']=={'A1':'mapped','A2':'mapped'}
    assert v['clinical_record']=='RCV006905138.1'
    assert v['genomic_allele_sequence_validated'] is False

def test_authentic_crosswalk_is_evidence_dependent(inputs):
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    panel['genomic_crosswalk_records'][0]['genomic_hgvs']='NC_000023.11:g.70039473C>T'
    result=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']['genomic_crosswalk']
    assert result['records'][0]['in_annotated_gene'] is False
    # The other half of the same explicit interval-membership rule: lower bound.
    panel['genomic_crosswalk_records'][0]['genomic_hgvs']='NC_000023.11:g.69616112C>T'
    result=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']['genomic_crosswalk']
    assert result['records'][0]['in_annotated_gene'] is False
    # Endpoints belong to the gene (1-based CLOSED bounds).
    for pos in (69616113, 70039472):
        panel['genomic_crosswalk_records'][0]['genomic_hgvs']=f'NC_000023.11:g.{pos}C>T'
        result=run_agent(inputs[0],inputs[1],curated_override=panel)['audit.json']['curated_concordance']['genomic_crosswalk']
        assert result['records'][0]['in_annotated_gene'] is True
    # An accession version mismatch must raise nonzero, not be treated as
    # merely outside the gene. Check version and chromosome identity alike.
    for bad in ('NC_000023.10:g.69957087C>T',
                'NC_000024.11:g.69957087C>T',
                'NC_000023:g.69957087C>T'):
        panel['genomic_crosswalk_records'][0]['genomic_hgvs']=bad
        with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
            run_agent(inputs[0],inputs[1],curated_override=panel)

def test_variant_rna_case_id_mismatches_are_rejected(inputs):
    """Both missing RNA cases and additional RNA-only cases violate the input contract."""
    cases={'cases':inputs[0]['cases'][:2]}
    chosen={c['case_id'] for c in cases['cases']}
    observations=[o for o in inputs[1]['observations'] if o['case_id'] in chosen]
    first=cases['cases'][0]['case_id']
    missing={'observations':[o for o in observations if o['case_id']==first]}
    assert set(chosen)-{first}
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(cases,missing)
    extra=json.loads(json.dumps(observations))
    extra.append({**extra[0], 'case_id':'qc_RNA_without_variant'})
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(cases,{'observations':extra})
    duplicated=copy.deepcopy(cases)
    duplicated['cases'].append(copy.deepcopy(duplicated['cases'][0]))
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(duplicated,{'observations':observations})


def test_curated_accession_assertion_schema_and_provenance(outputs):
    """Verify all contract keys, ordered unique source arrays and input-derived records."""
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    expected={r['variant_id']:r for r in panel['variants']}
    for output in outputs:
        c=output['audit.json']['curated_concordance']
        assert set(c)=={'panel_id','source_type','clinical_interpretation_inferred','variants','genomic_crosswalk'}
        assert c['panel_id']==panel['panel_id'] and c['source_type']==panel['source_type']
        assert c['clinical_interpretation_inferred'] is False
        assert type(c['variants']) is list
        assert [v['variant_id'] for v in c['variants']]==sorted(expected)
        for v in c['variants']:
            assert set(v)=={'variant_id','protein_change','source_accessions','reported_assertions','assertion_concordance','isoforms'}
            assert type(v['variant_id']) is str and type(v['protein_change']) is str
            inp=expected[v['variant_id']]
            assert v['protein_change']==inp['protein_change']
            records=inp['classification_observations']
            accessions=sorted({x['record'] for x in records})
            assertions=sorted({x['classification'] for x in records})
            assert accessions and assertions
            assert v['source_accessions']==accessions
            assert v['reported_assertions']==assertions
            assert all(type(x) is str and x for x in v['source_accessions']+v['reported_assertions'])
            assert v['assertion_concordance']==('discordant' if len(assertions)>1 else 'single_class')
            assert set(v['isoforms'])=={'A1','A2'}
            for iso,record in v['isoforms'].items():
                assert set(record)=={'status','position','reference_residue','alternate_residue','altered_protein_sha256'}
                assert record['status'] in {'mapped','reference_mismatch','isoform_absent'}
                if record['position'] is not None:
                    assert type(record['position']) is int and record['position']>0
                for key in ('reference_residue','alternate_residue'):
                    assert record[key] is None or (type(record[key]) is str and len(record[key])==1)
                if record['status']=='mapped':
                    digest=record['altered_protein_sha256']
                    assert type(digest) is str and len(digest)==64
                    assert all(ch in '0123456789abcdef' for ch in digest)
                    assert type(record['position']) is int
                else:
                    assert record['altered_protein_sha256'] is None
                    if record['status']=='isoform_absent':
                        assert record['position'] is None
                        assert record['reference_residue'] is None
                        assert record['alternate_residue'] is None


def test_curated_crosswalk_exact_schema_uniqueness_and_deduplication(inputs,outputs):
    """Crosswalk schema and duplicate-source normalisation follow the brief."""
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    context=panel['genome_context']
    source_rows={r['variant_id']:r for r in panel['genomic_crosswalk_records']}
    for output in outputs:
        cross=output['audit.json']['curated_concordance']['genomic_crosswalk']
        assert set(cross)=={'assembly','chromosome_accession','gene_orientation','gene_interval_1based_closed','simulated_workbench_orientation','actual_chromosome_sequence_supplied','records'}
        assert cross['assembly']==context['assembly']
        assert cross['chromosome_accession']==context['chromosome_accession']
        assert cross['gene_orientation']==context['gene_orientation']
        assert cross['gene_interval_1based_closed']==context['gene_interval_1based_closed']
        assert cross['simulated_workbench_orientation']=='-'
        assert cross['actual_chromosome_sequence_supplied'] is False
        assert type(cross['records']) is list
        ids=[r['variant_id'] for r in cross['records']]
        assert ids==sorted(source_rows) and len(ids)==len(set(ids))
        variants={v['variant_id']:v for v in output['audit.json']['curated_concordance']['variants']}
        for row in cross['records']:
            assert set(row)=={'variant_id','genomic_hgvs','in_annotated_gene','transcript_hgvs','isoform_status','clinical_record','genomic_allele_sequence_validated'}
            source=source_rows[row['variant_id']]
            assert row['genomic_hgvs']==source['genomic_hgvs']
            assert type(row['in_annotated_gene']) is bool
            assert set(row['transcript_hgvs'])=={'A1','A2'}
            assert row['transcript_hgvs']==source['transcript_hgvs']
            assert set(row['isoform_status'])=={'A1','A2'}
            assert row['isoform_status']=={i:variants[row['variant_id']]['isoforms'][i]['status'] for i in ('A1','A2')}
            assert type(row['clinical_record']) is type(source['clinical_record'])
            assert row['clinical_record']==source['clinical_record']
            assert row['genomic_allele_sequence_validated'] is False
    # Repeated source observations must not create duplicate output accessions.
    changed=json.loads(json.dumps(panel))
    changed['variants'][0]['classification_observations'].append(dict(changed['variants'][0]['classification_observations'][0]))
    got=run_agent({'cases':inputs[0]['cases'][:1]},
                  {'observations':[r for r in inputs[1]['observations'] if r['case_id']==inputs[0]['cases'][0]['case_id']]},
                  curated_override=changed)['audit.json']['curated_concordance']
    baseline=outputs[0]['audit.json']['curated_concordance']
    assert got==baseline
    # A duplicated crosswalk identifier must be rejected, not silently retained.
    changed['genomic_crosswalk_records'].append(dict(changed['genomic_crosswalk_records'][0]))
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent({'cases':inputs[0]['cases'][:1]},
                  {'observations':[r for r in inputs[1]['observations'] if r['case_id']==inputs[0]['cases'][0]['case_id']]},
                  curated_override=changed)


def test_reference_protein_anchors_are_mutation_sensitive(inputs):
    """A changed FASTA residue or synthetic CDS must invalidate protein concordance."""
    cases={'cases':[inputs[0]['cases'][0]]}
    cid=cases['cases'][0]['case_id']
    obs={'observations':[row for row in inputs[1]['observations'] if row['case_id']==cid]}
    fasta=(ROOT/'assets/NP_001390.1.faa').read_text()
    lines=fasta.splitlines()
    protein=''.join(lines[1:])
    assert len(protein)==391
    replacement='A' if protein[152]!='A' else 'V'
    edited=protein[:152]+replacement+protein[153:]
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(cases,obs,fasta_override=lines[0]+'\n'+edited+'\n')
    reference=load_exact(ROOT/'assets/reference.json')
    tx=reference['transcripts'][0]
    # On the negative strand, the last genomic base of the first exon
    # is the first CDS base. Altering it disrupts the initiation codon.
    genome=reference['genome']
    first_exon=reference['transcripts'][0]['exons'][0]
    position=first_exon[1]-1
    replacement=next(base for base in 'ACGT' if base!=genome[position])
    changed=json.loads(json.dumps(reference))
    changed['genome']=genome[:position]+replacement+genome[position+1:]
    assert changed is not None
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(cases,obs,reference_override=changed)


def test_transcript_hgvs_coordinate_and_allele_reconciliation(inputs):
    """Alleles and codon positions must agree with the reported protein variant."""
    cases={'cases':[inputs[0]['cases'][0]]}
    cid=cases['cases'][0]['case_id']
    obs={'observations':[row for row in inputs[1]['observations'] if row['case_id']==cid]}
    panel=load_exact(ROOT/'assets/public_eda_variants.json')
    for replacement in ('NM_001399.5:c.458C>T','NM_001399.5:c.457G>A',
                        'NM_001399.5:c.460C>T',
                        'NM_001399.4:c.457C>T'):
        changed=json.loads(json.dumps(panel))
        changed['genomic_crosswalk_records'][0]['transcript_hgvs']['A1']=replacement
        with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
            run_agent(cases,obs,curated_override=changed)
    for replacement in ('NM_001005609.2:c.458C>T',
                        'NM_001005609.2:c.457G>A',
                        'NM_001005609.2:c.460C>T',
                        'NM_001005609.1:c.457C>T'):
        changed=json.loads(json.dumps(panel))
        changed['genomic_crosswalk_records'][0]['transcript_hgvs']['A2']=replacement
        with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
            run_agent(cases,obs,curated_override=changed)
    changed=json.loads(json.dumps(panel))
    changed['genomic_crosswalk_records'][0]['genomic_hgvs']='NC_000023.11:g.69957087C>A'
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(cases,obs,curated_override=changed)



def test_projection_mismatch_requires_independently_consistent_cds_anchor(inputs,outputs):
    """A reference-provided transcript CDS anchor must be reconciled with the
    genomic projection. This direct challenge detects an omitted mismatch
    reason without requiring import of untrusted agent modules.
    """
    reference=load_exact(ROOT/'assets/reference.json')
    case=next(c for c in inputs[0]['cases'] if c['case_id']=='m01')
    # Also include a genuinely intronic-only edit. A discordant CDS anchor
    # must not invalidate variants with no exonic projection.
    intron_pos=2000
    genomic_base=reference['genome'][intron_pos]
    intronic={'case_id':'qc_intron_independent_cds', 'variants':[
        {'start':intron_pos,'ref':genomic_base,'alt':'C' if genomic_base!='C' else 'A'}]}
    assert all(not (a <= intron_pos < b) for tx in reference['transcripts'] for a,b in tx['exons'])
    cases={'cases':[case,intronic]}
    obs={'observations':[
        r for r in inputs[1]['observations'] if r['case_id']=='m01'] + [
        {**r,'case_id':intronic['case_id']} for r in inputs[1]['observations'] if r['case_id']=='m01']}
    strand=reference['strand']
    assert strand=='-'
    complements=str.maketrans('ACGT','TGCA')
    for tx in reference['transcripts']:
        tx['cds_sequence']=''.join(reference['genome'][a:b].translate(complements)[::-1]
                                  for a,b in tx['exons'])
    baseline=run_agent(cases,obs,reference_override=reference)['coding.json']
    expected=[r for r in outputs[0]['coding.json']['records'] if r['case_id']=='m01']
    assert [r for r in baseline['records'] if r['case_id']=='m01']==expected, 'Matching CDS anchors must preserve valid projection'

    # Introduce an independent transcript annotation disagreement without
    # changing the genomic variant, genomic reference or protein anchor.
    changed=copy.deepcopy(reference)
    cds=changed['transcripts'][0]['cds_sequence']
    changed['transcripts'][0]['cds_sequence']=cds[:20]+('T' if cds[20]!='T' else 'A')+cds[21:]
    actual=run_agent(cases,obs,reference_override=changed)['coding.json']
    records={r['transcript_id']:r for r in actual['records'] if r['case_id']=='m01'}
    flagged=records[changed['transcripts'][0]['transcript_id']]
    assert flagged['status']=='unresolved' and flagged['product'] is None
    assert flagged['reason_codes']==['projection_mismatch']
    other=changed['transcripts'][1]['transcript_id']
    assert records[other]==next(r for r in expected if r['transcript_id']==other)
    intronic_results=[r for r in actual['records'] if r['case_id']==intronic['case_id']]
    assert len(intronic_results)==2
    assert all(r['status']=='resolved' and r['reason_codes']==[]
               and r['product']['consequence']=='no_change' for r in intronic_results)


def test_stop_lost_even_when_later_inframe_stop_remains(inputs):
    """Remove the native final CDS stop after extending the synthetic 3' UTR
    with a second in-frame stop; the later stop must not mask stop_lost.
    The pre-edit reference still encodes the authentic EDA-A1/A2 peptides.
    """
    reference=load_exact(ROOT/'assets/reference.json')
    assert reference['strand']=='-'
    terminal_starts={tx['exons'][-1][0] for tx in reference['transcripts']}
    assert len(terminal_starts)==1
    start=terminal_starts.pop()
    genome=reference['genome']
    complements=str.maketrans('ACGT','TGCA')
    revcomp=lambda s:s.translate(complements)[::-1]
    terminal_ref=genome[start:start+3]
    assert revcomp(terminal_ref)=='TAA'
    # On the negative strand, lower genomic coordinates extend the CDS 3' end.
    # The added two codons are alanine then TGA (a downstream stop).
    extension=revcomp('GCTTGA')
    reference['genome']=genome[:start-6]+extension+genome[start:]
    for tx in reference['transcripts']:
        tx['exons'][-1][0]-=6
        tx['cds_end']+=6
    assert all(tx['cds_end']==sum(b-a for a,b in tx['exons']) for tx in reference['transcripts'])
    cid='qc_downstream_stop'
    # The second challenge deletes native TAA and the next GCT; the downstream
    # TGA lands at exactly the native stop's old codon index. Stop INDEX alone
    # must not be used as a proxy for preservation of the NATIVE termination.
    moved='qc_native_stop_deleted_with_stop_relocation'
    moved_start=start-3
    inserted='qc_intact_native_stop_preceded_by_new_stop'
    moved_ref=reference['genome'][moved_start:moved_start+6]
    assert revcomp(moved_ref)=='TAAGCT'
    cases={'cases':[
        {'case_id':cid,'variants':[{'start':start,'ref':terminal_ref,'alt':revcomp('CAA')}]},
        {'case_id':moved,'variants':[{'start':moved_start,'ref':moved_ref,'alt':''}]},
        {'case_id':inserted,'variants':[{'start':start+3,'ref':'','alt':revcomp('TAA')}]},
    ]}
    obs={'observations':[{'case_id':case_id,'group':group,'umi':f'{group}_u1',
                          'signal':'shared','mapping_quality':42,'base_quality':38}
                         for case_id in (cid,moved,inserted) for group in ('case','control')]}
    records=run_agent(cases,obs,reference_override=reference)['coding.json']['records']
    assert len(records)==6
    for r in records:
        assert r['status']=='resolved' and r['reason_codes']==[]
        product=r['product']
        assert product['consequence']==('inframe_insertion' if r['case_id']==inserted else 'stop_lost')
        assert product['stop_observed'] is True, 'The downstream stop must still be detected'
        assert product['premature_stop'] is (r['case_id']==inserted)
        if r['case_id']==inserted:
            assert product['altered_cds_length']==product['reference_cds_length']+3
            assert product['altered_protein_length']==product['reference_protein_length']
            assert product['altered_protein_sha256']==product['reference_protein_sha256']
        elif r['case_id']==cid:
            assert product['altered_protein_length']==product['reference_protein_length']+2
            assert product['first_changed_residue']==product['reference_protein_length']+1
            assert product['altered_cds_length']==product['reference_cds_length']
            assert product['altered_protein_sha256']!=product['reference_protein_sha256']
        else:
            assert product['altered_protein_length']==product['reference_protein_length']
            assert product['altered_cds_length']==product['reference_cds_length']-6
            assert product['altered_protein_sha256']==product['reference_protein_sha256']


def test_multi_exon_reason_codes_are_exhaustive(inputs):
    """A reference-spanning edit must report every applicable ambiguity reason."""
    ref=load_exact(ROOT/'assets/reference.json')
    exons=ref['transcripts'][0]['exons']
    # Find exons close in genomic coordinate even though transcript order is reversed.
    exons=sorted(exons)
    pair=next(((a,b) for a,b in zip(exons,exons[1:]) if a[1]<b[0]),None)
    assert pair is not None
    start=pair[0][1]-1;end=pair[1][0]+1
    cases={'cases':[{'case_id':'qc_multiexon','variants':[{'start':start,'ref':ref['genome'][start:end],'alt':'A'}]}]}
    obs={'observations':[]}
    # A minimal matching observation set satisfies the shared-input contract.
    # No molecules means the inference will be explicitly insufficient.
    for group in ('case','control'):
        obs['observations'].append({'case_id':'qc_multiexon','group':group,'umi':'u_'+group,
                                   'signal':'shared','mapping_quality':40,'base_quality':35})
    records=run_agent(cases,obs)['coding.json']['records']
    main=next(r for r in records if r['transcript_id']==ref['transcripts'][0]['transcript_id'])
    assert main['status']=='unresolved' and main['product'] is None
    assert {'multiple_exon_overlap','partial_exon_overlap'} <= set(main['reason_codes'])
    assert main['reason_codes']==sorted(set(main['reason_codes']))


def test_replicated_batch_assay_ground_truth(outputs):
    for output,tag in zip(outputs,('sealed','hidden')):
        got=output['batch_inference.json']
        truth={x['case_id']:x for x in load_exact(ROOT/(tag+'_batch_truth.json'))['records']}
        assert set(got)=={'assay_id','records'}
        assert got['assay_id']=='EDA-batch-assay-v1'
        assert [r['case_id'] for r in got['records']]==sorted(truth)
        # Designed positive/negative controls must discriminate for every case;
        # a method-free half-fraction with [0.35,0.65] cannot satisfy even one
        # complete case record (both conditions are graded together).
        assert all(max(abs(t[g+'_a1_fraction']-.5) for g in ('case','control'))>=.18
                   for t in truth.values())
        # Per-instance null check: naive raw pooling with a fixed +/-0.11
        # interval cannot satisfy either the point-error or coverage gate on
        # every individual biological case.  This is tested on planted truth,
        # not the reference estimator's predictions.
        assay=load_exact(ROOT/(tag+'_batch_assay.json'))
        raw_by={}
        for m in assay['measurements']:
            pair=raw_by.setdefault((m['case_id'],m['group']),[0,0])
            pair[0]+=m['counts']['A1_unique']
            pair[1]+=m['counts']['A2_unique']
        for cid,t in truth.items():
            raw_errors=[]
            for g in ('case','control'):
                unique=raw_by[(cid,g)]
                raw=unique[0]/sum(unique)
                raw_errors.append(abs(raw-t[g+'_a1_fraction']))
            assert max(raw_errors)>.11, (tag,cid,'raw pooling wrongly solves this case')
        abs_errors=[]; intervals=0; widths=[]
        for r in got['records']:
            assert set(r)=={'case_id','case','control','delta_a1_fraction','shift_class'}
            t=truth[r['case_id']]
            for group in ('case','control'):
                panel=r[group]
                assert set(panel)=={'a1_fraction','interval_95','informative_molecules','batches_used','quality_flag'}
                assert panel['quality_flag']=='estimated'
                assert type(panel['batches_used']) is int and 2 <= panel['batches_used'] <= 3
                assert type(panel['informative_molecules']) is int and panel['informative_molecules']>=120*panel['batches_used']
                assert type(panel['a1_fraction']) in (float,int) and type(panel['a1_fraction']) is not bool
                assert_six_decimal_number(panel['a1_fraction'])
                real=t[group+'_a1_fraction'];v=panel['a1_fraction']
                abs_errors.append(abs(v-real))
                ci=panel['interval_95']
                assert type(ci) is list and len(ci)==2
                assert all(type(q) in (int,float) and math.isfinite(q) for q in ci)
                for bound in ci:assert_six_decimal_number(bound)
                assert all(0<=q<=1 for q in ci)
                assert ci[0]<=v<=ci[1]
                intervals+=ci[0]<=real<=ci[1]
                widths.append(ci[1]-ci[0])
            delta=round(r['case']['a1_fraction']-r['control']['a1_fraction'],6)
            assert r['delta_a1_fraction']==delta
            assert r['shift_class']==('A1_shift' if delta>=.12 else 'A2_shift' if delta<=-.12 else 'no_clear_shift')
        assert max(abs_errors)<.065
        assert sum(abs_errors)/len(abs_errors)<=.035
        assert intervals>=math.ceil(.85*len(abs_errors))
        assert max(widths)<.45 and sum(widths)/len(widths)<.32

def test_batch_assay_input_order_is_invariant(inputs):
    data=load_exact(ROOT/'sealed_batch_assay.json')
    result=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    data['calibrators'].reverse();data['measurements'].reverse()
    again=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    assert result==again

def test_batch_specific_calibration_is_used(inputs):
    data=load_exact(ROOT/'sealed_batch_assay.json')
    base=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    for r in data['calibrators']:
        if r['batch_id']=='lane_B':
            c=r['counts'];shift=min(35,c['A2_unique']//3)
            c['A1_unique']+=shift;c['A2_unique']-=shift
    changed=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    diffs=[abs(x['case']['a1_fraction']-y['case']['a1_fraction']) for x,y in zip(base['records'],changed['records'])]
    assert max(diffs)>.004

def test_batch_corrupted_calibration_replicate_is_robust(inputs):
    data=load_exact(ROOT/'sealed_batch_assay.json')
    base=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    for r in data['calibrators']:
        if r['known_a1_fraction']==.5 and r['replicate']==0:
            a=r['counts']['A1_unique'];b=r['counts']['A2_unique']
            r['counts']['A1_unique']=a+b;r['counts']['A2_unique']=0
    changed=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    errors=[abs(x['case']['a1_fraction']-y['case']['a1_fraction']) for x,y in zip(base['records'],changed['records'])]
    assert max(errors)<.03

def test_batch_two_informative_lanes_are_sufficient(inputs):
    """Two valid independent batches satisfy the documented identifiability rule."""
    assay=load_exact(ROOT/'sealed_batch_assay.json')
    assay['calibrators']=[r for r in assay['calibrators'] if r['batch_id']!='lane_C']
    assay['measurements']=[r for r in assay['measurements'] if r['batch_id']!='lane_C']
    result=run_agent(inputs[0],inputs[1],batch_assay_override=assay)['batch_inference.json']
    assert len(result['records'])==len(inputs[0]['cases'])
    for r in result['records']:
        for group in ('case','control'):
            panel=r[group]
            assert panel['quality_flag']=='estimated'
            assert type(panel['batches_used']) is int and panel['batches_used']==2
            assert panel['a1_fraction'] is not None and panel['interval_95'] is not None
            assert_six_decimal_number(panel['a1_fraction'])
    # A lane with only three distinct standard fractions is not calibrated,
    # even when the observed measurement count is large.
    standards=load_exact(ROOT/'sealed_batch_assay.json')
    standards['calibrators']=[r for r in standards['calibrators'] if not (
        r['batch_id']=='lane_A' and r['known_a1_fraction'] in (0.0,1.0))]
    result=run_agent(inputs[0],inputs[1],batch_assay_override=standards)['batch_inference.json']
    assert all(r[g]['batches_used']==2 for r in result['records'] for g in ('case','control'))
    # Molecules in the 119-unique-read lane must not count as a fourth
    # estimate (or be included in informative_molecules).
    minimum=load_exact(ROOT/'sealed_batch_assay.json')
    affected=('m01','case','lane_A')
    rows=[r for r in minimum['measurements'] if (r['case_id'],r['group'],r['batch_id'])==affected]
    assert rows
    for r in rows:
        r['counts']['A1_unique']=0
        r['counts']['A2_unique']=0
    rows[0]['counts']['A1_unique']=119
    result=run_agent(inputs[0],inputs[1],batch_assay_override=minimum)['batch_inference.json']
    panel=next(r for r in result['records'] if r['case_id']=='m01')['case']
    assert panel['quality_flag']=='estimated' and panel['batches_used']==2
    expected=sum(r['counts']['A1_unique']+r['counts']['A2_unique'] for r in minimum['measurements']
                 if r['case_id']=='m01' and r['group']=='case' and r['batch_id']!='lane_A')
    assert panel['informative_molecules']==expected


def test_batch_rounding_counts_and_unidentifiable_schema(inputs,outputs):
    """Assert numerical precision, identity/interval semantics and honest molecule accounting."""
    for result,kind in zip(outputs,('sealed','hidden')):
        assay=load_exact(ROOT/(kind+'_batch_assay.json'))
        molecule_totals={}
        for row in assay['measurements']:
            key=(row['case_id'],row['group'])
            molecule_totals[key]=molecule_totals.get(key,0)+row['counts']['A1_unique']+row['counts']['A2_unique']
        for r in result['batch_inference.json']['records']:
            assert type(r['case_id']) is str and r['case_id']
            for group in ('case','control'):
                v=r[group]
                assert set(v)=={'a1_fraction','interval_95','informative_molecules','batches_used','quality_flag'}
                assert type(v['informative_molecules']) is int and v['informative_molecules']>=0
                assert type(v['batches_used']) is int and 0<=v['batches_used']<=3
                assert v['informative_molecules']==molecule_totals[(r['case_id'],group)]
                assert v['quality_flag'] in ('estimated','unidentifiable')
                if v['quality_flag']=='estimated':
                    assert 2<=v['batches_used']<=3
                    assert_six_decimal_number(v['a1_fraction'])
                    assert type(v['interval_95']) is list and len(v['interval_95'])==2
                    assert all(0<=x<=1 for x in v['interval_95'])
                    assert all(assert_six_decimal_number(x) is None for x in v['interval_95'])
                    assert v['interval_95'][0]<=v['a1_fraction']<=v['interval_95'][1]
                else:
                    assert v['a1_fraction'] is None and v['interval_95'] is None
                    assert v['batches_used']<2
            a=r['case']['a1_fraction'];b=r['control']['a1_fraction']
            if a is None or b is None:
                assert r['delta_a1_fraction'] is None and r['shift_class']=='unidentifiable'
            else:
                assert_six_decimal_number(r['delta_a1_fraction'])
                assert r['delta_a1_fraction']==round(a-b,6)
                delta=r['delta_a1_fraction']
                assert r['shift_class']==('A1_shift' if delta>=.12 else 'A2_shift' if delta<=-.12 else 'no_clear_shift')


def test_batch_neutral_controls_require_measurement_response(inputs,outputs):
    """Every graded case must respond to changed biological case RNA evidence."""
    for kind,cases,obs in (('sealed',inputs[0],inputs[1]),('hidden',inputs[2],inputs[3])):
        assay=load_exact(ROOT/(kind+'_batch_assay.json'))
        truth={r['case_id']:r for r in load_exact(ROOT/(kind+'_batch_truth.json'))['records']}
        baseline=outputs[0 if kind=='sealed' else 1]['batch_inference.json']
        altered=copy.deepcopy(assay)
        for r in altered['measurements']:
            if r['group']=='case':
                cnt=r['counts'];n=cnt['A1_unique']+cnt['A2_unique']
                moved=min(cnt['A2_unique'],round(.30*n))
                cnt['A1_unique']+=moved;cnt['A2_unique']-=moved
        observed=run_agent(cases,obs,batch_assay_override=altered)['batch_inference.json']
        first={r['case_id']:r for r in baseline['records']}
        second={r['case_id']:r for r in observed['records']}
        assert set(first)==set(second)==set(truth)
        for cid in first:
            assert abs(truth[cid]['case_a1_fraction']-truth[cid]['control_a1_fraction'])>=.13-1e-6
            for group in ('case','control'):
                before=first[cid][group];after=second[cid][group]
                assert before['quality_flag']==after['quality_flag']=='estimated'
                assert set(after)=={'a1_fraction','interval_95','informative_molecules',
                                    'batches_used','quality_flag'}
                assert type(after['informative_molecules']) is int
                assert type(after['batches_used']) is int
                assert 2<=after['batches_used']<=3
                assert_six_decimal_number(after['a1_fraction'])
                interval=after['interval_95']
                assert type(interval) is list and len(interval)==2
                assert all(assert_six_decimal_number(x) is None for x in interval)
                assert 0<=interval[0]<=after['a1_fraction']<=interval[1]<=1
                assert interval[1]-interval[0]<.45
                if group=='case':
                    # A flat/copied answer cannot respond in every independent case.
                    assert after['a1_fraction']>before['a1_fraction']+.02
                else:
                    # Unchanged control evidence must remain biologically accurate;
                    # shared-uncertainty methods may modify interval widths.
                    assert after['informative_molecules']==before['informative_molecules']
                    assert after['batches_used']==before['batches_used']
                    assert abs(after['a1_fraction']-truth[cid]['control_a1_fraction'])<.065
            a=second[cid]['case']['a1_fraction'];b=second[cid]['control']['a1_fraction']
            delta=round(a-b,6)
            assert second[cid]['delta_a1_fraction']==delta
            assert second[cid]['shift_class']==('A1_shift' if delta>=.12 else
                'A2_shift' if delta<=-.12 else 'no_clear_shift')


def test_batch_missing_ids_and_single_lane_are_unidentifiable(inputs):
    """Missing sample IDs are invalid; one informative lane gives null estimates."""
    assay=load_exact(ROOT/'sealed_batch_assay.json')
    invalid=copy.deepcopy(assay)
    invalid['measurements']=[r for r in invalid['measurements'] if r['case_id']!='m01']
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],batch_assay_override=invalid)
    alone=copy.deepcopy(assay)
    alone['calibrators']=[r for r in alone['calibrators'] if r['batch_id']=='lane_A']
    alone['measurements']=[r for r in alone['measurements'] if r['batch_id']=='lane_A']
    out=run_agent(inputs[0],inputs[1],batch_assay_override=alone)['batch_inference.json']
    for row in out['records']:
        assert row['delta_a1_fraction'] is None and row['shift_class']=='unidentifiable'
        for group in ('case','control'):
            panel=row[group]
            assert panel['quality_flag']=='unidentifiable'
            assert panel['a1_fraction'] is None and panel['interval_95'] is None
            assert type(panel['batches_used']) is int and panel['batches_used']<=1


def test_batch_unidentifiable_and_cross_sample_id_validation(inputs):
    data=load_exact(ROOT/'sealed_batch_assay.json')
    for r in data['calibrators']:
        if r['batch_id'] in ('lane_A','lane_C'):
            r['counts']['A1_unique']=0
            r['counts']['A2_unique']=400
    got=run_agent(inputs[0],inputs[1],batch_assay_override=data)['batch_inference.json']
    assert all(r['shift_class']=='unidentifiable' for r in got['records'])
    assert all(r[g]['a1_fraction'] is None and r[g]['quality_flag']=='unidentifiable' for r in got['records'] for g in ('case','control'))
    # Data changes that introduce a new case must be rejected by the full CLI.
    data=load_exact(ROOT/'sealed_batch_assay.json')
    data['measurements'][0]['case_id']='unmatched-case'
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],batch_assay_override=data)
    invalid_count=load_exact(ROOT/'sealed_batch_assay.json')
    invalid_count['measurements'][0]['counts']['A1_unique']=-1
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],batch_assay_override=invalid_count)
    invalid_fraction=load_exact(ROOT/'sealed_batch_assay.json')
    invalid_fraction['calibrators'][0]['known_a1_fraction']=1.25
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],batch_assay_override=invalid_fraction)
    invalid_panel=load_exact(ROOT/'sealed_batch_assay.json')
    invalid_panel['measurements'][0]['group']='other'
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],batch_assay_override=invalid_panel)


def test_phased_multivariate_mixture_truth_and_intervals(outputs):
    """Independent seeded truth for two disjoint experiments, not oracle values."""
    for tag,result in zip(('sealed','hidden'),outputs):
        assay=load_exact(ROOT/(tag+'_phased_assay.json'))
        truth={r['case_id']:r for r in load_exact(ROOT/(tag+'_phased_truth.json'))['records']}
        got=result['phased_rna.json']
        assert set(got)=={'assay_id','records'}
        assert got['assay_id']==assay['assay_id']
        rows=got['records']
        # The sealed/hidden panels have fully replicated pure and mixed
        # standards in each lane. Eligibility therefore follows the per-panel
        # >=150 informative-molecule condition, independent of solver outputs.
        by_panel={}
        for measurement in assay['measurements']:
            k=(measurement['case_id'],measurement['group'])
            informative=sum(n for name,n in measurement['counts'].items() if name!='unassigned')
            by_panel.setdefault(k,{}).setdefault(measurement['batch_id'],0)
            by_panel[k][measurement['batch_id']]+=informative
        assert len(rows)==len(truth)
        assert [r['case_id'] for r in rows]==sorted(truth)
        # Verify each haplotype independently. Pooling hap0 and hap1
        # errors would conceal an inaccurate allele behind the other one.
        error={'hap0':[],'hap1':[]}
        covered={'hap0':0,'hap1':0}
        width={'hap0':[],'hap1':[]}
        for row in rows:
            t=truth[row['case_id']]
            targets=[t[group][hap+'_a1_fraction']
                     for group in ('case','control') for hap in ('hap0','hap1')]
            assert max(abs(v-.5) for v in targets)>=.17, (
                'Fixture must discriminate every case from data-independent central intervals')
            assert set(row)=={'case_id','case','control','variant_linked_shift','linkage_class'}
            for grp in ('case','control'):
                panel=row[grp]
                assert set(panel)=={'hap0_a1_fraction','hap1_a1_fraction','allelic_contrast','hap0_interval_95','hap1_interval_95','informative_molecules','batches_used','quality_flag'}
                eligible=[n for n in by_panel[(row['case_id'],grp)].values() if n>=150]
                assert type(panel['informative_molecules']) is int
                assert panel['informative_molecules']==sum(eligible), 'Retained informative molecule total must be exact'
                assert type(panel['batches_used']) is int
                assert panel['batches_used']==len(eligible), 'Count all and only qualified calibrated lanes'
                assert len(eligible)>=2
                assert panel['quality_flag']=='estimated'
                assert panel['allelic_contrast']==round(panel['hap1_a1_fraction']-panel['hap0_a1_fraction'],6)
                for hap in ('hap0','hap1'):
                    f=panel[hap+'_a1_fraction']
                    assert isinstance(f,(float,int)) and not isinstance(f,bool) and math.isfinite(f) and 0<=f<=1
                    assert Decimal(str(f))==Decimal(str(f)).quantize(Decimal('.000001'))
                    iv=panel[hap+'_interval_95']
                    assert type(iv) is list and len(iv)==2
                    assert all(type(x) in (int,float) and math.isfinite(x) for x in iv)
                    assert all(Decimal(str(x))==Decimal(str(x)).quantize(Decimal('.000001')) for x in iv)
                    assert 0<=iv[0]<=f<=iv[1]<=1
                    target=t[grp][hap+'_a1_fraction']
                    error[hap].append(abs(f-target))
                    covered[hap]+=int(iv[0]<=target<=iv[1])
                    width[hap].append(iv[1]-iv[0])
            c=row['case']['allelic_contrast'];d=row['control']['allelic_contrast']
            assert row['variant_linked_shift']==round(c-d,6)
            classification='hap1_A1_increased' if c-d>=.18 else 'hap1_A1_decreased' if c-d<=-.18 else 'no_clear_linkage_shift'
            assert row['linkage_class']==classification
        for hap in ('hap0','hap1'):
            # Robust to one stochastic near-miss on an otherwise accurate
            # independent deconvolution. Strict 0.08 remains required on all
            # but at most one panel per haplotype and fixture; the exception
            # cannot exceed 0.09. The aggregate MAE and unexpanded interval
            # coverage/width constraints remain separate hard gates.
            assert sum(error[hap])/len(error[hap])<=.04, f'{tag}: {hap} MAE'
            assert max(error[hap])<.09, f'{tag}: {hap} worst-point hard limit'
            assert sum(e>=.08 for e in error[hap])<=1, f'{tag}: {hap} repeated errors >=0.08'
            assert covered[hap]>=math.ceil(.85*len(error[hap])), f'{tag}: {hap} coverage'
            assert sum(width[hap])/len(width[hap])<.28, f'{tag}: {hap} mean interval width'
            assert max(width[hap])<.45, f'{tag}: {hap} max interval width'
        # A single interval copied onto each case/group/haplotype cannot
        # masquerade as data-calibrated uncertainty, even where all its
        # endpoints happen to contain the withheld truth.
        bounds={tuple(row[group][hap+'_interval_95']) for row in rows
                for group in ('case','control') for hap in ('hap0','hap1')}
        assert len(bounds)>=len(rows), 'Intervals must reflect individual measurements' 


def test_phased_mixture_input_permutations(inputs,outputs):
    import random
    rng=random.Random(10391)
    data=load_exact(ROOT/'sealed_phased_assay.json')
    data['standards']=list(reversed(data['standards']))
    data['measurements']=list(reversed(data['measurements']))
    for r in data['measurements']:
        r['counts']=dict(reversed(list(r['counts'].items())))
    got=run_agent(inputs[0],inputs[1],phased_assay_override=data)['phased_rna.json']
    assert got==outputs[0]['phased_rna.json']


def test_phased_two_valid_batches_suffice(inputs):
    data=load_exact(ROOT/'sealed_phased_assay.json')
    data['standards']=[r for r in data['standards'] if r['batch_id']!='lane_C']
    data['measurements']=[r for r in data['measurements'] if r['batch_id']!='lane_C']
    got=run_agent(inputs[0],inputs[1],phased_assay_override=data)['phased_rna.json']
    for r in got['records']:
        for g in ('case','control'):
            assert r[g]['batches_used']==2 and r[g]['quality_flag']=='estimated'
            assert r[g]['hap0_a1_fraction'] is not None
            assert r[g]['hap1_a1_fraction'] is not None


def test_phased_calibration_mix_response(inputs,outputs):
    data=load_exact(ROOT/'sealed_phased_assay.json')
    # Change *known* standard mixture assignments, not just observables.
    # A correctly calibrated deconvolution cannot ignore those ratios.
    for r in data['standards']:
        if r['batch_id'] in ('lane_A','lane_B') and r['standard_id'] not in ('pure_0','pure_1','pure_2','pure_3'):
            original=r['known_product_fractions']
            changed={k:v*(1.35 if k.endswith('A1') else .65) for k,v in original.items()}
            total=sum(changed.values())
            r['known_product_fractions']={k:v/total for k,v in changed.items()}
    result=run_agent(inputs[0],inputs[1],phased_assay_override=data)['phased_rna.json']
    deltas=[abs(r['case']['hap1_a1_fraction']-s['case']['hap1_a1_fraction']) for r,s in zip(result['records'],outputs[0]['phased_rna.json']['records'])]
    assert max(deltas)>.004,'Mixed-standard calibration must influence estimates'


def test_phased_joint_allele_splice_read_response(inputs,outputs):
    # A fixed haplotype answer or case=control copy must fail EVERY graded case.
    assay=load_exact(ROOT/'sealed_phased_assay.json')
    for row in assay['measurements']:
        if row['group']=='case':
            moved=min(80,row['counts']['h1_A2'])
            row['counts']['h1_A2']-=moved
            row['counts']['h1_A1']+=moved
    changed=run_agent(inputs[0],inputs[1],phased_assay_override=assay)['phased_rna.json']
    base={r['case_id']:r for r in outputs[0]['phased_rna.json']['records']}
    after={r['case_id']:r for r in changed['records']}
    assert set(base)==set(after)
    for cid in base:
        first=base[cid]['case'];second=after[cid]['case']
        assert first['quality_flag']==second['quality_flag']=='estimated'
        assert second['hap1_a1_fraction']>first['hap1_a1_fraction']+.03


def test_phased_single_lane_and_missing_ids_are_not_guessed(inputs):
    data=load_exact(ROOT/'sealed_phased_assay.json')
    one=copy.deepcopy(data)
    one['standards']=[r for r in one['standards'] if r['batch_id']=='lane_A']
    one['measurements']=[r for r in one['measurements'] if r['batch_id']=='lane_A']
    result=run_agent(inputs[0],inputs[1],phased_assay_override=one)['phased_rna.json']
    for r in result['records']:
        assert r['variant_linked_shift'] is None and r['linkage_class']=='unidentifiable'
        for g in ('case','control'):
            assert r[g]['quality_flag']=='unidentifiable'
            assert r[g]['hap0_a1_fraction'] is None
            assert r[g]['hap1_a1_fraction'] is None
            assert_phased_no_call(r[g])
            assert r[g]['batches_used']<=1
    invalid=copy.deepcopy(data)
    invalid['measurements']=[r for r in invalid['measurements'] if r['case_id']!='m03']
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],phased_assay_override=invalid)
    dup=copy.deepcopy(data)
    dup['measurements'].append(copy.deepcopy(dup['measurements'][0]))
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],phased_assay_override=dup)


def test_phased_single_contaminated_standard_robustness(inputs,outputs):
    data=load_exact(ROOT/'sealed_phased_assay.json')
    r=next(r for r in data['standards'] if r['batch_id']=='lane_A' and r['standard_id']=='a1_rich' and r['replicate_id']=='r1')
    counts=r['counts'];total=sum(counts.values())
    counts.update({k:0 for k in counts})
    counts['unphased_A2']=total
    got=run_agent(inputs[0],inputs[1],phased_assay_override=data)['phased_rna.json']
    diffs=[abs(r['case']['hap0_a1_fraction']-s['case']['hap0_a1_fraction']) for r,s in zip(got['records'],outputs[0]['phased_rna.json']['records'])]
    assert max(diffs)<.05


def test_phased_pooled_and_fixed_answer_baselines_do_not_solve():
    """Explicit baseline audit: raw allele-category counts and constant 0.5."""
    for tag in ('sealed','hidden'):
        data=load_exact(ROOT/(tag+'_phased_assay.json'))
        truth={r['case_id']:r for r in load_exact(ROOT/(tag+'_phased_truth.json'))['records']}
        grouped={}
        for r in data['measurements']:
            key=(r['case_id'],r['group'])
            out=grouped.setdefault(key,{'hap0': [0,0],'hap1':[0,0]})
            c=r['counts']
            for hap,tag2 in [('hap0','h0'),('hap1','h1')]:
                out[hap][0]+=c[tag2+'_A1'];out[hap][1]+=c[tag2+'_A2']
        naive={'hap0':[],'hap1':[]};flat={'hap0':[],'hap1':[]}
        for (cid,group),res in grouped.items():
            for hap in ('hap0','hap1'):
                t=truth[cid][group][hap+'_a1_fraction']
                x=res[hap];naive[hap].append(abs(x[0]/sum(x)-t));flat[hap].append(abs(.5-t))
        # Demonstrate that neither simplistic method meets the quality
        # threshold on BOTH haplotypes of each independent experiment.
        assert max(sum(naive[h])/len(naive[h]) for h in ('hap0','hap1'))>.04
        assert all(sum(flat[h])/len(flat[h])>.10 for h in ('hap0','hap1'))
        # Do not include "free" individual graded cases where the raw-read
        # shortcut clears all four per-case checks without doing calibration.
        for cid in truth:
            case_errors=[]
            for grp in ('case','control'):
                for hap in ('hap0','hap1'):
                    counts=grouped[(cid,grp)][hap]
                    raw=counts[0]/sum(counts)
                    target=truth[cid][grp][hap+'_a1_fraction']
                    case_errors.append(abs(raw-target))
            assert max(case_errors)>=.08 or sum(case_errors)/len(case_errors)>.04, (
                f'{tag} {cid}: uncalibrated raw reads solve all four case checks')


def test_phased_duplicate_calibration_or_bad_molecule_types_rejected(inputs):
    assay=load_exact(ROOT/'sealed_phased_assay.json')
    duplicate=copy.deepcopy(assay)
    duplicate['standards'].append(copy.deepcopy(duplicate['standards'][0]))
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],phased_assay_override=duplicate)
    wrong_type=copy.deepcopy(assay)
    wrong_type['measurements'][0]['counts']['h0_A1']=3.0
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):
        run_agent(inputs[0],inputs[1],phased_assay_override=wrong_type)


def test_phased_absent_mixed_standards_is_unidentifiable(inputs):
    assay=load_exact(ROOT/'sealed_phased_assay.json')
    assay['standards']=[r for r in assay['standards'] if r['standard_id'].startswith('pure_')]
    got=run_agent(inputs[0],inputs[1],phased_assay_override=assay)['phased_rna.json']
    for row in got['records']:
        assert row['variant_linked_shift'] is None
        for group in ('case','control'):
            assert_phased_no_call(row[group])
            assert row[group]['batches_used']<2
    # Six distinct IDs advertising one recipe do not meet the requirement
    # of four distinct known compositions; pure controls do not fix this.
    duplicate_recipe=load_exact(ROOT/'sealed_phased_assay.json')
    first=next(r['known_product_fractions'] for r in duplicate_recipe['standards']
               if not r['standard_id'].startswith('pure_'))
    for row in duplicate_recipe['standards']:
        if not row['standard_id'].startswith('pure_'):
            row['known_product_fractions']=dict(first)
    got=run_agent(inputs[0],inputs[1],phased_assay_override=duplicate_recipe)['phased_rna.json']
    for row in got['records']:
        for group in ('case','control'):
            assert_phased_no_call(row[group])
    # Three-positive mixtures must be usable without either four-positive
    # recipe. With four distinct recipes, panels retain calibrated batches.
    three_product=load_exact(ROOT/'sealed_phased_assay.json')
    three_product['standards']=[r for r in three_product['standards']
        if r['standard_id'].startswith('pure_') or r['standard_id'] in
        {'a1_rich','a2_rich','hap0_rich','hap1_rich'}]
    got=run_agent(inputs[0],inputs[1],phased_assay_override=three_product)['phased_rna.json']
    assert all(row[g]['quality_flag']=='estimated' and row[g]['batches_used']>=2
               for row in got['records'] for g in ('case','control'))


def test_phased_two_informative_batches_can_still_be_unidentifiable(inputs):
    """Absent hap0 information must not be confused with too few good batches.

    Retain the valid mixed calibrators and two measurement batches. Replace
    only one case panel with essentially hap1-specific reads. Both batches
    remain informative by the stated molecule and calibration criteria, but
    hap0 A1/A2 cannot be estimated reliably from this evidence.
    """
    data=load_exact(ROOT/'sealed_phased_assay.json')
    target=[r for r in data['measurements'] if r['case_id']=='m03' and r['group']=='case'
            and sum(r['counts'].values())-r['counts']['unassigned']>=150]
    assert len({r['batch_id'] for r in target})==2
    expected=0
    for row in target:
        counts=row['counts']
        total=sum(counts.values())
        row['counts']={key:(total-40 if key=='h1_A1' else 40 if key=='unassigned' else 0)
                       for key in counts}
        expected+=total-40
    result=run_agent(inputs[0],inputs[1],phased_assay_override=data)['phased_rna.json']
    record=next(r for r in result['records'] if r['case_id']=='m03')
    panel=record['case']
    assert_phased_no_call(panel)
    assert type(panel['batches_used']) is int and panel['batches_used']==2
    assert type(panel['informative_molecules']) is int and panel['informative_molecules']==expected
    for field in ('hap0_a1_fraction','hap1_a1_fraction','allelic_contrast',
                  'hap0_interval_95','hap1_interval_95'):
        assert panel[field] is None
    assert record['control']['quality_flag']=='estimated'
    assert record['variant_linked_shift'] is None
    assert record['linkage_class']=='unidentifiable'


# Additional independent paired-donor cross-assay scientific challenge.
def test_longitudinal_paired_truth_and_calibrated_intervals(outputs):
    for tag,out in zip(('sealed','hidden'),outputs):
        truth={r['case_id']:r for r in load_exact(ROOT/(tag+'_longitudinal_truth.json'))['records']}
        observed_assay=load_exact(ROOT/(tag+'_longitudinal_assay.json'))
        molecule_totals={}
        for lib in observed_assay['measurements']:
            key=(lib['case_id'],lib['donor_id'])
            molecule_totals[key]=molecule_totals.get(key,0)+sum(
                n for signal,n in lib['counts'].items() if signal!='unassigned')
        actual=out['longitudinal.json']
        assert set(actual)=={'assay_id','records'}
        assert actual['assay_id']=='EDA-paired-perturbation-v1'
        assert [r['case_id'] for r in actual['records']]==sorted(truth)
        assert len(actual['records'])==len(truth)
        errors=[];d_errors=[];cover=[];d_cover=[];widths=[];d_widths=[]
        mean_values=set()
        for r in actual['records']:
            assert set(r)=={'case_id','donors','mean_treatment_time_effect','interval_95','effect_class','donors_used','sample_count','quality_flag'}
            assert r['quality_flag']=='estimated' and r['donors_used']==6 and r['sample_count']==24
            assert type(r['donors_used']) is int and type(r['sample_count']) is int
            assert [d['donor_id'] for d in r['donors']]==[f'D{k:02d}' for k in range(1,7)]
            true=truth[r['case_id']]
            f=r['mean_treatment_time_effect']
            assert type(f) in (float,int) and math.isfinite(f) and -4<=f<=4
            assert Decimal(str(f))==Decimal(str(f)).quantize(Decimal('.000001'))
            assert abs(f-sum(d['treatment_time_effect'] for d in r['donors'])/6)<.000002
            errors.append(abs(f-true['mean_treatment_time_effect']))
            mean_values.add(f)
            iv=r['interval_95']
            assert type(iv) is list and len(iv)==2 and -4<=iv[0]<=f<=iv[1]<=4
            assert all(type(x) in (int,float) and math.isfinite(x) for x in iv)
            assert all(Decimal(str(x))==Decimal(str(x)).quantize(Decimal('.000001')) for x in iv)
            widths.append(iv[1]-iv[0]);cover.append(iv[0]<=true['mean_treatment_time_effect']<=iv[1])
            expected='hap1_treatment_increase' if f>=.12 else 'hap1_treatment_decrease' if f<=-.12 else 'no_clear_interaction'
            assert r['effect_class']==expected
            target={d['donor_id']:d for d in true['donors']}
            for d in r['donors']:
                assert set(d)=={'donor_id','treatment_time_effect','interval_95','informative_molecules','quality_flag'}
                assert d['quality_flag']=='estimated' and type(d['informative_molecules']) is int and d['informative_molecules']>=2000
                assert d['informative_molecules']==molecule_totals[(r['case_id'],d['donor_id'])]
                est=d['treatment_time_effect'];d_iv=d['interval_95']
                assert type(est) in (float,int) and math.isfinite(est)
                assert type(d_iv) is list and len(d_iv)==2 and -4<=d_iv[0]<=est<=d_iv[1]<=4
                assert all(type(x) in (int,float) and math.isfinite(x) for x in d_iv)
                assert all(Decimal(str(x))==Decimal(str(x)).quantize(Decimal('.000001')) for x in d_iv)
                assert Decimal(str(est))==Decimal(str(est)).quantize(Decimal('.000001'))
                d_errors.append(abs(est-target[d['donor_id']]['treatment_time_effect']))
                d_cover.append(d_iv[0]<=target[d['donor_id']]['treatment_time_effect']<=d_iv[1]);d_widths.append(d_iv[1]-d_iv[0])
        # Broad statistical boundaries with visible headroom for correct solvers;
        # *separate* donor and case metrics prevent pooled-mean shortcuts.
        assert len(mean_values)>=4
        # Every graded case defeats a constant zero effect and its fixed
        # +/-0.15 case interval; accuracy is still checked against truth.
        assert all(abs(t['mean_treatment_time_effect'])>.15 for t in truth.values())
        # Cases are prospectively designed so fixed +/-0.25 guesses fail the
        # published case-error point threshold before counterfactual tests.
        assert all(abs(abs(t['mean_treatment_time_effect'])-.25)>.105 for t in truth.values())
        # Every graded donor must also require inference; the zero-centred
        # shortcut [-.15,.15] must fail even before aggregate error gates.
        assert all(abs(d['treatment_time_effect'])>.15 for t in truth.values() for d in t['donors'])
        assert sum(errors)/len(errors)<.06
        # Nearest-rank floor of 0.9*(n-1), zero-based (lower-order statistic).
        assert sorted(errors)[int(.9*(len(errors)-1))]<.105
        assert sum(d_errors)/len(d_errors)<.115
        assert sum(cover)/len(cover)>=.80
        assert sum(d_cover)/len(d_cover)>=.75
        assert sum(widths)/len(widths)<.31
        assert sum(d_widths)/len(d_widths)<.48
        assert len({tuple(r['interval_95']) for r in actual['records']})>=len(actual['records'])//3


def test_longitudinal_measurement_order_invariance(inputs,outputs):
    assay=load_exact(ROOT/'sealed_longitudinal_assay.json')
    assay['measurements'].reverse()
    after=run_agent(inputs[0],inputs[1],longitudinal_assay_override=assay)
    assert after['longitudinal.json']==outputs[0]['longitudinal.json']


def test_longitudinal_reacts_to_shifted_treatment_junction_calls(inputs,outputs):
    # Perturb EVERY sealed case in a single agent invocation. This requires
    # case-specific treatment/time inference, not a fixed value or global mean.
    assay=load_exact(ROOT/'sealed_longitudinal_assay.json')
    for row in assay['measurements']:
        if row['condition']=='treated' and row['time']=='followup':
            row['counts']['h1_A1'],row['counts']['h1_A2']=row['counts']['h1_A2'],row['counts']['h1_A1']
            row['counts']['unphased_A1'],row['counts']['unphased_A2']=row['counts']['unphased_A2'],row['counts']['unphased_A1']
    changed=run_agent(inputs[0],inputs[1],longitudinal_assay_override=assay)['longitudinal.json']
    first={r['case_id']:r for r in outputs[0]['longitudinal.json']['records']}
    second={r['case_id']:r for r in changed['records']}
    assert set(first)==set(second)
    for cid in first:
        a=first[cid];b=second[cid]
        assert a['quality_flag']==b['quality_flag']=='estimated'
        assert type(b['mean_treatment_time_effect']) in (float,int)
        assert math.isfinite(b['mean_treatment_time_effect'])
        assert abs(b['mean_treatment_time_effect']-a['mean_treatment_time_effect'])>.03, (
            'Every case must respond materially to altered treated-followup junction evidence')


def test_longitudinal_missing_donor_and_insufficient_case(inputs):
    assay=load_exact(ROOT/'sealed_longitudinal_assay.json')
    cid=sorted({r['case_id'] for r in assay['measurements']})[0]
    donor='D01'
    assay['measurements']=[r for r in assay['measurements'] if not (r['case_id']==cid and r['donor_id']==donor and r['time']=='followup' and r['condition']=='treated')]
    result=run_agent(inputs[0],inputs[1],longitudinal_assay_override=assay)['longitudinal.json']
    row=next(r for r in result['records'] if r['case_id']==cid)
    assert row['quality_flag']=='estimated' and row['donors_used']==5 and row['sample_count']==23
    first=next(d for d in row['donors'] if d['donor_id']==donor)
    assert first['quality_flag']=='insufficient_data' and first['treatment_time_effect'] is None and first['interval_95'] is None
    true_total=sum(sum(n for signal,n in lib['counts'].items() if signal!='unassigned')
                   for lib in assay['measurements'] if lib['case_id']==cid and lib['donor_id']==donor)
    assert type(first['informative_molecules']) is int and first['informative_molecules']==true_total
    assay['measurements']=[r for r in assay['measurements'] if not (r['case_id']==cid and r['donor_id'] in ('D02','D03') and r['time']=='followup' and r['condition']=='treated')]
    result=run_agent(inputs[0],inputs[1],longitudinal_assay_override=assay)['longitudinal.json']
    row=next(r for r in result['records'] if r['case_id']==cid)
    assert row['donors_used']==3 and row['sample_count']==21 and row['quality_flag']=='insufficient_data'
    assert row['mean_treatment_time_effect'] is None and row['interval_95'] is None and row['effect_class']=='insufficient_data'


def test_longitudinal_invalid_library_identifiers_and_counts(inputs):
    assay=load_exact(ROOT/'sealed_longitudinal_assay.json')
    duplicate=copy.deepcopy(assay)
    duplicate['measurements'].append({**duplicate['measurements'][0],'replicate_id':'library2'})
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],longitudinal_assay_override=duplicate)
    bad=copy.deepcopy(assay)
    bad['measurements'][0]['counts']['h0_A1']=9.5
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],longitudinal_assay_override=bad)
    extra=copy.deepcopy(assay)
    extra['measurements'][0]['batch_id']='unknown_lane'
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],longitudinal_assay_override=extra)


def test_longitudinal_raw_batch_pooling_is_insufficient(outputs):
    # Scientific baseline: untreated read-category fractions are not a matched,
    # batch-corrected treatment-time effect. This check contains no solver truth.
    for tag,out in zip(('sealed','hidden'),outputs):
        raw=load_exact(ROOT/(tag+'_longitudinal_assay.json'))['measurements']
        truth={r['case_id']:r for r in load_exact(ROOT/(tag+'_longitudinal_truth.json'))['records']}
        panels={}
        for r in raw:
            c=r['counts']
            h0=c['h0_A1']/(c['h0_A1']+c['h0_A2'])
            h1=c['h1_A1']/(c['h1_A1']+c['h1_A2'])
            panels.setdefault((r['case_id'],r['donor_id']),{})[(r['condition'],r['time'])]=h1-h0
        from collections import defaultdict
        raw_est=defaultdict(list)
        for (cid,donor),v in panels.items():
            raw_est[cid].append(v['treated','followup']-v['treated','baseline']-v['untreated','followup']+v['untreated','baseline'])
        errors=[abs(sum(v)/len(v)-truth[cid]['mean_treatment_time_effect']) for cid,v in raw_est.items()]
        ours={r['case_id']:r['mean_treatment_time_effect'] for r in out['longitudinal.json']['records']}
        calibrated=[abs(ours[cid]-truth[cid]['mean_treatment_time_effect']) for cid in raw_est]
        # Only the published scientific MAE and 90th-percentile gates are
        # used; there is no hidden reference-vs-baseline margin.
        raw_mae=sum(errors)/len(errors)
        raw_p90=sorted(errors)[int(.9*(len(errors)-1))]
        assert raw_mae>=.06 or raw_p90>=.105, (tag,raw_mae,raw_p90)


def test_longitudinal_case_ids_match_dna(inputs):
    assay=load_exact(ROOT/'sealed_longitudinal_assay.json')
    wrong=copy.deepcopy(assay)
    wrong['measurements'][0]['case_id']='not_a_case'
    with pytest.raises(WorkflowNonzeroExit, match='Workflow exited nonzero'):run_agent(inputs[0],inputs[1],longitudinal_assay_override=wrong)


def test_longitudinal_depends_on_independent_mixed_standard_calibration(inputs,outputs):
    # Holding donor read counts fixed, perturb all lane_B mixed calibration
    # readouts. A fixed correction factor or raw-count approach cannot match.
    standards=load_exact(ROOT/'sealed_phased_assay.json')
    for row in standards['standards']:
        known=row['known_product_fractions']
        if row['batch_id']=='lane_B' and sum(v>0 for v in known.values())>=3:
            counts=row['counts']
            counts['h1_A1'],counts['h1_A2']=counts['h1_A2'],counts['h1_A1']
            counts['unphased_A1'],counts['unphased_A2']=counts['unphased_A2'],counts['unphased_A1']
    altered=run_agent(inputs[0],inputs[1],phased_assay_override=standards)['longitudinal.json']['records']
    previous=outputs[0]['longitudinal.json']['records']
    changes=[abs(x['mean_treatment_time_effect']-y['mean_treatment_time_effect']) for x,y in zip(altered,previous)]
    assert len(changes)==len(previous)
    assert any(x>1e-6 for x in changes), 'Altered calibration must change at least one rounded longitudinal effect'


# Independent exact posterior validation for the long-range signed phase graph.
def independent_phase_posterior(assay):
    from itertools import product
    rates={}
    for row in assay['calibration']:
        c=row['counts']; key=(row['batch_id'],row['known_relation'])
        prev=rates.setdefault(key,[0,0]);prev[0]+=c['cis'];prev[1]+=c['trans']
    rates={key:(a+.5)/(a+b+1) for key,(a,b) in rates.items()}
    result={}
    for rec in assay['cases']:
        nodes=rec['loci']; idx={s:i for i,s in enumerate(nodes)}
        edges=[]
        for e in rec['edges']:
            c=e['counts']; q0=rates[e['batch_id'],'same'];q1=rates[e['batch_id'],'opposite']
            score0=c['cis']*math.log(q0)+c['trans']*math.log(1-q0)
            score1=c['cis']*math.log(q1)+c['trans']*math.log(1-q1)
            edges.append((idx[e['left']],idx[e['right']],score0,score1))
        pool=[]
        for tail in product([0,1],repeat=len(nodes)-1):
            configuration=(0,)+tail
            score=sum(s1 if configuration[a]!=configuration[b] else s0 for a,b,s0,s1 in edges)
            pool.append((configuration,score))
        maxscore=max(s for _,s in pool)
        unnorm=[math.exp(s-maxscore) for _,s in pool]; z=sum(unnorm)
        marg=[sum(w*bits[j] for (bits,_),w in zip(pool,unnorm))/z for j in range(len(nodes))]
        result[rec['case_id']]=(marg,1.0/z)
    return result


def assert_phase_graph_output_precision(panel):
    """Enforce the agent-visible six-decimal float/nullable integer contract.

    JSON booleans compare equal to integers in Python and unrounded numbers can
    be numerically close to the oracle, so both cases need explicit type/scale
    checks rather than approximate scientific comparisons alone.
    """
    for row in panel['records']:
        posterior=row['map_posterior']
        assert type(posterior) is float, 'MAP probability must be a JSON float'
        assert_six_decimal_number(posterior)
        for site in row['loci']:
            p=site['p_opposite_anchor']
            assert type(p) is float, 'Marginal probability must be a JSON float'
            assert_six_decimal_number(p)
            bit=site['phase_bit']
            assert bit is None or (type(bit) is int and bit in (0,1)), 'Phase bit must be 0, 1, or null, never bool'


def test_phase_graph_scientific_posterior_and_schema(outputs):
    for kind,out in [('sealed',outputs[0]),('hidden',outputs[1])]:
        assay=load_exact(ROOT/f'{kind}_phase_graph_assay.json')
        oracle=independent_phase_posterior(assay)
        panel=out['phase_graph.json']
        assert_phase_graph_output_precision(panel)
        assert set(panel)=={'assay_id','records'} and panel['assay_id']==assay['assay_id']
        assert [r['case_id'] for r in panel['records']]==sorted(oracle)
        truth={r['case_id']:r['phase_bits'] for r in load_exact(ROOT/f'{kind}_phase_graph_truth.json')['records']}
        errors=[];brier=[];sharp=0;uncertain=0
        bycase={r['case_id']:r for r in assay['cases']}
        for row in panel['records']:
            assert set(row)=={'case_id','anchor','loci','map_posterior','quality_flag'}
            names=bycase[row['case_id']]['loci']
            assert row['anchor']==names[0] and [s['locus_id'] for s in row['loci']]==names
            assert type(row['map_posterior']) in (float,int) and 0<=row['map_posterior']<=1
            expected,pm=oracle[row['case_id']]
            assert abs(row['map_posterior']-pm)<0.025
            seen=[]
            for k,(site,refp) in enumerate(zip(row['loci'],expected)):
                assert set(site)=={'locus_id','phase_bit','p_opposite_anchor'}
                v=site['p_opposite_anchor'];assert type(v) in (float,int) and math.isfinite(v) and 0<=v<=1
                assert abs(v-refp)<0.03,(row['case_id'],site['locus_id'],v,refp)
                errors.append(abs(v-refp));brier.append((v-truth[row['case_id']][k])**2)
                expected_bit=0 if refp<=.05 else (1 if refp>=.95 else None)
                assert site['phase_bit'] is expected_bit or site['phase_bit']==expected_bit
                if expected_bit is None:uncertain+=1
                else:sharp+=1
                seen.append(expected_bit)
            assert row['quality_flag']==('resolved' if all(v is not None for v in seen) else 'partial')
        assert sum(errors)/len(errors)<.01
        assert sum(brier)/len(brier)<.08
        assert sharp>100 and uncertain>=5


def test_phase_graph_requires_real_global_inference(outputs):
    for output in outputs:
        rows=output['phase_graph.json']['records']
        confident=[v for row in rows for v in row['loci'] if v['phase_bit'] is not None]
        # A fixed all-zero phase guess is wrong for many true sites.
        assert sum(v['phase_bit']==1 for v in confident)>len(confident)*.23
        # The corpus includes uncertainty: fixed 0.5 and fixed intervals also fail.
        assert any(row['quality_flag']=='partial' for row in rows)
        assert max(row['map_posterior'] for row in rows)>.95
        assert min(row['map_posterior'] for row in rows)<.80
        # Each graded panel must discriminate against the no-inference shortcut
        # of predicting bit 0, probability 0, and MAP probability 1 everywhere.
        # A nearly certain all-zero case would otherwise pass its own oracle
        # comparison even though the corpus-level test remains difficult.
        for row in rows:
            assert max(s['p_opposite_anchor'] for s in row['loci'][1:])>.5, (row['case_id'], 'non-discriminating phase panel')


def test_phase_graph_input_order_invariance(inputs,outputs):
    a=load_exact(ROOT/'sealed_phase_graph_assay.json')
    a['cases'].reverse();a['calibration'].reverse()
    for r in a['cases']:r['edges'].reverse()
    changed=run_agent(inputs[0],inputs[1],phase_graph_override=a)['phase_graph.json']
    before=outputs[0]['phase_graph.json']
    assert [r['case_id'] for r in changed['records']]==[r['case_id'] for r in before['records']]
    for x,y in zip(changed['records'],before['records']):
        for v,w in zip(x['loci'],y['loci']):
            assert abs(v['p_opposite_anchor']-w['p_opposite_anchor'])<.00002
            assert v['phase_bit']==w['phase_bit']
        assert abs(x['map_posterior']-y['map_posterior'])<.00002


def test_phase_graph_responds_to_calibration_and_edge_data(inputs,outputs):
    a=load_exact(ROOT/'sealed_phase_graph_assay.json')
    first=sorted(a['cases'],key=lambda x:x['case_id'])[0]['case_id']
    before={r['case_id']:r for r in outputs[0]['phase_graph.json']['records']}[first]
    # Reverse a whole batch's calibration relation; keep the observations fixed.
    for row in a['calibration']:
        if row['batch_id']=='lane_A':row['known_relation']='opposite' if row['known_relation']=='same' else 'same'
    changed=run_agent(inputs[0],inputs[1],phase_graph_override=a)['phase_graph.json']
    delta=sum(abs(v['p_opposite_anchor']-w['p_opposite_anchor']) for v,w in zip(before['loci'],next(r for r in changed['records'] if r['case_id']==first)['loci']))
    assert delta>.10,'The allele graph must depend on lane-specific calibration'
    b=load_exact(ROOT/'sealed_phase_graph_assay.json')
    for edge in b['cases'][0]['edges']:
        edge['counts']['cis'],edge['counts']['trans']=edge['counts']['trans'],edge['counts']['cis']
    moved=run_agent(inputs[0],inputs[1],phase_graph_override=b)['phase_graph.json']
    orig={r['case_id']:r for r in outputs[0]['phase_graph.json']['records']}[b['cases'][0]['case_id']]
    new=next(r for r in moved['records'] if r['case_id']==b['cases'][0]['case_id'])
    assert any(abs(v['p_opposite_anchor']-w['p_opposite_anchor'])>.03 for v,w in zip(orig['loci'],new['loci']))


def test_phase_graph_unidentifiable_panels_and_bad_inputs(inputs):
    a=load_exact(ROOT/'sealed_phase_graph_assay.json')
    target=a['cases'][0]
    target['edges']=[]
    after=run_agent(inputs[0],inputs[1],phase_graph_override=a)['phase_graph.json']
    rec=next(r for r in after['records'] if r['case_id']==target['case_id'])
    assert rec['quality_flag']=='partial'
    assert rec['loci'][0]['phase_bit']==0 and rec['loci'][0]['p_opposite_anchor']==0
    assert abs(rec['map_posterior']-2**(1-len(rec['loci'])))<.00001
    assert all(x['phase_bit'] is None and x['p_opposite_anchor']==.5 for x in rec['loci'][1:])
    bad=load_exact(ROOT/'sealed_phase_graph_assay.json')
    bad['cases'][0]['edges'][0]['counts']['cis']=-3
    with pytest.raises(WorkflowNonzeroExit):
        run_agent(inputs[0],inputs[1],phase_graph_override=bad)
    bad=load_exact(ROOT/'sealed_phase_graph_assay.json')
    bad['cases'][0]['edges'][0]['batch_id']='unknown_batch'
    with pytest.raises(WorkflowNonzeroExit):
        run_agent(inputs[0],inputs[1],phase_graph_override=bad)
    # Cardinality-only 9-locus rejection: discard edges touching a removed
    # locus so a dangling endpoint cannot cause the nonzero exit instead.
    cardinality=load_exact(ROOT/'sealed_phase_graph_assay.json')
    first=cardinality['cases'][0]
    first['loci']=first['loci'][:9]
    retained=set(first['loci'])
    first['edges']=[e for e in first['edges']
                    if e['left'] in retained and e['right'] in retained]
    assert first['edges'] and all(e['left'] in retained and e['right'] in retained
                                  for e in first['edges'])
    with pytest.raises(WorkflowNonzeroExit):
        run_agent(inputs[0],inputs[1],phase_graph_override=cardinality)
    # An over-limit 13-locus input has no dangling edges or invalid batches;
    # only the documented maximum locus count is violated.
    too_many=load_exact(ROOT/'sealed_phase_graph_assay.json')
    first=too_many['cases'][0]
    first['loci'] += [f'extra_valid_locus_{i}' for i in range(13-len(first['loci']))]
    assert len(first['loci'])==13 and len(set(first['loci']))==13
    assert all(e['left'] in first['loci'] and e['right'] in first['loci'] for e in first['edges'])
    with pytest.raises(WorkflowNonzeroExit):
        run_agent(inputs[0],inputs[1],phase_graph_override=too_many)
    duplicate=load_exact(ROOT/'sealed_phase_graph_assay.json')
    duplicate['cases'][1]['edges'][0]['replicate_id']=duplicate['cases'][0]['edges'][0]['replicate_id']
    with pytest.raises(WorkflowNonzeroExit):
        run_agent(inputs[0],inputs[1],phase_graph_override=duplicate)
