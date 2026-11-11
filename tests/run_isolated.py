#!/usr/bin/env python3
"""Run one complete trusted test manifest in bounded subprocess groups.

Splitting pytest processes avoids sporadic long-lived runner stalls; it does
not skip or weaken any test. Every result is merged into a complete JUnit file.
"""
import ast
import copy
import os
import tempfile
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

HERE=Path(__file__).resolve().parent
SUITE=HERE/'test_workflow.py'
GROUP_SIZES=(20,13,1,1,4,9,11,8,5)

def list_tests():
    tree=ast.parse(SUITE.read_text())
    names=[x.name for x in tree.body if isinstance(x,ast.FunctionDef) and x.name.startswith('test_')]
    if len(names)!=sum(GROUP_SIZES) or len(set(names))!=len(names):
        raise RuntimeError('Unexpected test manifest')
    return names

def run(destination):
    names=list_tests()
    destination=Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    root=ET.Element('testsuites')
    suite=ET.SubElement(root,'testsuite',name='eda_scientific_verifier')
    counts={'tests':0,'failures':0,'errors':0,'skipped':0}
    started=time.monotonic();index=0;all_ok=True
    # New private cache per verifier attempt, never shared with the agent or
    # a different submission. Avoids 10 repeated sealed/hidden baseline runs.
    cache_dir=Path(tempfile.mkdtemp(prefix='run_',dir=str(destination.parent)))
    cache_dir.chmod(0o700)
    subprocess_env=dict(os.environ,EDA_VERIFIER_RUN_DIR=str(cache_dir),PYTEST_DISABLE_PLUGIN_AUTOLOAD='1')
    for number,size in enumerate(GROUP_SIZES,1):
        group=names[index:index+size];index+=size
        xml=destination.parent/f'junit_part_{number}.xml'
        xml.unlink(missing_ok=True)  # Never reuse a stale group report
        cmd=[sys.executable,'-m','pytest','-q','-p','no:cacheprovider',
              *[str(SUITE)+'::'+x for x in group],'--junitxml',str(xml)]
        status=0;reason=None
        try:
            completed=subprocess.run(cmd,cwd=str(HERE),env=subprocess_env,timeout=950)
            status=completed.returncode
        except subprocess.TimeoutExpired:
            status=124;reason=f'Isolated group {number} exceeded 950 seconds (fail-safe)'
        if xml.exists():
            try:
                tree=ET.parse(xml);base=tree.getroot()
                subs=[base] if base.tag=='testsuite' else list(base.findall('testsuite'))
                if len(subs)!=1:raise ValueError('JUnit must have one suite')
                found=subs[0]
                cases=found.findall('testcase')
                if len(cases)!=size:raise ValueError('JUnit group has wrong test count')
                for c in cases:suite.append(copy.deepcopy(c))
            except (ET.ParseError,ValueError) as exc:
                reason=f'Invalid group {number} JUnit: {exc}'
        else:
            reason=f'No JUnit XML for group {number}'
        if status or reason:
            all_ok=False
            # Keep evaluating other science stages after a failed assertion;
            # never award credit for a failed group.
            group_cases=cases if xml.exists() and not reason else []
            has_reported_error=any(c.find('failure') is not None or c.find('error') is not None for c in group_cases)
            if reason or (status and not has_reported_error):
                message=reason or f'Group {number} exited {status} without a failed JUnit testcase'
                tc=ET.SubElement(suite,'testcase',name=f'isolated_group_{number}_runner',time='0')
                ET.SubElement(tc,'error',message=message).text=message
    # Enforce the *documented* global limit on actual agent executions rather
    # than counting pytest test functions or arbitrarily partitioning budgets.
    try:
        invocations=[float(line) for line in (cache_dir/'invocations_seconds.log').read_text().splitlines()]
        agent_seconds=sum(invocations)
        (destination.parent/'invocation_accounting.json').write_text(__import__('json').dumps({
            'actual_invocations':len(invocations), 'agent_seconds':round(agent_seconds,6),
            'maximum_invocations':100,'maximum_agent_seconds':900,
            'per_invocation_limit_seconds':35},indent=2)+'\n')
        if len(invocations)>100 or agent_seconds>900:
            all_ok=False
            diagnostic=f'Agent execution budget exceeded: {len(invocations)} invocations, {agent_seconds:.2f}s'
            tc=ET.SubElement(suite,'testcase',name='aggregate_execution_budget',time='0')
            ET.SubElement(tc,'error',message=diagnostic).text=diagnostic
    except (FileNotFoundError,ValueError) as exc:
        all_ok=False
        diagnostic=f'No trustworthy invocation accounting: {exc}'
        tc=ET.SubElement(suite,'testcase',name='aggregate_execution_accounting',time='0')
        ET.SubElement(tc,'error',message=diagnostic).text=diagnostic
    cases=suite.findall('testcase')
    counts['tests']=len(cases)
    counts['errors']=sum(bool(c.findall('error')) for c in cases)
    counts['failures']=sum(bool(c.findall('failure')) for c in cases)
    counts['skipped']=sum(bool(c.findall('skipped')) for c in cases)
    suite.attrib.update({k:str(v) for k,v in counts.items()})
    suite.set('time',f'{time.monotonic()-started:.4f}')
    ET.ElementTree(root).write(destination,encoding='utf-8',xml_declaration=True)
    return 0 if all_ok and counts=={'tests':sum(GROUP_SIZES),'failures':0,'errors':0,'skipped':0} else 1

if __name__=='__main__':
    if len(sys.argv)!=2:sys.exit(2)
    try:sys.exit(run(sys.argv[1]))
    except Exception as exc:
        dest=Path(sys.argv[1]);dest.parent.mkdir(parents=True,exist_ok=True)
        root=ET.Element('testsuite',name='eda_scientific_verifier',tests='1',failures='0',errors='1',skipped='0')
        case=ET.SubElement(root,'testcase',name='verifier_runner_exception',time='0')
        ET.SubElement(case,'error',message=str(exc)).text=repr(exc)
        ET.ElementTree(root).write(dest,encoding='utf-8',xml_declaration=True)
        sys.exit(1)
