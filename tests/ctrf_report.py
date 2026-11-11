#!/usr/bin/env python3
"""Write schema-complete CTRF on successful and failed verifier runs.

The shell EXIT trap calls this tool unconditionally. Missing or invalid JUnit,
missing expected reporting data and failure markers create explicit failed tests.
Reward decisions never depend on this report or its contents.
"""
import json
import math
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path


def _elapsed_ms(raw):
    try:
        v = float(raw)
        if not math.isfinite(v) or v < 0:
            raise ValueError('negative/nonfinite duration')
        return round(v * 1000)
    except (ValueError, TypeError):
        raise ValueError('Invalid JUnit testcase duration: '+repr(raw))


def _fallback_start(path, stop):
    if path is not None and path.is_file():
        try:
            value = int(path.read_text().strip())
            if 0 < value <= stop:
                return value
        except (ValueError, OSError):
            pass
    return stop


def create(junit_path, report_path, marker_path, start_marker=None):
    stop = int(time.time() * 1000)
    start = _fallback_start(start_marker, stop)
    cases = []
    diagnostics = []
    try:
        root = ET.parse(junit_path).getroot()
        suites = [root] if root.tag == 'testsuite' else root.findall('testsuite')
        if not suites:
            raise ValueError('JUnit did not contain a testsuite')
        for suite in suites:
            for test in suite.findall('testcase'):
                duration = _elapsed_ms(test.get('time', '0'))
                bad = test.find('failure')
                err = test.find('error')
                skipped = test.find('skipped')
                row = {'name':test.get('name','unknown'), 'duration':duration}
                if bad is not None or err is not None:
                    node = bad if bad is not None else err
                    row.update(status='failed',message=(node.get('message') or (node.text or ''))[:2000])
                elif skipped is not None:
                    row.update(status='skipped')
                else:
                    row.update(status='passed')
                cases.append(row)
            declared = int(suite.get('tests', -1))
            if declared != len(suite.findall('testcase')):
                diagnostics.append('JUnit declared-count mismatch')
        if not cases:
            diagnostics.append('JUnit contains no test cases')
    except Exception as exc:
        diagnostics.append('Missing or invalid JUnit result: '+str(exc))
    if marker_path.is_file():
        marker = marker_path.read_text().strip()
        if marker:
            diagnostics.append(marker[:2000])
    for i, message in enumerate(diagnostics, 1):
        cases.append({'name':f'verifier_execution_error_{i}', 'status':'failed',
                      'duration':0, 'message':message})
    counts = {v:sum(x['status']==v for x in cases) for v in ('passed','failed','skipped')}
    result={'results':{
        'tool':{'name':'pytest'},
        'summary':{'tests':len(cases),'passed':counts['passed'],
                   'failed':counts['failed'],'skipped':counts['skipped'],
                   'pending':0,'other':0,'start':start,'stop':max(stop,start)},
        'tests':cases}}
    report_path.parent.mkdir(parents=True,exist_ok=True)
    report_path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

if __name__=='__main__':
    create(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),
           Path(sys.argv[4]) if len(sys.argv)>4 else None)
