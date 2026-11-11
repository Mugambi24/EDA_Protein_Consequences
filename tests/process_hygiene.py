#!/usr/bin/env python3
"""Trusted verifier-only cleanup of detached agent descendants.

All agent runs are demoted to uid/gid 65534 in this dedicated container.
This tool runs as root in the separate verifier container, never in the agent
container. It terminates *remaining* 65534 processes, including descendants
that escaped an immediate child's process group via setsid().
"""
import os
import signal
import sys
import time
from pathlib import Path

AGENT_UID = 65534


def active_agent_pids():
    active = set()
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            pid = int(entry.name)
            status = (entry / 'status').read_text()
            lines = {}
            for line in status.splitlines():
                if line.startswith(('Uid:', 'State:')):
                    k, value = line.split(':', 1)
                    lines[k] = value.strip()
            uids = [int(v) for v in lines['Uid'].split()[:4]]
            state = lines['State'].split()[0]
            if AGENT_UID in uids and state not in ('Z', 'X', 'x'):
                active.add(pid)
        except (OSError, ValueError, KeyError, IndexError, ProcessLookupError):
            # A short-lived process may disappear while its procfs entry is read.
            continue
    return active


def clean_agent_processes(seconds=5.0):
    if os.geteuid() != 0 or not Path('/proc').is_dir():
        print('Verifier process cleanup requires root and /proc', file=sys.stderr)
        return False
    deadline = time.monotonic() + seconds
    clean_scans = 0
    while time.monotonic() < deadline:
        active = active_agent_pids()
        if not active:
            clean_scans += 1
            if clean_scans >= 3:
                return True
        else:
            clean_scans = 0
            for pid in active:
                try:
                    # Check identity a second time to limit pid-reuse races.
                    if pid in active_agent_pids():
                        os.kill(pid, signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    pass
        time.sleep(.08)
    survivors = sorted(active_agent_pids())
    if survivors:
        print('Live untrusted agent processes remain:', survivors, file=sys.stderr)
    return not survivors


if __name__ == '__main__':
    sys.exit(0 if clean_agent_processes() else 1)
