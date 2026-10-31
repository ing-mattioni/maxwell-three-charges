"""Run the reviewed tests and certificates, retaining logs and source hashes.

Usage: python reproduce.py --workers 16 [--output NEW_DIRECTORY]
The output directory must not already exist. No repository files are replaced.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCES = ['certify.py', 'base_triangle.py', 'verify.py', 'interval.py', 'ad.py',
           'domains.py', 'provenance.py', 'closedform.py', 'reproduce.py',
           'tests/test_enclosures.py', 'tests/test_rigorous_arithmetic.py']


def hashes():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'certificates'/datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    result = dict(source_sha256=hashes(), python=sys.version, numpy=np.__version__, runs=[])
    jobs = [
        ('arithmetic', ['tests/test_rigorous_arithmetic.py'], None),
        ('enclosures', ['tests/test_enclosures.py'], None),
        ('n12', ['certify.py', '1.0', str(args.workers), '60'], 'cert_s1.0.json'),
        ('n10', ['certify.py', '1.0', str(args.workers), '80', '10'], 'cert_s1.0_n10.json'),
        ('base64', ['base_triangle.py', '64'], 'cert_base_triangle.json'),
        ('base48', ['base_triangle.py', '48'], 'cert_base_triangle.json'),
    ]
    for label, params, filename in jobs:
        directory = output/label
        directory.mkdir()
        command = [sys.executable, '-u', str(ROOT/params[0]), *params[1:]]
        print('Starting', label, flush=True)
        start = time.perf_counter()
        with (directory/'run.log').open('w', encoding='utf-8') as log:
            process = subprocess.run(command, cwd=directory, env=env,
                                     stdout=log, stderr=subprocess.STDOUT)
        row = dict(label=label, command=command, returncode=process.returncode,
                   seconds=time.perf_counter()-start)
        result['runs'].append(row)
        result['source_unchanged'] = result['source_sha256'] == hashes()
        valid = process.returncode == 0 and result['source_unchanged']
        if filename and (directory/filename).exists():
            cert = json.loads((directory/filename).read_text())
            row['certificate'] = cert
            valid = valid and cert['unresolved'] == 0
            for name, digest in cert['source_sha256'].items():
                valid = valid and result['source_sha256'][name] == digest
            if label in ('n12', 'n10'):
                tag = 's1.0' if label == 'n12' else 's1.0_n10'
                pending = np.load(directory/f'unresolved_{tag}.npy', allow_pickle=False)
                valid = valid and pending.shape == (0, 8)
                st, margins = cert['stats'], cert['worst_certified']
                closed = sum(st[k] for k in ('disc_struct','disc_outside','ver_int','ver_bdry'))
                valid = valid and st['boxes'] == 2*closed-cert['initial_boxes']
                valid = valid and margins['jac'] < 0 < margins['Ih']
                valid = valid and margins['turn'] < .5 and margins['W'] < .5
        elif filename:
            valid = False
        row['validated'] = bool(valid)
        (output/'reproduction.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        if not valid:
            raise SystemExit('Validation failed; inspect '+str(directory))
        print('Passed', label, 'in', round(row['seconds'], 2), 'seconds', flush=True)
    print('All tests and certificates passed:', output, flush=True)


if __name__ == '__main__':
    main()
