"""Small provenance record shared by the computational certificates."""
import hashlib
import platform
import sys
from pathlib import Path
import numpy as np


def metadata(source_files, **parameters):
    root = Path(__file__).resolve().parent
    return dict(
        schema_version=2,
        parameters=parameters,
        coverage='outward-rounded exact rational simplex cells',
        runtime=dict(python=sys.version, numpy=np.__version__, platform=platform.platform()),
        source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                       for name in sorted(set(source_files) | {'provenance.py'})},
    )
