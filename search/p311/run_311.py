"""Run the original CP-SAT peeling search for LS[77](3,4,311)."""
import argparse
import json
from pathlib import Path
from affine_conjecture import peel
from verify_affine_certificate import verify

parser = argparse.ArgumentParser()
parser.add_argument('--seconds', type=float, default=150)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--output', default='new_completion_311.json')
args = parser.parse_args()
result = peel(311, seconds=args.seconds, seed=args.seed)
path = Path(args.output)
if path.exists():
    raise FileExistsError(f'Refusing to overwrite {path}')
path.write_text(json.dumps(result, indent=2) + '\n')
print({k: v for k, v in result.items() if k != 'base_blocks'})
if 'base_blocks' in result:
    print(verify(path))
