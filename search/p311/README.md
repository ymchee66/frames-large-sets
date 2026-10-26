# Search for the LS[77](3,4,311)

The large set LS[77](3,4,311) in `data/affine_large_sets.txt` was found with these programs.
They implement the same peeling method as `search/peelbt.py`, with one difference: the first
design is prescribed. It is the union of the orbit of {0,1,2,3} under AGL(1,311) and of the orbits
of the 4-subsets that contain the mean of their four elements. Each further design is found with
the OR-Tools CP-SAT solver, the search backtracks when no design can be found, and the remaining
orbits form the last design.

```
python3 -m pip install ortools
python3 run_311.py --seconds 150 --seed 1
```

The random seed and OR-Tools version of the original run were not recorded, so a run need not
reproduce the same large set or running time. A run that reaches its time limit is inconclusive;
try another seed or a larger time budget. `run_311.py` writes the large set found to a JSON file
and checks it with `verify_affine_certificate.py`, which uses only the Python standard library.
