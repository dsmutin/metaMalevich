# low75half

Same simulation, databases, assembly, and reprofiling as `examples/low75`, on every other genus in that accession table. Seventy-five is odd, so 38 genera are kept and 37 are dropped. Both species of a kept genus stay: one in the metagenome and one in the Kraken2 and Kaiju databases.

FASTA files come from `python examples/low75/run_example.py --stage download`. This script does not download them again.

```bash
python examples/low75half/run_example.py
```

The read draw is a new `random.Random(42).lognormvariate(0, 1.5)` allocation of 100000 paired fragments across these genomes, not a slice of the 75-genome abundance table.
