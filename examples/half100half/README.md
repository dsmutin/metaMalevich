# half100half

Every other family in `examples/high100`, scored at family rank. One hundred families is even, so this keeps 50. Both genomes of a kept family stay: one simulated, one in the Kraken2 and Kaiju databases.

```bash
python examples/half100half/run_example.py
```

Run `examples/high100/run_example.py --stage download` first. The read draw is a new seed-42 lognormal allocation of 100000 paired fragments, not a slice of the 100-family table.
