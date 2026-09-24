# low75

Seventy-five genera: 25 bacteria, 25 archaea, and 25 viruses. Each genus contributes two species. The lower species taxid is simulated. The higher species taxid is the only sequence from that genus in the Kraken2 and Kaiju databases. Genomes, reads, and indexes are not committed.

Accessions were taken from the NCBI RefSeq `assembly_summary.txt` files dated 2026-09-23 (bacteria and archaea Last-Modified 2026-09-23, viral the same download) and NCBI taxdump `nodes.dmp` / `names.dmp`. The rule is in `select_accessions.py`. The table in this directory is the pin. Re-running the selector on a newer summary can change accessions.

The read budget matches `examples/heldout_genera`: `random.Random(42).lognormvariate(0, 1.5)` and 100000 paired fragments across all simulated genomes.

```bash
python examples/low75/run_example.py
```

`examples/low75half/` keeps every other genus in this table (38 of 75, because 75 is odd) and reuses these FASTA files.
