# high100

One hundred families: 25 bacteria, 25 archaea, 25 viruses, and 25 eukaryotes (fungi and protozoa). Each family contributes two assemblies. The lower species taxid is simulated. The higher species taxid is the Kraken2 and Kaiju database genome.

In each domain, 13 families use two species of one genus and 12 use one species from each of two genera. Twenty-five is odd, so the split is 13 and 12 rather than 12.5. Scores roll to family rank. A call of the database species counts as the simulated family when that species is in the same family.

Eukaryotic assemblies are the smallest complete or chromosome RefSeq records that still form those pairs. The smallest pair sums to about 4.7 Mb and the largest to about 82 Mb. Larger families were not used.

Accessions are pinned from NCBI RefSeq `assembly_summary.txt` files (bacteria, archaea, and viral downloaded 2026-09-24; fungi and protozoa the same day) and the local NCBI taxdump. `select_accessions.py` is the rule. The TSV is the pin.

```bash
python examples/high100/run_example.py
```

The read budget is again 100000 paired fragments, seed 42. `examples/half100half/` keeps every other family (50 of 100).
