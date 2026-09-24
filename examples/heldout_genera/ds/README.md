# Where the scored hypotheses are wrong

`analyze.py` reads the samovar10 benchmark and the assignment tables under each example's `work/reprofile/assignments/`. It writes TSV summaries, Altair charts, and a findings note per dataset. Half-strain charts are written to `examples/half_strains/ds/`. Re-run from the repository with `PYTHONPATH=src`:

```bash
python examples/heldout_genera/ds/analyze.py
```

Soft calls are `probability_sum`. Error rates below are base-weighted.

On samovar10, 1.97 Mb of the wrong bases are Enterococcus sp. DIV2432 called as *E. faecalis*, and the neighbour vote matches that wrong call. Nine other nodes (0.85 Mb) have a neighbour vote equal to the truth; every graph smoother that rescues them also harms other nodes. `label_drop` and `logistic_drop` change no argmax.

On the held-out genera, 6423 of 6439 nodes are isolated (24 edges). Wrong bases are unclassified or a same-genus swap, led by Pseudomonas fluorescens. The half-strain graph is the same shape: 10035 of 10072 nodes are isolated. Neighbour rules cannot move those calls.

Five follow-up hypotheses come from that split: lock a confident call, adopt a neighbour only when every neighbour agrees and the call is uncertain, drop edges whose confident endpoints disagree before leakage, fill an uncertain node with the plurality species of its genus, and keep a close same-genus second label instead of collapsing it.
