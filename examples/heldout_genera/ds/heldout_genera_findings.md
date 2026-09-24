# heldout_genera error properties

Soft calls are `probability_sum`. Base rates weight contig length.
Scored nodes: 6439. Wrong nodes: 3381 (0.525). Wrong bases: 0.526.

## Neighbour situation

- isolated: n=6423, bases=2614287, base error=0.526, same-genus wrong bases=513939, unclassified bases=649527
- vote_matches_call: n=7, bases=5845, base error=0.000, same-genus wrong bases=0, unclassified bases=0
- vote_matches_truth: n=2, bases=2575, base error=1.000, same-genus wrong bases=2575, unclassified bases=0
- vote_matches_wrong_call: n=3, bases=4721, base error=1.000, same-genus wrong bases=562, unclassified bases=0
- vote_other: n=4, bases=4913, base error=0.128, same-genus wrong bases=628, unclassified bases=0

## Edges

- same_truth: 16 edges, prediction disagree 0.250
- different_truth: 8 edges, prediction disagree 0.750

## Rescue versus probability_sum

- lca: rescued 0 (0 bp), harmed 1676 (706351 bp), net -706351 bp
- initial_colouring: rescued 44 (19345 bp), harmed 52 (23933 bp), net -4588 bp
- bayesian_edge: rescued 2 (2575 bp), harmed 2 (3923 bp), net -1348 bp
- label_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- logistic_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- gated_neighbour: rescued 2 (2575 bp), harmed 1 (362 bp), net 2213 bp
- leakage_flipped: rescued 2 (2575 bp), harmed 1 (362 bp), net 2213 bp
- edge_union: rescued 2 (2575 bp), harmed 0 (0 bp), net 2575 bp
- leakage: rescued 2 (2575 bp), harmed 0 (0 bp), net 2575 bp

## Largest confusions

- Pseudomonas fluorescens → unclassified: 394939 bp
- Pseudomonas fluorescens → Pseudomonas protegens: 188276 bp
- Pseudomonas putida → unclassified: 81061 bp
- Pseudomonas fluorescens → Pseudomonas putida: 72051 bp
- Shigella boydii → unclassified: 67078 bp
- Shigella boydii → Escherichia marmotae: 54442 bp
- Pseudomonas fluorescens → Pseudomonas aeruginosa: 48254 bp
- Shigella boydii → Shigella flexneri: 46555 bp
