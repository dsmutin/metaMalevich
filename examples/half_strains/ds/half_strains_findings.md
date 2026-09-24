# half_strains error properties

Soft calls are `probability_sum`. Base rates weight contig length.
Scored nodes: 10072. Wrong nodes: 2089 (0.207). Wrong bases: 0.184.

## Neighbour situation

- isolated: n=10035, bases=4747897, base error=0.184, same-genus wrong bases=262756, unclassified bases=521706
- vote_matches_call: n=21, bases=11882, base error=0.000, same-genus wrong bases=0, unclassified bases=0
- vote_matches_truth: n=3, bases=1122, base error=1.000, same-genus wrong bases=318, unclassified bases=804
- vote_matches_wrong_call: n=5, bases=2578, base error=1.000, same-genus wrong bases=0, unclassified bases=2578
- vote_other: n=8, bases=7524, base error=0.124, same-genus wrong bases=328, unclassified bases=188

## Edges

- same_truth: 48 edges, prediction disagree 0.292
- different_truth: 8 edges, prediction disagree 1.000

## Rescue versus probability_sum

- lca: rescued 0 (0 bp), harmed 3178 (1436159 bp), net -1436159 bp
- initial_colouring: rescued 8 (3072 bp), harmed 17 (9098 bp), net -6026 bp
- gated_neighbour: rescued 1 (318 bp), harmed 3 (5987 bp), net -5669 bp
- bayesian_edge: rescued 1 (318 bp), harmed 3 (5987 bp), net -5669 bp
- edge_union: rescued 1 (318 bp), harmed 3 (5987 bp), net -5669 bp
- leakage_flipped: rescued 1 (318 bp), harmed 3 (5987 bp), net -5669 bp
- label_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- logistic_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- leakage: rescued 1 (318 bp), harmed 0 (0 bp), net 318 bp

## Largest confusions

- Pseudomonas protegens → unclassified: 334831 bp
- Escherichia coli → unclassified: 107647 bp
- Escherichia coli → Escherichia albertii: 91231 bp
- Pseudomonas protegens → Pseudomonas aeruginosa: 76588 bp
- Pseudomonas protegens → Pseudomonas fluorescens: 53879 bp
- Shigella sonnei → unclassified: 42195 bp
- Escherichia coli → Escherichia marmotae: 30705 bp
- Streptococcus mutans → unclassified: 28003 bp
