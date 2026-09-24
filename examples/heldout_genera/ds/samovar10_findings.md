# samovar10 error properties

Soft calls are `probability_sum`. Base rates weight contig length.
Scored nodes: 399. Wrong nodes: 30 (0.075). Wrong bases: 0.073.

## Neighbour situation

- vote_matches_call: n=350, bases=34706240, base error=0.000, same-genus wrong bases=0, unclassified bases=0
- vote_matches_truth: n=9, bases=853044, base error=1.000, same-genus wrong bases=453044, unclassified bases=400000
- vote_matches_wrong_call: n=20, bases=1973006, base error=1.000, same-genus wrong bases=1973006, unclassified bases=0
- vote_other: n=20, bases=1753848, base error=0.017, same-genus wrong bases=0, unclassified bases=0

## Edges

- same_truth: 3944 edges, prediction disagree 0.059
- different_truth: 506 edges, prediction disagree 0.486

## Rescue versus probability_sum

- lca: rescued 0 (0 bp), harmed 292 (28735540 bp), net -28735540 bp
- bayesian_edge: rescued 3 (300000 bp), harmed 8 (754761 bp), net -454761 bp
- gated_neighbour: rescued 9 (853044 bp), harmed 12 (1154761 bp), net -301717 bp
- initial_colouring: rescued 0 (0 bp), harmed 4 (259872 bp), net -259872 bp
- edge_union: rescued 3 (300000 bp), harmed 6 (554761 bp), net -254761 bp
- leakage: rescued 3 (300000 bp), harmed 6 (554761 bp), net -254761 bp
- cut_then_leak: rescued 3 (300000 bp), harmed 6 (554761 bp), net -254761 bp
- confident_lock: rescued 3 (300000 bp), harmed 4 (400000 bp), net -100000 bp
- leakage_flipped: rescued 9 (853044 bp), harmed 10 (863557 bp), net -10513 bp
- label_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- logistic_drop: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- unanimous_rescue: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- margin_mixture: rescued 0 (0 bp), harmed 0 (0 bp), net 0 bp
- genus_plurality: rescued 4 (400000 bp), harmed 0 (0 bp), net 400000 bp

## Largest confusions

- Enterococcus sp. DIV2432 → Enterococcus faecalis: 1973006 bp
- Enterobacter ludwigii → unclassified: 300000 bp
- Enterococcus sp. DIV0180 → Enterococcus sp. DIV0176: 100000 bp
- Enterobacter ludwigii → Enterobacter asburiae: 100000 bp
- Enterobacter ludwigii → Enterobacter cloacae complex sp. FDA-CDC-AR_0132: 100000 bp
- Enterobacter ludwigii → Enterobacter hormaechei: 100000 bp
- Flagellimonas sp. → unclassified: 100000 bp
- Enterococcus sp. DIV0180 → Enterococcus sp. DIV0210g: 53044 bp
