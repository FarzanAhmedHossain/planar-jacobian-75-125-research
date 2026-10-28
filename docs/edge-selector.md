# Corrected-power selector for the quartic F2 edge

This note records the exact selector identities for the normalized quartic edge carrier used in the `(75,125)` F2 model. These are restrictions on the residual edge polynomial, not an exclusion of the original degree pair.

The script `src/f2_75_125_power_relation_edge_certificate.py` checks the candidate faces by exact rational arithmetic. The paired result is in `results/f2_75_125_power_relation_edge_certificate.json`.

## Exact outcomes

- For the `2+2` root pattern, the edge parameter satisfies `alpha^2 - 3 alpha + 1 = 0`.
- For the `2+1+1` pattern, the parameters satisfy `alpha + beta = 3` and `alpha beta = 3`; equivalently the selected parameter satisfies `t^2 - 3t + 3 = 0`.
- The candidates pass the selector equations. Therefore the selector does not eliminate either pattern.

The result is a necessary-condition certificate inside the stated normalized F2 model. It is not a proof that these parameters extend to a full Keller pair, and it is not a contradiction. Later coefficient layers and a scope-correct bridge to the original pair are still required for a global exclusion.
