# functions/link_budget/

Vectorised RF link-budget engine.

| File | Purpose |
|---|---|
| `link_calc_matrix.m` | Vectorised link budget: FSPL, P.618 attenuation, antenna gain, SINR, PFD. |
| `beam_gain_and_interference.m` | Per-UE carrier and interference PSDs via KD-tree nearest-beam lookup. |
| `calculate_throughput_matrix.m` | SINR → throughput (modified Shannon, with bandwidth sharing). |
| `PFD_calc.m` | Reverse-solve required Tx power for an ITU power-flux-density target. |
