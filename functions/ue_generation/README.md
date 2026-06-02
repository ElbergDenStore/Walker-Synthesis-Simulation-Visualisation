# functions/ue_generation/

User-equipment (UE) grid generators. Each returns flattened lat/lon vectors.

| File | Purpose |
|---|---|
| `generate_equal_area_ues.m` | Equal-area UE grid (accounts for true spherical area). |
| `generate_population_ues.m` | Population-weighted UE grid (WorldPop). Also returns total population. |
| `generate_random_ues.m` | Uniform-random UE grid. |
