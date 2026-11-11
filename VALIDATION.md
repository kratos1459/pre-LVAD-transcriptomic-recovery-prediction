# Rebuild vs manuscript (7 bulk cohorts, 1,108 samples)

| Quantity | Manuscript | Rebuild |
|---|---|---|
| Consensus HF genes (FDR<1e-6, abs g>0.8) | 607 (394 up / 213 down) | 624 (382 / 242) |
| Top genes | FCN3, AQP3, RARRES1, HTRA1, EVL, ARHGAP1, GRB14... | FCN3, AQP3, RARRES1, LUM, GRB14, HMOX2, HMGN2, HTRA1, SCUBE2 |
| LOCO AUC, 20 genes, L2 | 0.981 | 0.968 |
| Batch x outcome (11943 NR / 10931 NR) | 32/37, 15/26 | 32/37, 15/26 (identical) |
| Elastic-net CV AUC (in-fold batch adj.) | 0.75 | 0.72 |
| Permutation p (500) | 0.004 | 0.016 |
| Held-out batch 11943 / 10931 (enet) | 0.67 / 0.75 | 0.70 / 0.82 |
| LRRN4CL / ATF3 Welch p | 0.004 / 0.003 | 0.004 / 0.003 |
| Genes at FDR<0.05 | 1 | 0 |
| Subtypes (K by PAC) | S2 immune-fibrotic, unrelated to recovery | K=2, unrelated to recovery (p=1.0) |
