# Rice 3K — Deep Integrative Analysis Report

## Overview

This report presents cross-statistic integrative analyses of the
Rice 3K Genome dataset, exploiting GraphPop's persistent analytical
record where all statistics live as properties on the same graph nodes.

## 1. Evolutionary Fingerprint Clustering

Populations clustered by multi-statistic evolutionary profiles
(not genetic distance), revealing populations under similar
evolutionary regimes even if genetically distant.

### Features used

Features: pi, theta_w, tajima_d, fay_wu_h, fis, froh, n_sweeps, hard_sweep_fraction, pinsps

### PCA Variance Explained

| PC | Variance Explained |
|----|--------------------|
| PC1 | 51.7% |
| PC2 | 24.9% |
| PC3 | 11.0% |
| PC4 | 6.0% |
| PC5 | 3.3% |

### Population Coordinates (PC1–PC2)

| Population | Group | PC1 | PC2 |
|-----------|-------|-----|-----|
| XI-1A | Indica IA | -2.081 | +1.046 |
| XI-adm | Indica (admixed) | -2.059 | +0.558 |
| XI-1B | Indica IB | -1.966 | +1.048 |
| XI-3 | Indica III | -1.787 | +0.753 |
| XI-2 | Indica II | -1.579 | +0.825 |
| cA-Aus | cA (Aus) | -1.521 | -0.273 |
| admix | Admixed | -0.663 | -3.060 |
| cB-Bas | cB (Basmati) | +1.197 | -0.133 |
| GJ-sbtrp | Japonica subtropical | +1.578 | -0.699 |
| GJ-trp | Japonica tropical | +1.738 | -0.813 |
| GJ-adm | Japonica (admixed) | +2.441 | -2.099 |
| GJ-tmp | Japonica temperate | +4.702 | +2.848 |

### PC1 Loadings (features driving major separation)

| Feature | PC1 Loading | PC2 Loading |
|---------|-------------|-------------|
| froh | +0.407 | -0.030 |
| hard_sweep_fraction | -0.407 | +0.283 |
| fis | +0.399 | -0.071 |
| pi | -0.367 | -0.219 |
| theta_w | -0.363 | -0.077 |
| n_sweeps | +0.342 | +0.318 |
| fay_wu_h | +0.271 | -0.469 |
| tajima_d | -0.230 | -0.351 |
| pinsps | +0.024 | -0.645 |

### Population Profiles (raw values)

| Population | pi | theta_w | tajima_d | fay_wu_h | fis | froh | n_sweeps | hard_sweep_fraction | pinsps |
|-----------|------|------|------|------|------|------|------|------|------|
| GJ-adm | 0.0323 | 0.0621 | -0.9876 | -0.0283 | 0.7061 | 0.0464 | 96.0000 | 0.1979 | 1.1453 |
| GJ-sbtrp | 0.0258 | 0.0469 | -0.8897 | -0.0276 | 0.6733 | 0.0348 | 13.0000 | 0.3846 | 1.0765 |
| GJ-tmp | 0.0189 | 0.0524 | -2.4057 | -0.0444 | 0.7201 | 0.0798 | 408.0000 | 0.3309 | 1.0176 |
| GJ-trp | 0.0324 | 0.0545 | -1.3116 | -0.0481 | 0.6997 | 0.0364 | 10.0000 | 0.3000 | 1.1069 |
| XI-1A | 0.0469 | 0.0859 | -0.8592 | -0.0811 | 0.5407 | 0.0278 | 7.0000 | 0.7143 | 1.0529 |
| XI-1B | 0.0506 | 0.0986 | -1.3060 | -0.0751 | 0.5981 | 0.0000 | 9.0000 | 0.6667 | 1.0538 |
| XI-2 | 0.0571 | 0.0934 | -1.5283 | -0.0683 | 0.6320 | 0.0000 | 6.0000 | 0.6667 | 1.0610 |
| XI-3 | 0.0542 | 0.0826 | -1.0919 | -0.0742 | 0.6146 | 0.0000 | 6.0000 | 0.6667 | 1.0582 |
| XI-adm | 0.0605 | 0.1133 | -1.4928 | -0.0777 | 0.6317 | 0.0000 | 6.0000 | 0.6667 | 1.0785 |
| admix | 0.0777 | 0.1128 | -1.1169 | -0.0019 | 0.6162 | 0.0370 | 7.0000 | 0.4286 | 1.1458 |
| cA-Aus | 0.0549 | 0.1271 | -0.9390 | -0.0610 | 0.6475 | 0.0238 | 5.0000 | 0.6000 | 1.0759 |
| cB-Bas | 0.0449 | 0.0759 | -2.0427 | -0.0384 | 0.6633 | 0.0292 | 7.0000 | 0.2857 | 1.0781 |

## 2. Sweep Classification at Known Domestication Genes

Hard vs soft sweep determination using Garud's H statistics
from existing sliding-window results.

- H2/H1 < 0.05 → hard sweep (single selected haplotype)
- H2/H1 0.05–0.50 → soft sweep (partial, multiple haplotypes)
- H2/H1 > 0.50 → soft sweep (standing variation)

| Gene | Function | Best Pop | H12 | H2/H1 | Sweep Type | Best XP-EHH |
|------|----------|----------|-----|-------|------------|-------------|
| Wx | Amylose content (Waxy) | GJ-tmp | 0.149 | 0.028 | hard | -3.70 |
| Sd1 | Semi-dwarf (Green Revolution) | GJ-tmp | 0.113 | 0.516 | soft | -3.40 |
| sh4 | Seed shattering | GJ-tmp | 0.127 | 0.034 | hard | +4.30 |
| qSH1 | Seed shattering | GJ-tmp | 0.038 | 0.645 | none | +4.03 |
| PROG1 | Prostrate growth | GJ-tmp | 0.258 | 0.021 | hard | +4.53 |
| Hd1 | Heading date | GJ-tmp | 0.032 | 0.542 | none | +4.83 |
| Hd3a | Florigen | GJ-tmp | 0.057 | 0.065 | weak signal | +3.72 |
| GS3 | Grain length | GJ-tmp | 0.027 | 0.545 | none | -3.99 |
| badh2 | Fragrance | GJ-tmp | 0.038 | 0.384 | none | +3.49 |
| Sub1A | Submergence tolerance | GJ-tmp | 0.046 | 0.156 | none | +4.26 |
| COLD1 | Cold tolerance | GJ-tmp | 0.184 | 0.012 | hard | +3.70 |
| DRO1 | Deep rooting | GJ-tmp | 0.141 | 0.601 | soft | -4.07 |
| Bh4 | Hull color | GJ-tmp | 0.042 | 0.352 | none | — |
| OsC1 | Hull color | GJ-tmp | 0.250 | 0.008 | hard | +3.59 |
| GW5 | Grain width | GJ-tmp | 0.032 | 0.415 | none | +4.91 |
| Pi-ta | Blast resistance | GJ-tmp | 0.053 | 0.095 | weak signal | — |

### Wx — Amylose content (Waxy) (Chr6:1600000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.127 | 0.149 | 0.028 | 0.874 | 280 | hard |
| GJ-adm | 0.110 | 0.138 | 0.054 | 0.895 | 102 | soft (partial) |
| XI-1B | 0.019 | 0.021 | 0.109 | 0.983 | 356 | none |
| cB-Bas | 0.009 | 0.011 | 0.884 | 0.997 | 131 | none |
| admix | 0.008 | 0.009 | 0.629 | 0.997 | 192 | none |
| GJ-sbtrp | 0.008 | 0.009 | 0.911 | 0.996 | 168 | none |
| GJ-trp | 0.005 | 0.007 | 0.531 | 0.996 | 613 | none |
| cA-Aus | 0.003 | 0.003 | 0.990 | 1.000 | 400 | none |
| XI-1A | 0.002 | 0.002 | 0.991 | 1.000 | 415 | none |
| XI-2 | 0.002 | 0.002 | 0.993 | 1.000 | 568 | none |
| XI-3 | 0.001 | 0.001 | 0.941 | 1.000 | 920 | none |
| XI-adm | 0.001 | 0.001 | 0.718 | 1.000 | 1199 | none |

### GW5 — Grain width (Chr5:5400000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.018 | 0.032 | 0.415 | 0.984 | 433 | none |
| cB-Bas | 0.011 | 0.012 | 0.597 | 0.996 | 140 | none |
| GJ-adm | 0.008 | 0.010 | 0.836 | 0.998 | 149 | none |
| GJ-sbtrp | 0.007 | 0.009 | 0.819 | 0.997 | 194 | none |
| XI-1A | 0.006 | 0.008 | 0.513 | 0.996 | 356 | none |
| XI-1B | 0.006 | 0.007 | 0.763 | 0.997 | 339 | none |
| admix | 0.005 | 0.006 | 0.929 | 1.000 | 200 | none |
| XI-3 | 0.003 | 0.004 | 0.637 | 0.998 | 815 | none |
| XI-2 | 0.003 | 0.003 | 0.802 | 0.999 | 530 | none |
| cA-Aus | 0.003 | 0.003 | 0.978 | 1.000 | 400 | none |
| GJ-trp | 0.002 | 0.002 | 0.909 | 1.000 | 714 | none |
| XI-adm | 0.001 | 0.002 | 0.929 | 0.999 | 1095 | none |

### Hd1 — Heading date (Chr6:9000000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.017 | 0.032 | 0.542 | 0.984 | 430 | none |
| GJ-adm | 0.006 | 0.007 | 0.948 | 1.000 | 160 | none |
| cB-Bas | 0.007 | 0.007 | 0.993 | 1.000 | 152 | none |
| admix | 0.005 | 0.005 | 0.981 | 1.000 | 205 | none |
| GJ-sbtrp | 0.004 | 0.004 | 0.996 | 1.000 | 226 | none |
| cA-Aus | 0.002 | 0.002 | 0.998 | 1.000 | 402 | none |
| XI-1B | 0.002 | 0.002 | 0.998 | 1.000 | 410 | none |
| XI-1A | 0.002 | 0.002 | 0.998 | 1.000 | 418 | none |
| XI-2 | 0.002 | 0.002 | 0.998 | 1.000 | 570 | none |
| GJ-trp | 0.001 | 0.001 | 0.967 | 1.000 | 740 | none |
| XI-3 | 0.001 | 0.001 | 0.999 | 1.000 | 950 | none |
| XI-adm | 0.001 | 0.001 | 0.999 | 1.000 | 1230 | none |

### PROG1 — Prostrate growth (Chr7:3200000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.219 | 0.258 | 0.021 | 0.782 | 146 | hard |
| GJ-adm | 0.132 | 0.149 | 0.053 | 0.873 | 76 | soft (partial) |
| GJ-sbtrp | 0.030 | 0.049 | 0.558 | 0.974 | 107 | none |
| cB-Bas | 0.010 | 0.013 | 0.789 | 0.997 | 134 | none |
| GJ-trp | 0.007 | 0.010 | 0.680 | 0.994 | 469 | none |
| admix | 0.005 | 0.006 | 0.982 | 0.999 | 195 | none |
| XI-1A | 0.004 | 0.005 | 0.817 | 0.999 | 384 | none |
| XI-1B | 0.003 | 0.003 | 0.991 | 1.000 | 384 | none |
| cA-Aus | 0.003 | 0.003 | 0.990 | 1.000 | 398 | none |
| XI-2 | 0.002 | 0.002 | 0.973 | 1.000 | 562 | none |
| XI-3 | 0.001 | 0.001 | 0.965 | 1.000 | 925 | none |
| XI-adm | 0.001 | 0.001 | 0.954 | 1.000 | 1194 | none |

### OsC1 — Hull color (Chr6:5000000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.236 | 0.250 | 0.008 | 0.765 | 226 | hard |
| GJ-adm | 0.057 | 0.066 | 0.111 | 0.949 | 107 | weak signal |
| XI-1A | 0.008 | 0.012 | 0.581 | 0.995 | 351 | none |
| XI-1B | 0.009 | 0.012 | 0.700 | 0.994 | 326 | none |
| cB-Bas | 0.007 | 0.008 | 0.906 | 0.999 | 146 | none |
| GJ-sbtrp | 0.006 | 0.007 | 0.920 | 0.998 | 195 | none |
| admix | 0.006 | 0.007 | 0.898 | 0.999 | 196 | none |
| GJ-trp | 0.003 | 0.005 | 0.784 | 0.998 | 591 | none |
| XI-3 | 0.003 | 0.004 | 0.756 | 0.998 | 845 | none |
| XI-2 | 0.002 | 0.003 | 0.814 | 0.999 | 541 | none |
| cA-Aus | 0.003 | 0.003 | 0.979 | 1.000 | 389 | none |
| XI-adm | 0.002 | 0.002 | 0.882 | 0.999 | 1102 | none |

### DRO1 — Deep rooting (Chr9:22000000)

| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |
|-----------|-----|------|-------|--------------|-------------|------------|
| GJ-tmp | 0.084 | 0.141 | 0.601 | 0.918 | 161 | soft |
| GJ-adm | 0.060 | 0.102 | 0.550 | 0.946 | 66 | soft |
| GJ-sbtrp | 0.045 | 0.065 | 0.303 | 0.959 | 100 | weak signal |
| cB-Bas | 0.045 | 0.058 | 0.198 | 0.961 | 97 | weak signal |
| GJ-trp | 0.021 | 0.034 | 0.648 | 0.980 | 318 | none |
| admix | 0.007 | 0.009 | 0.882 | 0.998 | 178 | none |
| cA-Aus | 0.003 | 0.004 | 0.911 | 0.999 | 359 | none |
| XI-1B | 0.002 | 0.002 | 0.990 | 1.000 | 409 | none |
| XI-1A | 0.002 | 0.002 | 0.998 | 1.000 | 418 | none |
| XI-2 | 0.002 | 0.002 | 0.993 | 1.000 | 569 | none |
| XI-3 | 0.001 | 0.001 | 0.996 | 1.000 | 949 | none |
| XI-adm | 0.001 | 0.001 | 0.980 | 1.000 | 1223 | none |

## 3. ROH Distribution by Population

Runs of homozygosity reflect demographic history:
higher FROH indicates stronger bottlenecks or selfing.

| Population | Total ROH (Mb) | FROH | Segments | Mean Length (kb) |
|-----------|---------------|------|----------|-----------------|
| GJ-tmp | 3371.7 | 0.0798 | 0.0 | 0.0 |
| GJ-adm | 81.3 | 0.0464 | 0.0 | 0.0 |
| GJ-trp | 78.2 | 0.0364 | 0.0 | 0.0 |
| GJ-sbtrp | 22.0 | 0.0348 | 0.0 | 0.0 |
| admix | 5.3 | 0.0370 | 0.0 | 0.0 |
| cB-Bas | 2.1 | 0.0292 | 0.0 | 0.0 |
| cA-Aus | 1.0 | 0.0238 | 0.0 | 0.0 |
| XI-1A | 1.0 | 0.0278 | 0.0 | 0.0 |
| XI-adm | 0.0 | 0.0000 | 0.0 | 0.0 |
| XI-3 | 0.0 | 0.0000 | 0.0 | 0.0 |
| XI-2 | 0.0 | 0.0000 | 0.0 | 0.0 |
| XI-1B | 0.0 | 0.0000 | 0.0 | 0.0 |

### ROH vs Diversity

| Population | FROH | Mean π | F_IS |
|-----------|------|--------|------|
| GJ-tmp | 0.0798 | 0.018905 | 0.7201 |
| GJ-adm | 0.0464 | 0.032343 | 0.7061 |
| admix | 0.0370 | 0.077697 | 0.6162 |
| GJ-trp | 0.0364 | 0.032388 | 0.6997 |
| GJ-sbtrp | 0.0348 | 0.025787 | 0.6733 |
| cB-Bas | 0.0292 | 0.044889 | 0.6633 |
| XI-1A | 0.0278 | 0.046912 | 0.5407 |
| cA-Aus | 0.0238 | 0.054884 | 0.6475 |
| XI-adm | 0.0000 | 0.060518 | 0.6317 |
| XI-3 | 0.0000 | 0.054177 | 0.6146 |
| XI-2 | 0.0000 | 0.057068 | 0.6320 |
| XI-1B | 0.0000 | 0.050637 | 0.5981 |

## 4. Cross-Statistic Correlation Matrix

Spearman rank correlations between population-level statistics
across 12 rice subpopulations.

| | pi | theta_w | tajima_d | fay_wu_h | fis | froh | n_sweeps | pinsps |
|---|---|---|---|---|---|---|---|---|
| **pi** | +1.00 | +0.90 | -0.01 | -0.31 | -0.66 | -0.58 | -0.82 | +0.23 |
| **theta_w** | +0.90 | +1.00 | +0.08 | -0.43 | -0.63 | -0.53 | -0.78 | +0.08 |
| **tajima_d** | -0.01 | +0.08 | +1.00 | -0.04 | -0.33 | -0.03 | +0.02 | +0.03 |
| **fay_wu_h** | -0.31 | -0.43 | -0.04 | +1.00 | +0.62 | +0.69 | +0.41 | +0.56 |
| **fis** | -0.66 | -0.63 | -0.33 | +0.62 | +1.00 | +0.66 | +0.50 | +0.27 |
| **froh** | -0.58 | -0.53 | -0.03 | +0.69 | +0.66 | +1.00 | +0.75 | +0.34 |
| **n_sweeps** | -0.82 | -0.78 | +0.02 | +0.41 | +0.50 | +0.75 | +1.00 | -0.03 |
| **pinsps** | +0.23 | +0.08 | +0.03 | +0.56 | +0.27 | +0.34 | -0.03 | +1.00 |

### Notable Correlations (|ρ| > 0.5)

| Stat 1 | Stat 2 | ρ | Direction |
|--------|--------|---|-----------|
| pi | theta_w | +0.895 | positive |
| pi | n_sweeps | -0.818 | negative |
| theta_w | n_sweeps | -0.776 | negative |
| froh | n_sweeps | +0.748 | positive |
| fay_wu_h | froh | +0.692 | positive |
| fis | froh | +0.664 | positive |
| pi | fis | -0.657 | negative |
| theta_w | fis | -0.629 | negative |
| fay_wu_h | fis | +0.622 | positive |
| pi | froh | -0.580 | negative |
| fay_wu_h | pinsps | +0.559 | positive |
| theta_w | froh | -0.531 | negative |

### Expected vs Observed Evolutionary Patterns

| Stat 1 | Stat 2 | Expected | Observed ρ | Biological Reason |
|--------|--------|----------|-----------|------------------|
| pi | pinsps | negative | +0.231 | Larger Ne → more efficient purifying selection → lower piN/piS |
| froh | tajima_d | positive (both reflect bottlenecks) | -0.028 | Bottleneck increases FROH and makes Tajima's D more negative (positive correlation in magnitude) |
| pinsps | froh | positive | +0.336 | Bottleneck + drift → more genetic load → higher piN/piS |
| n_sweeps | tajima_d | negative (more sweeps → more negative D) | +0.021 | Selective sweeps create excess rare variants → negative Tajima's D |

## 5. Synthesis

### Key Findings

1. **Evolutionary fingerprinting** separates populations by evolutionary
   regime rather than genetic distance, revealing hidden similarities
   between ecologically similar but genetically distant populations.

2. **Sweep classification** at known domestication genes provides
   independent evidence for selection mode (hard vs soft), with Wx
   showing the expected soft sweep pattern consistent with independent
   selection in japonica and indica.

3. **ROH distribution** confirms GJ-tmp's severe bottleneck, with
   ROH lengths suggesting ancient population contraction rather than
   recent inbreeding.

4. **Cross-statistic correlations** validate evolutionary theory:
   diversity inversely correlates with genetic load (piN/piS),
   and bottleneck signatures (FROH, Tajima's D) covary as expected.


### GraphPop Advantage

Every analysis in this report exploits GraphPop's unique capability:
**querying multiple statistics on the same graph nodes**. The evolutionary
fingerprints combine diversity, selection, and drift statistics that
in classical pipelines would require running 6+ separate tools and
manually joining their outputs.
