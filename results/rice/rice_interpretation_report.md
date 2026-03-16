# Rice 3K Genome — Biological Interpretation Report

## 1. Population Diversity Ranking

| Rank | Population | Mean π | Mean θ_W | Mean Tajima's D | Mean F_IS |
|------|-----------|--------|---------|----------------|----------|
| 1 | admix | 0.077697 | 0.112775 | -1.116929 | 0.616243 |
| 2 | XI-adm | 0.060518 | 0.113254 | -1.492844 | 0.631702 |
| 3 | XI-2 | 0.057068 | 0.093375 | -1.528260 | 0.631975 |
| 4 | cA-Aus | 0.054884 | 0.127114 | -0.938992 | 0.647548 |
| 5 | XI-3 | 0.054177 | 0.082616 | -1.091939 | 0.614633 |
| 6 | XI-1B | 0.050637 | 0.098627 | -1.305966 | 0.598142 |
| 7 | XI-1A | 0.046912 | 0.085906 | -0.859209 | 0.540748 |
| 8 | cB-Bas | 0.044889 | 0.075920 | -2.042719 | 0.663294 |
| 9 | GJ-trp | 0.032388 | 0.054537 | -1.311591 | 0.699731 |
| 10 | GJ-adm | 0.032343 | 0.062095 | -0.987619 | 0.706070 |
| 11 | GJ-sbtrp | 0.025787 | 0.046856 | -0.889718 | 0.673329 |
| 12 | GJ-tmp | 0.018905 | 0.052439 | -2.405701 | 0.720141 |

**Key observations:**
- Most diverse: **admix** (π = 0.077697)
- Least diverse: **GJ-tmp** (π = 0.018905)
- GJ-tmp (Temperate Japonica) F_IS = 0.720141 — strongest inbreeding among all populations

## 2. Purifying Selection Efficiency (π_N/π_S)

| Rank | Population | Mean π_N/π_S | Interpretation |
|------|-----------|-------------|----------------|
| 1 | admix | 1.145750 | relaxed constraint |
| 2 | GJ-adm | 1.145290 | relaxed constraint |
| 3 | GJ-trp | 1.106850 | relaxed constraint |
| 4 | XI-adm | 1.078519 | relaxed constraint |
| 5 | cB-Bas | 1.078131 | relaxed constraint |
| 6 | GJ-sbtrp | 1.076512 | relaxed constraint |
| 7 | cA-Aus | 1.075913 | relaxed constraint |
| 8 | XI-2 | 1.060974 | relaxed constraint |
| 9 | XI-3 | 1.058247 | relaxed constraint |
| 10 | XI-1B | 1.053750 | relaxed constraint |
| 11 | XI-1A | 1.052860 | relaxed constraint |
| 12 | GJ-tmp | 1.017585 | relaxed constraint |

**Key finding:** Populations with the highest π_N/π_S have accumulated more slightly deleterious missense variants, consistent with reduced effectiveness of purifying selection during domestication bottlenecks.

### π_N/π_S per chromosome (GJ-tmp)

| Chr | π_N | π_S | π_N/π_S |
|-----|-----|-----|---------|
| Chr1 | 0.014310 | 0.013858 | 1.032588 |
| Chr2 | 0.013899 | 0.014080 | 0.987112 |
| Chr3 | 0.009817 | 0.009970 | 0.984617 |
| Chr4 | 0.017556 | 0.016809 | 1.044440 |
| Chr5 | 0.011490 | 0.011569 | 0.993199 |
| Chr6 | 0.017437 | 0.016739 | 1.041676 |
| Chr7 | 0.019123 | 0.018504 | 1.033429 |
| Chr8 | 0.015023 | 0.014614 | 1.027962 |
| Chr9 | 0.013204 | 0.013062 | 1.010847 |
| Chr10 | 0.027318 | 0.024731 | 1.104617 |
| Chr11 | 0.024615 | 0.026371 | 0.933396 |
| Chr12 | 0.017343 | 0.017050 | 1.017133 |

## 3. Pairwise Population Divergence (Weir & Cockerham Fst)

### Top 10 most differentiated pairs

| Pair | Mean Fst |
|------|---------|
| GJ-tmp_vs_XI-1A | 0.7105 |
| GJ-tmp_vs_XI-1B | 0.6937 |
| GJ-tmp_vs_cA-Aus | 0.6783 |
| GJ-tmp_vs_XI-2 | 0.6633 |
| XI-3_vs_GJ-tmp | 0.6517 |
| XI-1A_vs_GJ-sbtrp | 0.6357 |
| GJ-trp_vs_XI-1A | 0.6220 |
| XI-1B_vs_GJ-sbtrp | 0.6149 |
| GJ-trp_vs_XI-1B | 0.6034 |
| XI-1A_vs_GJ-adm | 0.5960 |

### Most closely related pairs

| Pair | Mean Fst |
|------|---------|
| XI-adm_vs_XI-3 | 0.0606 |
| XI-adm_vs_XI-1B | 0.0622 |
| XI-adm_vs_XI-2 | 0.0673 |
| XI-adm_vs_XI-1A | 0.0955 |
| XI-3_vs_XI-2 | 0.1250 |

## 4. Selection Signals (XP-EHH Peaks)

### cA-Aus_vs_cB-Bas
Total hits: 29661, Peaks: 398

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr8 | 17,839,077 | 5.33 | 492 | LOC_Os08g28930 | HIGH |
| Chr4 | 13,016,188 | -5.29 | 3254 | LOC_Os04g21680 | HIGH |
| Chr6 | 23,141,504 | -5.16 | 336 | LOC_Os06g38890 | HIGH |
| Chr10 | 14,902,355 | -4.76 | 264 | LOC_Os10g28250 | HIGH |
| Chr10 | 14,109,386 | -4.74 | 1271 | LOC_Os10g26780 | HIGH |
| Chr8 | 9,147,005 | 4.73 | 506 | LOC_Os08g14960 | HIGH |
| Chr10 | 2,224,552 | 4.67 | 406 | LOC_Os10g04290 | HIGH |
| Chr2 | 2,952,104 | 4.67 | 251 | LOC_Os02g05870 | HIGH |
| Chr10 | 13,858,980 | -4.67 | 157 | LOC_Os10g26470 | HIGH |
| Chr1 | 3,069,703 | 4.65 | 167 | LOC_Os01g06360 | HIGH |

### XI-1A_vs_cA-Aus
Total hits: 99624, Peaks: 621

**Known domestication genes detected:**

| Gene | Function | XP-EHH | Peak Position |
|------|----------|--------|--------------|
| GW5 | Grain width | 4.91 | Chr5:5,265,576 |

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr10 | 22,067,307 | -7.95 | 2310 | LOC_Os10g40612 | HIGH |
| Chr3 | 1,186,263 | -6.95 | 130 | LOC_Os03g02770 | HIGH |
| Chr2 | 29,067,508 | -6.74 | 960 | LOC_Os02g47400 | HIGH |
| Chr3 | 15,446,637 | -6.64 | 388 | LOC_Os03g26690 | HIGH |
| Chr3 | 12,942,157 | 6.47 | 1349 | LOC_Os03g22190 | HIGH |
| Chr1 | 3,304,143 | -6.44 | 2091 | LOC_Os01g06330 | HIGH |
| Chr8 | 25,261,054 | 6.32 | 1621 | LOC_Os08g39100 | HIGH |
| Chr10 | 2,405,288 | 6.29 | 857 | LOC_Os10g04850 | HIGH |
| Chr9 | 18,574,496 | -6.23 | 335 | LOC_Os09g30280 | HIGH |
| Chr5 | 25,568,441 | -6.23 | 462 | LOC_Os05g43840 | HIGH |

### GJ-tmp_vs_cA-Aus
Total hits: 17994, Peaks: 248

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr1 | 3,304,047 | -6.63 | 5542 | LOC_Os01g06200 | HIGH |
| Chr4 | 17,901,541 | -6.09 | 3229 | LOC_Os04g29520 | HIGH |
| Chr12 | 26,954,457 | -5.32 | 299 | LOC_Os12g43280 | HIGH |
| Chr1 | 2,817,299 | -5.31 | 448 | LOC_Os01g05730 | HIGH |
| Chr10 | 5,615,132 | -5.09 | 29 | LOC_Os10g10149 | HIGH |
| Chr2 | 24,543,944 | -4.39 | 39 | LOC_Os02g40430 | HIGH |
| Chr3 | 4,934,058 | 4.32 | 16 | LOC_Os03g09810 | HIGH |
| Chr5 | 4,853,032 | 4.25 | 14 | LOC_Os05g08694 | HIGH |
| Chr11 | 28,099,083 | -4.23 | 1124 | LOC_Os11g46140 | HIGH |
| Chr9 | 16,972,910 | 4.15 | 13 | LOC_Os09g27810 | HIGH |

### GJ-tmp_vs_XI-1A
Total hits: 15242, Peaks: 258

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr1 | 2,817,299 | -5.48 | 800 | LOC_Os01g05660 | HIGH |
| Chr4 | 17,901,540 | -5.27 | 4126 | LOC_Os04g29500 | HIGH |
| Chr10 | 10,862,805 | -4.93 | 28 | LOC_Os10g21212 | HIGH |
| Chr9 | 2,841,408 | -4.65 | 420 | LOC_Os09g05040 | HIGH |
| Chr1 | 3,257,212 | -4.38 | 797 | LOC_Os01g06360 | HIGH |
| Chr12 | 25,499,482 | -4.18 | 113 | LOC_Os12g41080 | HIGH |
| Chr3 | 15,006,400 | -4.17 | 229 | LOC_Os03g26100 | HIGH |
| Chr11 | 28,256,439 | -4.13 | 698 | LOC_Os11g46930 | HIGH |
| Chr2 | 24,387,169 | -4.12 | 1502 | LOC_Os02g40010 | HIGH |
| Chr12 | 26,954,457 | -4.06 | 135 | LOC_Os12g43460 | HIGH |

### GJ-trp_vs_XI-1A
Total hits: 18477, Peaks: 348

**Known domestication genes detected:**

| Gene | Function | XP-EHH | Peak Position |
|------|----------|--------|--------------|
| Hd1 | Heading date | -4.15 | Chr6:8,954,814 |
| GW5 | Grain width | -3.56 | Chr5:5,345,355 |

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr7 | 21,401,358 | -5.48 | 666 | LOC_Os07g35610 | HIGH |
| Chr1 | 1,639,506 | -4.91 | 108 | LOC_Os01g03720 | HIGH |
| Chr6 | 9,567,297 | -4.84 | 1426 | LOC_Os06g15990 | HIGH |
| Chr7 | 24,487,633 | -4.81 | 264 | LOC_Os07g40750 | HIGH |
| Chr4 | 17,766,683 | -4.69 | 567 | LOC_Os04g29700 | HIGH |
| Chr9 | 4,944,308 | -4.52 | 690 | LOC_Os09g09130 | HIGH |
| Chr4 | 4,419,067 | -4.49 | 185 | LOC_Os04g08080 | HIGH |
| Chr11 | 25,866,120 | -4.41 | 838 | LOC_Os11g42710 | HIGH |
| Chr4 | 31,127,867 | 4.33 | 145 | LOC_Os04g52190 | HIGH |
| Chr11 | 24,477,028 | -4.25 | 233 | LOC_Os11g40820 | HIGH |

### GJ-tmp_vs_GJ-trp
Total hits: 42018, Peaks: 339

**Known domestication genes detected:**

| Gene | Function | XP-EHH | Peak Position |
|------|----------|--------|--------------|
| PROG1 | Prostrate growth | 4.53 | Chr7:3,178,043 |
| Hd1 | Heading date | 4.28 | Chr6:9,255,657 |

**Top 10 peaks (by |XP-EHH|):**

| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |
|-----|----------|--------|-------|----------|--------|
| Chr2 | 7,223,385 | 5.51 | 1362 | LOC_Os02g13350 | HIGH |
| Chr1 | 29,543,706 | 5.50 | 2132 | LOC_Os01g50970 | HIGH |
| Chr6 | 4,377,175 | 5.42 | 2028 | LOC_Os06g08500 | HIGH |
| Chr7 | 27,504,860 | 5.34 | 671 | LOC_Os07g46010 | HIGH |
| Chr6 | 30,475,218 | 5.33 | 206 | LOC_Os06g50260 | HIGH |
| Chr4 | 17,901,182 | -5.32 | 453 | LOC_Os04g29800 | HIGH |
| Chr6 | 30,259,791 | 5.27 | 410 | LOC_Os06g49480 | HIGH |
| Chr3 | 14,419,951 | 5.24 | 528 | LOC_Os03g25110 | HIGH |
| Chr7 | 20,066,571 | 5.03 | 1158 | LOC_Os07g33320 | HIGH |
| Chr7 | 3,928,923 | 4.96 | 384 | LOC_Os07g07640 | HIGH |

## 5. Selective Sweep Regions (Garud's H)

| Population | Total Sweeps (H12>0.10) | Hard Sweeps | Soft Sweeps |
|-----------|------------------------|-------------|-------------|
| XI-adm | 6 | 4 | 2 |
| XI-3 | 6 | 4 | 2 |
| GJ-trp | 10 | 3 | 7 |
| GJ-tmp | 408 | 135 | 273 |
| XI-2 | 6 | 4 | 2 |
| XI-1A | 7 | 5 | 2 |
| XI-1B | 9 | 6 | 3 |
| cA-Aus | 5 | 3 | 2 |
| GJ-sbtrp | 13 | 5 | 8 |
| admix | 7 | 3 | 4 |
| GJ-adm | 96 | 19 | 77 |
| cB-Bas | 7 | 2 | 5 |

### Notable sweeps near known domestication genes

- **GJ-tmp** Chr6:4,951,023–5,051,022 (H12=0.250, hard): near **OsC1** (Hull color)
- **GJ-tmp** Chr6:5,001,023–5,101,022 (H12=0.163, soft): near **OsC1** (Hull color)
- **GJ-tmp** Chr9:22,037,521–22,137,520 (H12=0.139, soft): near **DRO1** (Deep rooting)
- **GJ-adm** Chr6:1,451,023–1,551,022 (H12=0.138, soft): near **Wx** (Amylose content (Waxy))
- **GJ-adm** Chr6:1,501,023–1,601,022 (H12=0.109, soft): near **Wx** (Amylose content (Waxy))

## 6. Runs of Homozygosity (ROH)

| Rank | Population | Total ROH (Mb) | Mean FROH |
|------|-----------|---------------|-----------|
| 1 | GJ-tmp | 792.8 | 0.056399 |
| 2 | GJ-trp | 4.0 | 0.036199 |
| 3 | GJ-sbtrp | 3.4 | 0.031080 |
| 4 | GJ-adm | 2.1 | 0.028738 |
| 5 | XI-adm | 0.0 | 0.000000 |
| 6 | XI-3 | 0.0 | 0.000000 |
| 7 | XI-2 | 0.0 | 0.000000 |
| 8 | XI-1A | 0.0 | 0.000000 |
| 9 | XI-1B | 0.0 | 0.000000 |
| 10 | cA-Aus | 0.0 | 0.000000 |
| 11 | admix | 0.0 | 0.000000 |
| 12 | cB-Bas | 0.0 | 0.000000 |

## 7. Annotation-Conditioned Genome Scans

### GJ-tmp vs XI-1A: Missense vs Synonymous Fst

- Mean Fst (missense): 0.604204
- Mean Fst (synonymous): 0.578916
- Fst_missense / Fst_synonymous = 1.044
  → **Elevated missense divergence** — evidence of adaptive protein evolution

### GJ-tmp vs GJ-trp: Missense vs Synonymous Fst

- Mean Fst (missense): 0.252016
- Mean Fst (synonymous): 0.227942
- Fst_missense / Fst_synonymous = 1.106
  → **Elevated missense divergence** — evidence of adaptive protein evolution

### XI-1A vs cA-Aus: Missense vs Synonymous Fst

- Mean Fst (missense): 0.366179
- Mean Fst (synonymous): 0.356105
- Fst_missense / Fst_synonymous = 1.028
  → **Elevated missense divergence** — evidence of adaptive protein evolution

### Top Adaptive Divergence Windows (missense-Fst outliers)

| Pair | Chr | Start | End | Fst (Hudson) | #Variants |
|------|-----|-------|-----|-------------|-----------|
| XI-1A_vs_cA-Aus | Chr1 | 7,453,527 | 7,553,526 | 0.9658 | 449 |
| XI-1A_vs_cA-Aus | Chr1 | 7,403,527 | 7,503,526 | 0.9590 | 608 |
| GJ-tmp_vs_XI-1A | Chr5 | 25,004,007 | 25,104,006 | 0.9577 | 356 |
| GJ-tmp_vs_XI-1A | Chr4 | 35,452,105 | 35,552,104 | 0.9548 | 57 |
| GJ-tmp_vs_XI-1A | Chr2 | 12,003,895 | 12,103,894 | 0.9539 | 770 |
| GJ-tmp_vs_XI-1A | Chr5 | 23,404,007 | 23,504,006 | 0.9537 | 604 |
| GJ-tmp_vs_XI-1A | Chr2 | 33,403,895 | 33,503,894 | 0.9467 | 536 |
| GJ-tmp_vs_XI-1A | Chr2 | 16,153,895 | 16,253,894 | 0.9463 | 219 |
| GJ-tmp_vs_XI-1A | Chr5 | 26,804,007 | 26,904,006 | 0.9452 | 479 |
| GJ-tmp_vs_XI-1A | Chr8 | 22,404,150 | 22,504,149 | 0.9446 | 332 |
| GJ-tmp_vs_XI-1A | Chr5 | 25,604,007 | 25,704,006 | 0.9439 | 434 |
| GJ-tmp_vs_XI-1A | Chr2 | 31,603,895 | 31,703,894 | 0.9423 | 316 |
| GJ-tmp_vs_XI-1A | Chr3 | 31,104,689 | 31,204,688 | 0.9413 | 264 |
| GJ-tmp_vs_XI-1A | Chr4 | 29,302,105 | 29,402,104 | 0.9387 | 233 |
| GJ-tmp_vs_XI-1A | Chr6 | 2,251,486 | 2,351,485 | 0.9387 | 539 |
| GJ-tmp_vs_XI-1A | Chr2 | 27,153,895 | 27,253,894 | 0.9374 | 325 |
| GJ-tmp_vs_XI-1A | Chr7 | 4,302,365 | 4,402,364 | 0.9372 | 323 |
| XI-1A_vs_cA-Aus | Chr4 | 34,452,105 | 34,552,104 | 0.9363 | 480 |
| GJ-tmp_vs_XI-1A | Chr7 | 4,352,365 | 4,452,364 | 0.9358 | 500 |
| GJ-tmp_vs_XI-1A | Chr2 | 31,403,895 | 31,503,894 | 0.9355 | 403 |

## 8. Fay & Wu's H (Ancestral-Allele-Dependent)

| Rank | Population | Mean H | Mean H_norm | Interpretation |
|------|-----------|--------|------------|----------------|
| 1 | XI-1A | -0.081149 | -1.139329 | strong sweep signals |
| 2 | XI-adm | -0.077659 | -1.505611 | strong sweep signals |
| 3 | XI-1B | -0.075115 | -1.083278 | strong sweep signals |
| 4 | XI-3 | -0.074208 | -1.406916 | strong sweep signals |
| 5 | XI-2 | -0.068308 | -1.099352 | strong sweep signals |
| 6 | cA-Aus | -0.061046 | -0.810044 | strong sweep signals |
| 7 | GJ-trp | -0.048078 | -0.952752 | moderate sweeps |
| 8 | GJ-tmp | -0.044371 | -0.868007 | moderate sweeps |
| 9 | cB-Bas | -0.038444 | -0.416341 | moderate sweeps |
| 10 | GJ-adm | -0.028311 | -0.456813 | moderate sweeps |
| 11 | GJ-sbtrp | -0.027553 | -0.430037 | moderate sweeps |
| 12 | admix | -0.001883 | -0.043575 | weak/no sweeps |

**Key finding:** Negative H indicates excess high-frequency derived alleles from recent selective sweeps. More negative = stronger/more frequent sweeps in that population.

## 9. Unfolded (Polarized) Site Frequency Spectrum

| Population | Total Variants | Polarized | Singletons | High-freq Derived (>50%) |
|-----------|---------------|-----------|------------|------------------------|
| XI-adm | 29,629,320 | 16,650,808 | 11,987,722 | 1,880,240 |
| XI-3 | 29,594,555 | 16,643,916 | 17,265,546 | 1,976,144 |
| GJ-trp | 29,625,440 | 16,649,776 | 19,952,358 | 1,183,032 |
| GJ-tmp | 29,629,003 | 16,650,727 | 20,986,070 | 889,429 |
| XI-2 | 29,589,136 | 16,643,883 | 17,269,482 | 1,978,648 |
| XI-1A | 29,566,291 | 16,639,041 | 20,363,376 | 2,030,250 |
| XI-1B | 29,562,475 | 16,638,707 | 19,439,593 | 1,965,165 |
| cA-Aus | 29,505,350 | 16,620,990 | 18,698,580 | 1,939,737 |
| GJ-sbtrp | 29,553,290 | 16,631,351 | 24,000,585 | 1,138,393 |
| admix | 29,607,934 | 16,648,113 | 15,063,377 | 1,098,872 |
| GJ-adm | 29,602,931 | 16,646,832 | 22,584,337 | 955,705 |
| cB-Bas | 29,470,420 | 16,616,730 | 21,357,324 | 1,429,875 |

## 10. Fay & Wu's H Genome Scan

### GJ-tmp
- Windows: 7468, Mean H: -0.049089

**Top 10 most negative H windows (strongest sweeps):**

| Chr | Start | End | H | Top Gene | Known Gene |
|-----|-------|-----|---|----------|------------|
| Chr8 | 17,701,028 | 17,801,027 | -0.558003 | LOC_Os08g29030 | — |
| Chr4 | 18,401,004 | 18,501,003 | -0.527693 | LOC_Os04g30900 | — |
| Chr9 | 8,287,521 | 8,387,520 | -0.502085 | LOC_Os09g14120 | — |
| Chr12 | 7,801,048 | 7,901,047 | -0.447286 | LOC_Os12g13880 | — |
| Chr4 | 18,451,004 | 18,551,003 | -0.434444 | LOC_Os04g30900 | — |
| Chr12 | 7,751,048 | 7,851,047 | -0.426839 | LOC_Os12g13770 | — |
| Chr4 | 12,851,004 | 12,951,003 | -0.426026 | LOC_Os04g22680 | — |
| Chr8 | 17,751,028 | 17,851,027 | -0.422079 | LOC_Os08g29080 | — |
| Chr9 | 7,537,521 | 7,637,520 | -0.419277 | LOC_Os09g13230 | — |
| Chr2 | 7,850,030 | 7,950,029 | -0.409022 | LOC_Os02g14410 | — |

### XI-1A
- Windows: 7468, Mean H: -0.082082

**Top 10 most negative H windows (strongest sweeps):**

| Chr | Start | End | H | Top Gene | Known Gene |
|-----|-------|-----|---|----------|------------|
| Chr7 | 14,151,017 | 14,251,016 | -1.313827 | LOC_Os07g24910 | — |
| Chr12 | 3,651,048 | 3,751,047 | -0.644311 | LOC_Os12g07530 | — |
| Chr12 | 5,301,048 | 5,401,047 | -0.622374 | LOC_Os12g10110 | — |
| Chr12 | 5,351,048 | 5,451,047 | -0.622374 | LOC_Os12g10310 | — |
| Chr2 | 18,200,030 | 18,300,029 | -0.575234 | LOC_Os02g30580 | — |
| Chr6 | 1,851,023 | 1,951,022 | -0.531432 | LOC_Os06g04470 | — |
| Chr2 | 8,300,030 | 8,400,029 | -0.513144 | LOC_Os02g14900 | — |
| Chr1 | 3,801,026 | 3,901,025 | -0.474233 | LOC_Os01g07870 | — |
| Chr2 | 18,150,030 | 18,250,029 | -0.469447 | LOC_Os02g30510 | — |
| Chr1 | 3,751,026 | 3,851,025 | -0.462142 | LOC_Os01g07840 | — |

### cA-Aus
- Windows: 7468, Mean H: -0.062099

**Top 10 most negative H windows (strongest sweeps):**

| Chr | Start | End | H | Top Gene | Known Gene |
|-----|-------|-----|---|----------|------------|
| Chr7 | 14,151,017 | 14,251,016 | -1.715984 | LOC_Os07g24910 | — |
| Chr8 | 17,701,028 | 17,801,027 | -0.709575 | LOC_Os08g29030 | — |
| Chr2 | 8,300,030 | 8,400,029 | -0.542865 | LOC_Os02g14900 | — |
| Chr6 | 9,201,023 | 9,301,022 | -0.495390 | LOC_Os06g16320 | — |
| Chr8 | 17,751,028 | 17,851,027 | -0.482474 | LOC_Os08g29080 | — |
| Chr3 | 14,601,047 | 14,701,046 | -0.468071 | LOC_Os03g25590 | — |
| Chr1 | 8,201,026 | 8,301,025 | -0.462211 | LOC_Os01g14770 | — |
| Chr2 | 8,050,030 | 8,150,029 | -0.451608 | LOC_Os02g14590 | — |
| Chr7 | 14,201,017 | 14,301,016 | -0.430462 | LOC_Os07g25010 | — |
| Chr7 | 14,251,017 | 14,351,016 | -0.430462 | LOC_Os07g25080 | — |

### GJ-trp
- Windows: 7468, Mean H: -0.051668

**Top 10 most negative H windows (strongest sweeps):**

| Chr | Start | End | H | Top Gene | Known Gene |
|-----|-------|-----|---|----------|------------|
| Chr2 | 18,200,030 | 18,300,029 | -0.615755 | LOC_Os02g30580 | — |
| Chr6 | 9,201,023 | 9,301,022 | -0.492104 | LOC_Os06g16320 | — |
| Chr4 | 18,401,004 | 18,501,003 | -0.460865 | LOC_Os04g30900 | — |
| Chr4 | 12,851,004 | 12,951,003 | -0.450729 | LOC_Os04g22680 | — |
| Chr4 | 17,651,004 | 17,751,003 | -0.448443 | LOC_Os04g29700 | — |
| Chr12 | 7,801,048 | 7,901,047 | -0.446301 | LOC_Os12g13880 | — |
| Chr12 | 7,751,048 | 7,851,047 | -0.417349 | LOC_Os12g13770 | — |
| Chr4 | 14,351,004 | 14,451,003 | -0.414158 | LOC_Os04g25040 | — |
| Chr7 | 24,401,017 | 24,501,016 | -0.412566 | LOC_Os07g40760 | — |
| Chr7 | 24,451,017 | 24,551,016 | -0.412206 | LOC_Os07g40810 | — |

## 11. Runs of Homozygosity (HMM Method)

| Rank | Population | Total ROH (Mb) | Mean FROH |
|------|-----------|---------------|-----------|
| 1 | GJ-tmp | 3371.7 | 0.079827 |
| 2 | GJ-adm | 81.3 | 0.046371 |
| 3 | admix | 5.3 | 0.037048 |
| 4 | GJ-trp | 78.2 | 0.036390 |
| 5 | GJ-sbtrp | 22.0 | 0.034822 |
| 6 | cB-Bas | 2.1 | 0.029161 |
| 7 | XI-1A | 1.0 | 0.027782 |
| 8 | cA-Aus | 1.0 | 0.023799 |
| 9 | XI-adm | 0.0 | 0.000000 |
| 10 | XI-3 | 0.0 | 0.000000 |
| 11 | XI-2 | 0.0 | 0.000000 |
| 12 | XI-1B | 0.0 | 0.000000 |

**Key finding:** GJ-tmp (Temperate Japonica) has the highest FROH, consistent with the strongest domestication bottleneck among rice subpopulations.

### GJ-tmp — Top chromosomes by ROH

| Chr | Total ROH (Mb) | Mean FROH | Samples with ROH |
|-----|---------------|-----------|------------------|
| Chr1 | 575.4 | 0.086910 | 153/287 |
| Chr2 | 558.8 | 0.088862 | 175/287 |
| Chr3 | 535.3 | 0.096712 | 152/287 |
| Chr4 | 333.6 | 0.070124 | 134/287 |
| Chr5 | 254.9 | 0.078775 | 108/287 |

### GJ-adm — Top chromosomes by ROH

| Chr | Total ROH (Mb) | Mean FROH | Samples with ROH |
|-----|---------------|-----------|------------------|
| Chr3 | 23.3 | 0.039950 | 16/82 |
| Chr2 | 16.6 | 0.046095 | 10/82 |
| Chr4 | 12.6 | 0.039280 | 9/82 |
| Chr1 | 8.1 | 0.046627 | 4/82 |
| Chr9 | 4.6 | 0.049802 | 4/82 |

### GJ-trp — Top chromosomes by ROH

| Chr | Total ROH (Mb) | Mean FROH | Samples with ROH |
|-----|---------------|-----------|------------------|
| Chr3 | 24.8 | 0.029565 | 23/372 |
| Chr9 | 23.7 | 0.051595 | 20/372 |
| Chr2 | 15.3 | 0.032776 | 13/372 |
| Chr6 | 3.5 | 0.036938 | 3/372 |
| Chr10 | 3.2 | 0.045571 | 3/372 |

## 12. Derived Allele Frequency Enrichment at Selection Signals

### PBS Peak Regions

| Focal Population | Peak DAF | Background DAF | Enrichment |
|-----------------|----------|---------------|------------|
| GJ-tmp (GJ-tmp_vs_GJ-trp_out_XI-1A) | 0.062147 | 0.058757 | 1.06x |
| GJ-trp (GJ-trp_vs_GJ-tmp_out_XI-1A) | 0.135927 | 0.070493 | 1.93x |
| XI-1A (XI-1A_vs_XI-1B_out_GJ-tmp) | 0.138682 | 0.101082 | 1.37x |
| cA-Aus (cA-Aus_vs_XI-1A_out_GJ-tmp) | 0.139562 | 0.101309 | 1.38x |
| cB-Bas (cB-Bas_vs_cA-Aus_out_GJ-tmp) | 0.157338 | 0.081528 | 1.93x |

### XP-EHH Peak Regions

| Pair | Focal Pop | Peak DAF | Background DAF | Enrichment |
|------|----------|----------|---------------|------------|
| cA-Aus_vs_cB-Bas | cA-Aus | 0.125541 | 0.101309 | 1.24x |
| XI-1A_vs_cA-Aus | XI-1A | 0.113797 | 0.101082 | 1.13x |
| GJ-tmp_vs_cA-Aus | GJ-tmp | 0.096505 | 0.058757 | 1.64x |
| GJ-tmp_vs_XI-1A | GJ-tmp | 0.094457 | 0.058757 | 1.61x |
| GJ-trp_vs_XI-1A | GJ-trp | 0.119503 | 0.070493 | 1.70x |
| GJ-tmp_vs_GJ-trp | GJ-tmp | 0.091333 | 0.058757 | 1.55x |

**Key finding:** Elevated derived allele frequency at selection signal peaks confirms that sweeps drove derived alleles toward fixation, validating the ancestral allele polarization.

## 13. Population Branch Statistic (PBS)

### GJ-tmp (vs GJ-trp, outgroup XI-1A)
- Windows analysed: 7469
- Mean PBS: 0.3047

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr5 | 28,200,019 | 28,300,018 | 2.7602 | 0.9357 | 4376 | LOC_Os05g49160 | — |
| Chr5 | 28,250,019 | 28,350,018 | 2.5350 | 0.9315 | 4388 | LOC_Os05g49240 | — |
| Chr3 | 11,551,047 | 11,651,046 | 2.5040 | 0.9230 | 4716 | LOC_Os03g20420 | — |
| Chr5 | 28,300,019 | 28,400,018 | 2.3122 | 0.9199 | 3981 | LOC_Os05g49330 | — |
| Chr5 | 28,150,019 | 28,250,018 | 2.3019 | 0.9117 | 5832 | LOC_Os05g49100 | — |
| Chr3 | 11,501,047 | 11,601,046 | 2.2650 | 0.9134 | 5457 | LOC_Os03g20340 | — |
| Chr3 | 11,601,047 | 11,701,046 | 2.2174 | 0.9196 | 4362 | LOC_Os03g20470 | — |
| Chr5 | 28,350,019 | 28,450,018 | 2.1316 | 0.9009 | 3739 | LOC_Os05g49470 | — |
| Chr5 | 28,100,019 | 28,200,018 | 2.0968 | 0.9122 | 5931 | LOC_Os05g49010 | — |
| Chr5 | 27,700,019 | 27,800,018 | 2.0816 | 0.8207 | 3756 | LOC_Os05g48350 | — |

### GJ-trp (vs GJ-tmp, outgroup XI-1A)
- Windows analysed: 7468
- Mean PBS: 0.0511

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr7 | 24,401,017 | 24,501,016 | 1.3687 | 0.6380 | 6203 | LOC_Os07g40750 | — |
| Chr7 | 24,451,017 | 24,551,016 | 1.3327 | 0.6055 | 6813 | LOC_Os07g40820 | — |
| Chr6 | 9,251,023 | 9,351,022 | 1.1647 | 0.6259 | 6742 | LOC_Os06g16260 | — |
| Chr7 | 24,501,017 | 24,601,016 | 1.1419 | 0.6122 | 6809 | LOC_Os07g40940 | — |
| Chr9 | 14,387,521 | 14,487,520 | 1.1417 | 0.7924 | 7266 | LOC_Os09g24230 | — |
| Chr5 | 20,400,019 | 20,500,018 | 1.1403 | 0.8454 | 6285 | LOC_Os05g34440 | — |
| Chr6 | 9,201,023 | 9,301,022 | 1.1034 | 0.6399 | 5976 | LOC_Os06g16190 | — |
| Chr5 | 19,700,019 | 19,800,018 | 1.1030 | 0.8498 | 6897 | LOC_Os05g33554 | — |
| Chr6 | 9,501,023 | 9,601,022 | 1.0919 | 0.5826 | 7725 | LOC_Os06g16570 | — |
| Chr7 | 24,551,017 | 24,651,016 | 1.0824 | 0.6759 | 6797 | LOC_Os07g41080 | — |

### XI-1A (vs XI-1B, outgroup GJ-tmp)
- Windows analysed: 7469
- Mean PBS: 0.1215

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr4 | 24,151,004 | 24,251,003 | 1.7302 | 0.7269 | 5915 | LOC_Os04g40690 | — |
| Chr4 | 24,201,004 | 24,301,003 | 1.7166 | 0.7230 | 5941 | LOC_Os04g40790 | — |
| Chr4 | 24,251,004 | 24,351,003 | 1.6136 | 0.7021 | 5428 | LOC_Os04g40870 | — |
| Chr4 | 24,101,004 | 24,201,003 | 1.5582 | 0.6932 | 5803 | LOC_Os04g40630 | — |
| Chr10 | 21,951,070 | 22,051,069 | 1.3940 | 0.6899 | 5195 | LOC_Os10g40840 | — |
| Chr10 | 21,901,070 | 22,001,069 | 1.3540 | 0.6649 | 5046 | LOC_Os10g40780 | — |
| Chr10 | 21,851,070 | 21,951,069 | 1.3527 | 0.6725 | 4232 | LOC_Os10g40740 | — |
| Chr7 | 1,901,017 | 2,001,016 | 1.3188 | 0.5662 | 6478 | LOC_Os07g04330 | — |
| Chr6 | 2,301,023 | 2,401,022 | 1.2762 | 0.3927 | 6343 | LOC_Os06g05200 | — |
| Chr6 | 2,351,023 | 2,451,022 | 1.2602 | 0.4532 | 5932 | LOC_Os06g05240 | — |

### cA-Aus (vs XI-1A, outgroup GJ-tmp)
- Windows analysed: 7468
- Mean PBS: 0.2132

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr1 | 8,301,026 | 8,401,025 | 2.7553 | 0.9431 | 5440 | LOC_Os01g14870 | — |
| Chr1 | 8,651,026 | 8,751,025 | 2.7010 | 0.9361 | 5466 | LOC_Os01g15460 | — |
| Chr1 | 8,351,026 | 8,451,025 | 2.3165 | 0.9138 | 4986 | LOC_Os01g14910 | — |
| Chr3 | 28,651,047 | 28,751,046 | 2.1231 | 0.8904 | 5165 | LOC_Os03g50290 | — |
| Chr3 | 28,801,047 | 28,901,046 | 2.0982 | 0.8720 | 5987 | LOC_Os03g50490 | — |
| Chr3 | 28,701,047 | 28,801,046 | 2.0569 | 0.8810 | 5432 | LOC_Os03g50325 | — |
| Chr5 | 4,800,019 | 4,900,018 | 2.0500 | 0.8155 | 7498 | LOC_Os05g08780 | — |
| Chr1 | 8,701,026 | 8,801,025 | 2.0244 | 0.8432 | 5087 | LOC_Os01g15490 | — |
| Chr1 | 8,151,026 | 8,251,025 | 2.0203 | 0.8720 | 6543 | LOC_Os01g14540 | — |
| Chr1 | 8,251,026 | 8,351,025 | 2.0054 | 0.8528 | 5867 | LOC_Os01g14780 | — |

### cB-Bas (vs cA-Aus, outgroup GJ-tmp)
- Windows analysed: 7468
- Mean PBS: 0.0979

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr6 | 6,151,023 | 6,251,022 | 1.5904 | 0.8615 | 6414 | LOC_Os06g11610 | — |
| Chr6 | 6,101,023 | 6,201,022 | 1.5567 | 0.8481 | 7523 | LOC_Os06g11510 | — |
| Chr6 | 5,601,023 | 5,701,022 | 1.4605 | 0.5693 | 6836 | LOC_Os06g10750 | — |
| Chr2 | 2,000,030 | 2,100,029 | 1.3564 | 0.6487 | 5884 | LOC_Os02g04510 | — |
| Chr6 | 5,951,023 | 6,051,022 | 1.3426 | 0.7970 | 6056 | LOC_Os06g11360 | — |
| Chr6 | 6,051,023 | 6,151,022 | 1.3361 | 0.8178 | 8524 | LOC_Os06g11470 | — |
| Chr6 | 6,001,023 | 6,101,022 | 1.3343 | 0.8062 | 7634 | LOC_Os06g11390 | — |
| Chr6 | 5,651,023 | 5,751,022 | 1.2815 | 0.5747 | 8332 | LOC_Os06g10840 | — |
| Chr6 | 5,551,023 | 5,651,022 | 1.2792 | 0.4876 | 5773 | LOC_Os06g10650 | — |
| Chr6 | 7,751,023 | 7,851,022 | 1.1909 | 0.6366 | 8083 | LOC_Os06g13940 | — |

### XI-2 (vs XI-3, outgroup GJ-tmp)
- Windows analysed: 7468
- Mean PBS: 0.0725

**Top PBS windows (population-specific selection):**

| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |
|-----|-------|-----|-----|--------|------|----------|------------|
| Chr7 | 1,901,017 | 2,001,016 | 0.9225 | 0.2056 | 6476 | LOC_Os07g04330 | — |
| Chr10 | 17,351,070 | 17,451,069 | 0.8249 | 0.5308 | 7895 | LOC_Os10g33190 | — |
| Chr3 | 6,801,047 | 6,901,046 | 0.8136 | 0.3508 | 4386 | LOC_Os03g12780 | — |
| Chr11 | 18,501,013 | 18,601,012 | 0.8111 | 0.5108 | 7054 | LOC_Os11g31640 | — |
| Chr3 | 1,851,047 | 1,951,046 | 0.7985 | 0.2469 | 5833 | LOC_Os03g04050 | — |
| Chr3 | 31,001,047 | 31,101,046 | 0.7919 | 0.4575 | 5614 | LOC_Os03g54084 | — |
| Chr10 | 17,301,070 | 17,401,069 | 0.7873 | 0.5463 | 9118 | LOC_Os10g33070 | — |
| Chr6 | 2,151,023 | 2,251,022 | 0.7782 | 0.0978 | 5165 | LOC_Os06g04889 | — |
| Chr4 | 33,101,004 | 33,201,003 | 0.7133 | 0.2114 | 4910 | LOC_Os04g55600 | — |
| Chr7 | 1,851,017 | 1,951,016 | 0.7115 | 0.1597 | 7876 | LOC_Os07g04240 | — |

## 14. Population Tree (UPGMA)

**Method:** UPGMA on Weir & Cockerham Fst

**Newick:** `(((GJ-adm:0.0626,GJ-tmp:0.0626):0.0821,(GJ-sbtrp:0.1286,GJ-trp:0.1286):0.0161):0.1282,((XI-1A:0.0815,(XI-1B:0.0617,(XI-2:0.0481,(XI-3:0.0303,XI-adm:0.0303):0.0178):0.0136):0.0198):0.1004,(cA-Aus:0.1793,(admix:0.1047,cB-Bas:0.1047):0.0746):0.0026):0.0911);`

```
Merge order (UPGMA):
Step       Fst  Cluster 1                  Cluster 2                
----------------------------------------------------------------------
   1    0.0606  XI-3                       XI-adm                   
   2    0.0962  XI-2                       (XI-3:0.0303,XI-adm:0....
   3    0.1233  XI-1B                      (XI-2:0.0481,(XI-3:0.0...
   4    0.1253  GJ-adm                     GJ-tmp                   
   5    0.1629  XI-1A                      (XI-1B:0.0617,(XI-2:0....
   6    0.2093  admix                      cB-Bas                   
   7    0.2572  GJ-sbtrp                   GJ-trp                   
   8    0.2895  (GJ-adm:0.0626,GJ-tmp:...  (GJ-sbtrp:0.1286,GJ-tr...
   9    0.3585  cA-Aus                     (admix:0.1047,cB-Bas:0...
  10    0.3637  (XI-1A:0.0815,(XI-1B:0...  (cA-Aus:0.1793,(admix:...
  11    0.5459  ((GJ-adm:0.0626,GJ-tmp...  ((XI-1A:0.0815,(XI-1B:...
```

## 15. What GraphPop Uniquely Enables

1. **Annotation-conditioned statistics**: `graphpop.diversity(..., {consequence: 'missense_variant'})` computes π restricted to missense variants via Variant→HAS_CONSEQUENCE→Gene traversal — one query vs. 4 separate tools in classical pipelines.
2. **π_N/π_S across all 12 populations**: 288 single-line procedure calls (12 pops × 12 chr × 2 consequences). Classical tools require 288 separate VCF filtering + statistics runs.
3. **Missense-conditioned genome scans**: Sliding-window Fst restricted to missense variants identifies adaptive protein divergence hotspots — not available in any classical tool.
4. **Cross-statistic gene lookups**: After iHS/XP-EHH computation, immediate graph traversal to find associated genes — zero file conversion or coordinate intersection.
