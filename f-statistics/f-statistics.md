# f-statistics

**Introduction**

 We calculated f-statistics by qp3Pop and f4 functions from the R library ADMIXTOOLS2 v2.0.0 using 1240K SNP set. For calculating both outgroup-f3 and f4 statistics, we used Mbuti, the central African rainforest hunter-gatherer population, as and outgroup. Genotype data of Mbuti is extracted from Simons Genome Diversity Project dataset.

 

## 1. Outgroup-f3

 We calculated f3(Mbuti; Jomon1, Jomon2) for all pairs of Jomon groups to measure shared genetic drift between two Jomon groups.

 The code below is written in R.
```
library(admixtools)

genotype_dir <- '<Path_to_file>/<Genotype_Data>'
Jomon_list <- c("Funadomari_3700BP", "Ikawazu_2600BP", "JpKa", "JpFu", "JpHi", 
                "JpKo", "JpOd", "Nagabaka_2800BP", "Nagabaka_4000BP",
                "Rokutsu_3500BP", "Yokjido_4000BP")
outgroup = "Mbuti.DG"

df_f3 <- subset(qp3pop(genotype_dir, outgroup, Jomon_list, Jomon_list), pop2!=pop3)

write.table(df_f3, file="outgroup_f3.txt", quote=F, sep="\t", row.names=F)
```

 As a result, you can get ```outgroup_f3.txt``` file like :
```
pop1    pop2    pop3    est     se      z       p       n
Mbuti.DG        Funadomari_3700BP       Ikawazu_2600BP  0.408609220358147       0.00381605376846017     107.076379199716        0       836968
Mbuti.DG        Funadomari_3700BP       JpFu    0.407105768426452       0.00408751118863638     99.597468884792 0       695324
Mbuti.DG        Funadomari_3700BP       JpHi    0.408557577617262       0.00418316797787465     97.6670264684992        0       460263
Mbuti.DG        Funadomari_3700BP       JpKa    0.402270832731516       0.0039058070339256      102.993012516342        0       1145010
Mbuti.DG        Funadomari_3700BP       JpKo    0.408735699364518       0.00370862765906623     110.212115353589        0       1113386
```

## 2. f4-statistics between Jomon groups

 We calculated f4(Mbuti, Jomon1; Jomon2, Jomon3) for all pairs of Jomon and Jomon-related groups to figure out the population structure within the Jomon-related groups.

 The code below is written in R.
```
library(admixtools)

genotype_dir <- '<Path_to_file>/<Genotype_Data>'
Jomon_list <- c("Funadomari_3700BP", "Ikawazu_2600BP", "JpKa", "JpFu", "JpHi", 
                "JpKo", "JpOd", "Nagabaka_2800BP", "Nagabaka_4000BP",
                "Rokutsu_3500BP", "Yokjido_4000BP", "Janghang_6700BP",
                "Yeondaedo_7000BP")
outgroup = "Mbuti.DG"

df_f4 <- subset(f4(genotype_dir, outgroup, Jomon_list, Jomon_list, Jomon_list), pop2!=pop3 & pop2!=pop4 & pop3!=pop4)

write.table(df_f4, file="f4_between_Jomon.txt", quote=F, sep="\t", row.names=F)
```

 As a result, you can get ```f4_between_Jomon.txt``` file like :
```
pop1    pop2    pop3    pop4    est     se      z       p       n
Mbuti   Funadomari_3700BP       Ikawazu_2600BP  JpFu    -0.000955540334683795   0.00055571895859203     -1.71946686343895       0.0855293963676676      532164
Mbuti   Funadomari_3700BP       Ikawazu_2600BP  JpHi    -6.30381610387181e-06   0.000561977229884355    -0.0112172091121361     0.991050149721505       401212
Mbuti   Funadomari_3700BP       Ikawazu_2600BP  JpKa    -0.00208450206175681    0.000565327251940911    -3.68724850005053       0.000226691958449509    834240
Mbuti   Funadomari_3700BP       Ikawazu_2600BP  JpKo    -0.000407204791121813   0.00046931690480545     -0.867654216058156      0.385583656833513       813616
Mbuti   Funadomari_3700BP       Ikawazu_2600BP  JpOd    0.000388503722492309    0.000430167890853527    0.90314440187865        0.366449261063249       827925
```

## 3. f4-statistics against world-wide populations

 We calculated f4(Mbuti, worldwide; Jomon 1, Jomon 2) for testing genetic symmetry between Jomon groups or searching for additional admixture sources over worldwide populations.

 The list of modern world-wide populations is uploaded under the file name "f4_worldwide_pops.txt".

 The code below is written in R.
```
library(admixtools)

genotype_dir <- '<Path_to_file>/<Genotype_Data>'
Jomon_list <- c("Funadomari_3700BP", "Ikawazu_2600BP", "JpKa", "JpFu", "JpHi", 
                "JpKo", "JpOd", "Nagabaka_2800BP", "Nagabaka_4000BP",
                "Rokutsu_3500BP", "Yokjido_4000BP")
worldwide_list <- read.table("f4_worldwide_pops.txt", header=F)$V1
outgroup = "Mbuti.DG"

df_f4 <- subset(f4(genotype_dir, outgroup, worldwide_list, Jomon_list, Jomon_list), pop1!=pop2, pop3!=pop4)

write.table(df_f4, file="f4_worldwide.txt", quote=F, sep="\t", row.names=F)
```

 As a result, you can get ```f4_worldwide.txt``` file like :
```
pop1 pop2 pop3 pop4 est se z p n
Mbuti.DG Abkhasian.DG Funadomari_3700BP Ikawazu_2600BP 6.98823363181106e-05 0.00032994983208434 0.211796853711529 0.832265523217229 816253
Mbuti.DG Adygei.DG Funadomari_3700BP Ikawazu_2600BP -0.000172555673989129 0.000326817631166444 -0.527987652848657 0.597507905752838 816210
Mbuti.DG Albanian.DG Funadomari_3700BP Ikawazu_2600BP 0.000551871506472877 0.000401746337019497 1.37368148908871 0.169540586075834 815377
Mbuti.DG Aleut.DG Funadomari_3700BP Ikawazu_2600BP -0.000589990783515519 0.000357062688364088 -1.65234510001202 0.0984642203270847 816291
Mbuti.DG Altaian.DG Funadomari_3700BP Ikawazu_2600BP 0.000377688968739299 0.000423004957505873 0.892871258451044 0.371926120385664 814857
```
