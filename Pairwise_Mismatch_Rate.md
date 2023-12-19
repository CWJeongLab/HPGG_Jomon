# Pairwise Mismatch Rate

**Introduction**

We calculate pairwise mismatch rate (pmr) between Jomon individuals and modern east Asian(Han, Japanese) using the autosomal SNPs in the 1240K panel.

Genotype data of modern individuals(modern_EAS.geno/snp/ind) is extracted from Simons Genome Diversity Project dataset. List of modern east Asian individuals is uploaded as 'modern_EAS.ind'.

## 1. Calculate PMR between Jomon individuals

**input & output file preset**

Input files are HPGG_Jomon.geno/snp/ind.
Output file prefix is as following.

```
fn1="HPGG_Jomon" # input file = HPGG_Jomon.geno/ind/snp
of1=${fn1}".pairwise.mismatch.txt"
tn1="temp1_"${of1}
```

**Extract only autosomal data**

```
paste ${fn1}.geno -d '' | paste ${fn1}.snp - | awk '{if ($2 <= 22) print $7}' | sed 's/./& /g' > ${tn1}_1.geno
cat ${fn1}.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${tn1}_1.ind
```

**Calculate pmr**

```
nind=($(wc -l ${tn1}_1.ind))
let n1=${nind}-1


echo 'ID1 ID2 nSNPs nmismatch pmismatch' | sed s/" "/"\t"/g > ${of1}
for i in $(seq 1 $n1); do
    iid1=($(head -${i} ${tn1}_1.ind | tail -1 | awk '{print $1}'))
    let n2=${i}+1
    for j in $(seq $n2 $nind); do
        iid2=($(head -${j} ${tn1}_1.ind | tail -1 | awk '{print $1}'))
        cut -d ' ' -f ${i},${j} ${tn1}_1.geno | awk '$1 != 9 && $2 != 9' > ${tn1}_2
        ntot=($(wc -l ${tn1}_2))
        nm=($(awk '$1 != $2' ${tn1}_2 | wc -l))
        pm=($(echo ${ntot}" "${nm} | awk '{printf("%.5f\n", $2/$1)}'))
        echo -e ${iid1}"\t"${iid2}"\t"${ntot}"\t"${nm}"\t"${pm} >> ${of1}
    done
    echo ${iid1}" is processed"
done


rm ${tn1}_*
```

As a result, you can get result file like:
```
ID1     ID2     nSNPs   nmismatch       pmismatch
FUN5    FUN23   1025092 182270  0.17781
FUN5    IK002   750510  144629  0.19271
FUN5    NAG019  574471  109662  0.19089
FUN5    TYJ001  160604  31368   0.19531
FUN5    NAG038  453882  87317   0.19238
```

## 2. Calculate PMR between modern individuals

**input & output file preset**

Input files are extracted from Simons Genome Diversity Project dataset and saved as modern_EAS.geno/snp/ind.
Output file prefix is as following.

```
fn1="modern_EAS" # input file = modern_EAS.geno/ind/snp
of1=${fn1}"_1.pairwise.mismatch.txt"
tn1="temp1_"${of1}
```

**Extract only autosomal data**

```
paste ${fn1}.geno -d '' | paste ${fn1}.snp - | awk '{if ($2 <= 22) print $7}' | sed 's/./& /g' > ${tn1}_1.geno
cat ${fn1}.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${tn1}_1.ind
```

**Calculate pmr**

```
nind=($(wc -l ${tn1}_1.ind))
let n1=${nind}-1


echo 'ID1 ID2 nSNPs nmismatch pmismatch' | sed s/" "/"\t"/g > ${of1}
for i in $(seq 1 $n1); do
    iid1=($(head -${i} ${tn1}_1.ind | tail -1 | awk '{print $1}'))
    let n2=${i}+1
    for j in $(seq $n2 $nind); do
        iid2=($(head -${j} ${tn1}_1.ind | tail -1 | awk '{print $1}'))
        cut -d ' ' -f ${i},${j} ${tn1}_1.geno | awk '$1 != 9 && $2 != 9' > ${tn1}_2
        ntot=($(wc -l ${tn1}_2))
	nm1_1=($(grep -c 1 ${tn1}_2))
	nm1=($(echo $nm1_1 | awk '{printf "%.1f", $1 / 2}'))
        nm2=($(grep -v 1 ${tn1}_2 | awk '$1 != $2' | wc -l))
	nm=($(echo $nm1 $nm2 | awk '{printf "%.1f", $1 + $2}'))
        pm=($(echo ${ntot}" "${nm} | awk '{printf("%.5f\n", $2/$1)}'))
        echo -e ${iid1}"\t"${iid2}"\t"${ntot}"\t"${nm}"\t"${pm} >> ${of1}
    done
    echo ${iid1}" is processed"
done


rm ${tn1}_*
```

As a result, you can get result file like:

```
ID1     ID2     nSNPs   nmismatch       pmismatch
A_Han-4.DG      S_Han-1.DG      1119407 274654.0        0.24536
A_Han-4.DG      S_Japanese-2.DG 1116376 274971.5        0.24631
A_Han-4.DG      B_Han-3.DG      1118523 273197.5        0.24425
A_Han-4.DG      S_Japanese-1.DG 1118453 274223.5        0.24518
A_Han-4.DG      S_Han-2.DG      1116990 272497.0        0.24396
```
