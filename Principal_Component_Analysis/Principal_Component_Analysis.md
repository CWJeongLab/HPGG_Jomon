# Principal Component Analysis

**Introduction**

We performed PCA using the smartpca v18140 from the Eigensoft v8.0.0 package using 2,077 present-day Eurasian individuals (149 populations) genotyped on the Affymetrix Axiom Genome-wide Human Origins 1 array. We used 593,124 autosomal SNPs of HumanOrigins SNP panel overlapping with the 1240K set. This SNP set is uploaded under the file name "[1240KHO.auto.snp](https://github.com/jkc0326/jomon/blob/main/Principal_Component_Analysis/1240KHO.auto.230629.snp)". 

Genotype data of modern individuals(modern_HO.geno/snp/ind) is extracted from Simons Genome Diversity Project dataset. "[modern.pops](https://github.com/jkc0326/jomon/blob/main/Principal_Component_Analysis/modern.pops)" file includes the list of 149 modern population.

We projected ancient Jomon-related individuals not included in PC calculation using the ‘lsqrproject: YES’ option.

For PCA of ancient Jomon-related individuals only, we first computed covariance between each pair of individuals using SNPs not missing in both individuals to avoid an issue due to varying level of genotype missingness across individuals using R v4.2.0. Then, we performed PCA on the calculated covariance matrix.

## 1. PCA projection on modern populations

**Prepare parameter file**

The parameter file format of smartpca is as follows:

```
genotypename: <Path_to_file>/<Genotype_Data>.geno
snpname: <Path_to_file>/<Genotype_Data>.snp
indivname: <Path_to_file>/<Genotype_Data>.ind
evecoutname: <Path_to_file>/<Output_file_prefix>.evec
evaloutname: <Path_to_file>/<Output_file_prefix>.eval
poplistname: <Path_to_file>/moderns.pops
altnormstype: NO
numoutevec: 20
numoutlieriter: 0
numoutlierevec: 0
outliersigmathresh: 6.0
numthreads: 8
qtmode: 0
lsqproject: YES
```

**Run smartpca**

After preparing parameter file, run smartpca with ```smartpca -p <Parameter_File> > <Output_file_prefix>.log```.

As a result, you can get ```*.eval``` and ```*.evec``` file which includes the information of eigenvector and eigenvalue, respectively.

```*.eval``` file is like:
```
  104.218952
   13.870599
    6.886297
    6.137649
    4.445190
```

```*.evec``` file is like:
```
#eigvals:   104.219    13.871     6.886     6.138     4.445     4.051     3.389     3.204     2.924     2.766     2.609     2.532     2.454     2.305     2.236     2.226     2.195     2.182     2.161     2.124
           IREJ-T006    -0.0155     -0.0100      0.0364      0.0021      0.0083     -0.0037      0.0002      0.0003     -0.0013     -0.0006     -0.0037      0.0038      0.0036      0.0062      0.0036     -0.0015      0.0031      0.0049      0.0012     -0.0006  Iran_Non-Zoroastrian_Fars
           IREJ-T009    -0.0162     -0.0089      0.0368      0.0029      0.0120     -0.0029      0.0003      0.0009     -0.0022      0.0032      0.0038      0.0015      0.0004      0.0052     -0.0011      0.0014     -0.0001     -0.0010      0.0037     -0.0015  Iran_Non-Zoroastrian_Fars
           IREJ-T022    -0.0163     -0.0070      0.0343      0.0028      0.0065     -0.0031      0.0003      0.0009     -0.0033     -0.0000      0.0022      0.0012     -0.0004     -0.0030     -0.0008      0.0012      0.0015      0.0007     -0.0022      0.0025  Iran_Non-Zoroastrian_Fars
           IREJ-T023    -0.0168     -0.0103      0.0384      0.0014      0.0059     -0.0042      0.0006     -0.0008     -0.0018      0.0007     -0.0012      0.0006     -0.0001      0.0066     -0.0001      0.0018      0.0011     -0.0042      0.0042     -0.0012  Iran_Non-Zoroastrian_Fars
           IREJ-T026    -0.0171     -0.0086      0.0383      0.0007      0.0064     -0.0019     -0.0005     -0.0019      0.0029      0.0002      0.0002     -0.0029      0.0035     -0.0009     -0.0009      0.0028      0.0013      0.0042      0.0026      0.0005  Iran_Non-Zoroastrian_Fars
```

## 2. PCA with only Jomon individuals

We performed PCA with only Jomon individuals in R v4.2.0. We used ```cov()``` function with ```use="pairwise.complete.obs"``` for computed covariance between each pair of Jomon individuals using SNPs not missing in both individuals, and ```eigen()``` function to calculate PCs.

We excluded populations with low Jomon ancestry proportions(Janghang_6700BP and Yeondaedo_7000BP).



The code below is written in R.
```
df_geno <- read.table('<Genotype_Data>.geno')
df_ind <- read.table('<Genotype_Data>.ind')

Mcov <- cov(df_geno, use = "pairwise.complete.obs", method = "pearson")
df_evals <- eigen(Mcov)$values
df_evecs <- as.matrix(cbind(unlist(eigen(Mcov)$vectors), df_ind$V3))
```
