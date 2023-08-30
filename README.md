# SharePro for colocalization analysis

SharePro is a command line tool for efficient and accurate colocalization. For analysis conducted in the [SharePro for colocalization analysis paper](https://doi.org/10.1101/2023.07.24.550431), please refer to [SharePro_coloc_analysis](https://github.com/zhwm/sharepro_coloc_analysis).

## Overview 

Colocalization analysis is a commonly used statistical procedure for assessing whether two or more traits share the same genetic signals identified in genome-wide association studies (GWAS). It is important for understanding the interplay between heritable traits.
SharePro takes marginal associations (z-scores) from GWAS summary statistics and linkage disequilibrium (LD) infromation as inputs, and infers posterior probability of colocalization. Unlike existing methods, SharePro takes an effect group-level approach to integrate LD modelling and colocalization assessment to account for multiple causal variants in colocalization analysis.

<p align="center">
  <img src="doc/SharePro_loc.png" alt="example image">
  <br>
  <em>Figure 1: SharePro overview.</em>
</p>

## Installation

SharePro was developed under Python 3.9.7 but should be compatible with other versions of Python 3. The following Python modules are required:

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [pandas](https://pandas.pydata.org/getpandas.html)

To install SharePro for colocalization analysis:

```
git clone https://github.com/zhwm/SharePro_coloc.git
cd SharePro_coloc
pip install -r requirements.txt 
``` 

To test the installation and display basic usages:
```
python sharepro_loc.py -h
```

## Input files

Example input files are included in the [dat](dat/) directory.

SharePro takes in a summary file as well as z-scores files and LD files as inputs.

1. The **summary file** contains two mandatory columns: names of z-score files and ld files. Multiple files are allowed and should be separated by comma. An example can be found at [dat/BMD_RSPO3.zld](dat/BMD_RSPO3.zld).

2. The **zscore files** contain two mandatory columns: variant IDs and z-scores. Here are examples for [BMD](dat/BMD_SH.txt) GWAS and [RSPO3](dat/RSPO3_SH.txt) pQTL zscore files.

3. The **LD files** contain Pearson correlation between variants. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/RSPO3.ld](dat/RSPO3.ld) and a working script for matching raw GWAS summary statistics and PLINK bim file is [provided](match_bim_ss.py).

## Usage examples

We use the colocalization analysis of a RSPO3 cis-pQTL and eBMD GWAS as an example. The [dat/](dat/) folder contains all files required for this example.
If you want to learn more about this locus, please refer to the [analysis repo](https://github.com/zhwm/SharePro_coloc_analysis/tree/main/dat).

We use `--zld` to indicate the path to the summary file and `--zdir` to indicate the path to the zscore files.
Additionally, we can specify the sample sizes of the eBMD GWAS study and RSPO3 pQTL study with `--N`.
Next, we specify path to save results with `--save`. The max number of causal signals is set to 10 using `--K`.

```
python sharepro_loc.py \
--zld dat/BMD_RSPO3.zld \
--zdir dat \
--N 426824 10708 \
--save res \
--verbose \
--K 10
```

## Output files

Here are the expected outputs from the example above:

```
**********************************************************************
* SharePro for accurate and efficient colocalization                 *
* Version 2.0.0                                                      *
* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *
**********************************************************************
LD list with 1 LD blocks loaded

processing BMD_SH.txt,RSPO3_SH.txt
**********************************************************************
Iteration-->0 . Likelihood: 1265.8 . KL_b: -7.0 . KL_c: -23.0 . KL_s: 0.7 . ELBO: 1236.4
**********************************************************************
Iteration-->1 . Likelihood: 1265.8 . KL_b: -7.0 . KL_c: -23.0 . KL_s: 0.7 . ELBO: 1236.4
**********************************************************************
Iteration-->0 . Likelihood: 1453.7 . KL_b: -25.1 . KL_c: -104.0 . KL_s: 40.6 . ELBO: 1365.2
**********************************************************************
Iteration-->1 . Likelihood: 1458.1 . KL_b: -23.0 . KL_c: -94.9 . KL_s: 39.2 . ELBO: 1379.5
**********************************************************************
Iteration-->2 . Likelihood: 1459.5 . KL_b: -22.4 . KL_c: -92.4 . KL_s: 39.1 . ELBO: 1383.8
**********************************************************************
Iteration-->3 . Likelihood: 1459.7 . KL_b: -22.4 . KL_c: -92.3 . KL_s: 39.1 . ELBO: 1384.1
**********************************************************************
Iteration-->4 . Likelihood: 1459.7 . KL_b: -22.4 . KL_c: -92.2 . KL_s: 39.1 . ELBO: 1384.1
Colocalization probability: 1.0
```

The results have been saved into the **colocalization summary** (cs) files with columns:
`cs` for variant representations in effect groups; 
`share` for colocalization probabilities;
`variantProb` for variant representation weights in effect groups;
`k` for the max number of causal signals specified in the model.

```
$> cat res/BMD_SH.txt_RSPO3_SH.txt.cs 
cs      share   variantProb     k
rs7741021/rs9482773     1.0     0.5381/0.4619   1
rs7741021/rs9482773     1.0     0.594/0.406     10
rs853974        1.0     1.0     10
rs577721086/rs72959041  0.0097  0.5313/0.4687   10
rs717796/rs7775814/rs1512450/rs7775090/rs1512449/rs910536/rs6902741/rs1569870/rs7738255/rs1080708/rs7756072     0.0001  0.1217/0.1167/0.1018/0.0879/0.0867/0.085/0.0785/0.0747/0.0705/0.0705/0.0639     10
rs11759578/rs77525683/rs73593094/rs73593068/rs73577838/rs9482768        0.0003  0.6213/0.2032/0.0666/0.031/0.0264/0.0131        10
rs2800728/rs727330/rs727331/rs2800720/rs2800721/rs10456964/rs2800733/rs2800727/rs2745351/rs2745356/rs2800732/rs2745355/rs727332/rs2800729/rs2800730/rs719728/rs2800719/rs2745354/rs2800718/rs2800723/rs2800722  0.0002  0.0688/0.0665/0.0658/0.0647/0.0608/0.0533/0.0521/0.0514/0.0463/0.0412/0.041/0.0404/0.0399/0.0394/0.0386/0.0379/0.034/0.0327/0.0324/0.0322/0.0322 10
```

## Versions
* Version 1.0.0 (2023/01/25) Initial release
* Version 2.0.0 (2023/08/25) Added K=1 in case of LD mismatch

## Citations

If you find SharePro for colocalization analysis useful, please cite:

[Wenmin Zhang, Tianyuan Lu, Robert Sladek, Yue Li, Hamed Najafabadi, Josée Dupuis. SharePro: an accurate and efficient genetic colocalization method accounting for multiple causal signals.](https://doi.org/10.1101/2023.07.24.550431)