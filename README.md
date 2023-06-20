# SharePro for colocalization analysis

SharePro is a command line tool for efficient and accurate colocalization. For analysis conducted in the SharePro for colocalization analysis paper, please refer to [SharePro_coloc_analysis](https://github.com/zhwm/sharepro_coloc_analysis).

## Overview 

Colocalization analysis is a commonly used statistical procedure for assessing whether two traits share the same genetic signals identified in genome-wide association studies (GWAS). It is important for understanding the interplay between heritable traits.
SharePro takes marginal associations (z-scores) from GWAS summary statistics and matched LD matrix as inputs, and infers posterior probability of colocalization. Different from the variant-level approach and the locus-level approach in existing methods, SharePro takes an effect group-level approach to model uncertainties originating from LD and multiple causal signals in colocalization analysis. 

<p align="center">
  <img src="doc/SharePro_loc.png" alt="example image">
  <br>
  <em>Figure 1: SharePro overview.</em>
</p>

## Installation

SharePro was developed under Python 3.9.7 environment but should be compatible with older versions of Python 3. The following Python modules are required:

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [pandas](https://pandas.pydata.org/getpandas.html)

To install SharePro for colocalization analysis:

```
git clone https://github.com/zhwm/SharePro_coloc.git
cd SharePro_coloc
pip install -r requirements.txt 
``` 

To test the installation and display basic usage:
```
python sharepro_loc.py -h
```

## Input files

Example input files are included in the [dat](dat/) directory.

SharePro takes in a summary file with path to z-score files and LD files, z-scores files, LD files as inputs.

1. **a summary file** contains two mandatory columns: names of z-score file and ld files. Multiple files are allowed and should be separated by comma. An example can be found at [dat/BMD_RSPO3.zld](dat/BMD_RSPO3.zld).

2. **zscore files** that contain two mandatory columns: variant IDs and z-scores. Here are examples for [BMD](dat/BMD_SH.txt) and [RSPO3](dat/RSPO3_SH.txt) pQTL zscore files.

3. **LD files** that contain correlation coefficient matrix. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/RSPO3.ld](dat/RSPO3.ld) and a working script for matching raw GWAS summary statistics and PLINK bim file is provided [here](match_bim_ss.py).

## Usage examples

Here we use the colocalization analysis of a cis RSPO3 pQTL locus and eBMD GWAS locus as an example. In the [dat/](dat/) folder we have provided all files needed to run this example.
If you want to learn more about this locus, please refer to the [analysis repo](https://github.com/zhwm/SharePro_coloc_analysis/tree/main/dat).

Here we use `--zld` to indicate path to the summary file and `--zdir` to indicate path to zscore files.
Additionally, we specify the sample sizes of both the eBMD GWAS study and RSPO3 pQTL study with `--N`.
We use `--save` to specify path to save result and `--prefix` to specify prefix of output files. We set the max number of causal signals as 10 with `--K`.

```
python3 sharepro_loc.py \
--zld dat/BMD_RSPO3.zld \
--zdir dat \
--N 426824 10708 \
--save res \
--prefix BMD_RSPO3 \
--verbose \
--K 10
```

## Output interpretation

With the example above, we can obtain the following output. The first section records evidence lower bound (ELBO) of our algorithm and to guarantee convergence, ELBO should always increase at each iteration.

The (total probabilities -1) gives us an estimate of total number of causal signals. In this example, we have identified 6 causal signals. Among them, the first two causal signals are shared while the additional causal signals are trait-specific.

```
**********************************************************************
* SharePro for accurate and efficient colocalization                 *
* Version 1.0.0                                                      *
* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *
**********************************************************************
LD list with 1 LD blocks loaded

processing BMD_SH.txt,RSPO3_SH.txt
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
**********************************************************************
Iteration-->5 . Likelihood: 1459.7 . KL_b: -22.4 . KL_c: -92.2 . KL_s: 39.1 . ELBO: 1384.1
Attainable coverage for effect groups: [1.   1.   1.   0.99 1.   1.   0.17 0.   0.   0.78]
The 0-th effect group contains effective variants:
causal variants: ['rs7741021', 'rs9482773']
variant probabilities for this effect group: [0.5941, 0.4059]
shared probability for this effect group: 1.0

The 1-th effect group contains effective variants:
causal variants: ['rs853974']
variant probabilities for this effect group: [1.0]
shared probability for this effect group: 1.0

The 2-th effect group contains effective variants:
causal variants: ['rs577721086', 'rs72959041']
variant probabilities for this effect group: [0.5313, 0.4687]
shared probability for this effect group: 0.0094

The 3-th effect group contains effective variants:
causal variants: ['rs717796', 'rs7775814', 'rs1512450', 'rs7775090', 'rs1512449', 'rs910536', 'rs6902741', 'rs1569870', 'rs7738255', 'rs1080708', 'rs7756072']
variant probabilities for this effect group: [0.1217, 0.1167, 0.1018, 0.0879, 0.0867, 0.085, 0.0785, 0.0747, 0.0706, 0.0705, 0.0639]
shared probability for this effect group: 0.0001

The 4-th effect group contains effective variants:
causal variants: ['rs11759578', 'rs77525683', 'rs73593094', 'rs73593068', 'rs73577838', 'rs9482768']
variant probabilities for this effect group: [0.6213, 0.2032, 0.0666, 0.031, 0.0264, 0.0131]
shared probability for this effect group: 0.0003

The 5-th effect group contains effective variants:
causal variants: ['rs2800728', 'rs727330', 'rs727331', 'rs2800720', 'rs2800721', 'rs10456964', 'rs2800733', 'rs2800727', 'rs2745351', 'rs2745356', 'rs2800732', 'rs2745355', 'rs727332', 'rs2800729', 'rs2800730', 'rs719728', 'rs2800719', 'rs2745354', 'rs2800718', 'rs2800723', 'rs2800722']
variant probabilities for this effect group: [0.0689, 0.0665, 0.0658, 0.0647, 0.0608, 0.0533, 0.0521, 0.0514, 0.0463, 0.0412, 0.041, 0.0404, 0.0399, 0.0394, 0.0386, 0.0379, 0.034, 0.0327, 0.0324, 0.0322, 0.0322]
shared probability for this effect group: 0.0002
```

## Output files

1. **colocalization summary** (cs) file contains seven columns: 
`cs` for variant representation in effect groups; 
`share` for colocalization probability;
`variantProb` for variant representation weight in effect groups.

```
$> cat res/BMD_SH.txt_RSPO3_SH.txt.cs 
cs      share   variantProb
rs7741021/rs9482773     1.0     0.5941/0.4059
rs853974        1.0     1.0
rs577721086/rs72959041  0.0094  0.5313/0.4687
rs717796/rs7775814/rs1512450/rs7775090/rs1512449/rs910536/rs6902741/rs1569870/rs7738255/rs1080708/rs7756072     0.0001  0.1217/0.1167/0.1018/0.0879/0.0867/0.085/0.0785/0.0747/0.0706/0.0705/0.0639
rs11759578/rs77525683/rs73593094/rs73593068/rs73577838/rs9482768        0.0003  0.6213/0.2032/0.0666/0.031/0.0264/0.0131
rs2800728/rs727330/rs727331/rs2800720/rs2800721/rs10456964/rs2800733/rs2800727/rs2745351/rs2745356/rs2800732/rs2745355/rs727332/rs2800729/rs2800730/rs719728/rs2800719/rs2745354/rs2800718/rs2800723/rs2800722  0.0002  0.0689/0.0665/0.0658/0.0647/0.0608/0.0533/0.0521/0.0514/0.0463/0.0412/0.041/0.0404/0.0399/0.0394/0.0386/0.0379/0.034/0.0327/0.0324/0.0322/0.0322
```

2. **hyperparamters summary** (h2) file adds two additional columns in the summary file to record the heritability and effect size variance estimates used in the colocalization algorithms.

```
$> cat BMD_RSPO3.h2              
z	ld	h2	varb
BMD_SH.txt,RSPO3_SH.txt	RSPO3.ld,RSPO3.ld	4.86e-03,7.84e-02	2.77e-04,6.76e-03
```
