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
python sharepro_loc.py \
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
Iteration-->0 . Likelihood: 1448.9 . KL_b: -29.6 . KL_c: -99.2 . KL_s: 42.4 . ELBO: 1362.5
**********************************************************************
Iteration-->1 . Likelihood: 1455.9 . KL_b: -29.5 . KL_c: -96.3 . KL_s: 39.6 . ELBO: 1369.7
**********************************************************************
Iteration-->2 . Likelihood: 1458.3 . KL_b: -30.0 . KL_c: -92.3 . KL_s: 39.3 . ELBO: 1375.3
**********************************************************************
Iteration-->3 . Likelihood: 1459.0 . KL_b: -30.4 . KL_c: -92.2 . KL_s: 39.2 . ELBO: 1375.7
**********************************************************************
Iteration-->4 . Likelihood: 1459.1 . KL_b: -30.5 . KL_c: -92.1 . KL_s: 39.2 . ELBO: 1375.7
**********************************************************************
Iteration-->5 . Likelihood: 1459.1 . KL_b: -30.5 . KL_c: -92.1 . KL_s: 39.2 . ELBO: 1375.7
total probabilities 6.9374
entropy of posterior distribution:[0.68 0.   0.69 2.69 1.37 3.06 7.68 7.68 7.68 7.68]
The 0-th effect group contains effective variants:
causal variants: ['rs853974']
variant probabilities for this effect group: [1.]
shared probability for this effect group: 1.0
specific probabilities for this effect group: [0.0, 0.0]
causal probabilities of this effect group for traits: [1.0, 1.0]
causal effect sizes for traits: [-0.0108, -0.0671]

The 1-th effect group contains effective variants:
causal variants: ['rs7741021', 'rs9482773']
variant probabilities for this effect group: [0.5895 0.4105]
shared probability for this effect group: 1.0
specific probabilities for this effect group: [0.0, 0.0]
causal probabilities of this effect group for traits: [1.0, 1.0]
causal effect sizes for traits: [-0.0702, -0.2367]

The 2-th effect group contains effective variants:
causal variants: ['rs577721086', 'rs72959041']
variant probabilities for this effect group: [0.5322 0.4677]
shared probability for this effect group: 0.0044
specific probabilities for this effect group: [0.0, 0.9956]
causal probabilities of this effect group for traits: [0.0044, 1.0]
causal effect sizes for traits: [-0.0044, -0.0822]

The 3-th effect group contains effective variants:
causal variants: ['rs11759578', 'rs77525683', 'rs73593094', 'rs73593068', 'rs73577838', 'rs75607329']
variant probabilities for this effect group: [0.591  0.2055 0.0724 0.0357 0.0301 0.0153]
shared probability for this effect group: 0.0001
specific probabilities for this effect group: [0.0, 0.9999]
causal probabilities of this effect group for traits: [0.0001, 1.0]
causal effect sizes for traits: [0.0015, -0.0651]

The 4-th effect group contains effective variants:
causal variants: ['rs717796', 'rs7775814', 'rs1512450', 'rs7775090', 'rs1512449', 'rs910536', 'rs6902741', 'rs1569870', 'rs7738255', 'rs1080708', 'rs7756072']
variant probabilities for this effect group: [0.1194 0.1146 0.1003 0.0873 0.0858 0.0846 0.0779 0.0743 0.0703 0.0702 0.0638]
shared probability for this effect group: 0.0
specific probabilities for this effect group: [1.0, 0.0]
causal probabilities of this effect group for traits: [1.0, 0.0]
causal effect sizes for traits: [0.0114, -0.0016]

The 5-th effect group contains effective variants:
causal variants: ['rs2800728', 'rs727330', 'rs727331', 'rs2800720', 'rs2800721', 'rs10456964', 'rs2800733', 'rs2800727', 'rs2745351', 'rs2800732', 'rs2745356', 'rs2745355', 'rs727332', 'rs2800729', 'rs2800730', 'rs719728', 'rs2800719', 'rs2745354', 'rs2800718', 'rs2800723', 'rs2800722']
variant probabilities for this effect group: [0.0668 0.065  0.0644 0.0634 0.0599 0.0527 0.0518 0.0512 0.0465 0.0415 0.0414 0.0407 0.0402 0.0399 0.0393 0.0386 0.0347 0.0335 0.0331 0.033  0.033 ]
shared probability for this effect group: 0.0
specific probabilities for this effect group: [1.0, 0.0]
causal probabilities of this effect group for traits: [1.0, 0.0]
causal effect sizes for traits: [0.0119, 0.0049]
```

## Output files

1. **colocalization summary** (cs) file contains seven columns: 
`cs` for variant representation in effect groups; 
`totalProb` for overall probability weights for effect groups; 
`share` for colocalization probability;
`beta` for causal effect sizes;
`specific` for trait specific probabilities;
`causalProb` for causal probabilities; 
`variantProb` for variant representation weight in effect groups.

```
$> cat BMD_SH.txt_RSPO3_SH.txt.cs 
cs	totalProb	share	beta	specific	causalProb	variantProb
rs853974	1.0	1.0	-0.0108,-0.0671	0.0,0.0	1.0,1.0	1.0
rs7741021/rs9482773	1.0	1.0	-0.0702,-0.2367	0.0,0.0	1.0,1.0	0.567/0.433
rs577721086/rs72959041	0.9999	0.0042	-0.0044,-0.0821	0.0,0.9958	0.0042,1.0	0.5327/0.4672
rs11759578/rs77525683/rs73593094/rs73593068/rs73577838/rs75607329	0.9501	0.0001	0.0015,-0.0651	0.0,0.9999	0.0001,1.0	0.592/0.2079/0.0723/0.034/0.0293/0.0146
rs717796/rs7775814/rs7775090/rs1512450/rs910536/rs1512449/rs6902741/rs1569870/rs1080708/rs7738255/rs7756072	0.9493	0.0	0.0114,-0.0017	1.0,0.0	1.0,0.0	0.1065/0.1/0.0985/0.0979/0.0953/0.0878/0.0796/0.0769/0.0706/0.0705/0.0657
rs727330/rs2800728/rs727331/rs2800720/rs2800721/rs2800727/rs10456964/rs2800733/rs2800730/rs2745351/rs2800729/rs2800732/rs719728/rs2745356/rs2745355/rs727332/rs2800719/rs2745354/rs2800718/rs10456965/rs2800723	0.9716	0.0	0.012,0.0042	1.0,0.0	1.0,0.0	0.0711/0.0696/0.0691/0.0687/0.0643/0.0574/0.0571/0.0468/0.0448/0.0447/0.0447/0.042/0.0404/0.0354/0.0348/0.0342/0.0308/0.0296/0.029/0.0287/0.0284

```

2. **hyperparamters summary** (h2) file adds two additional columns in the summary file to record the heritability and effect size variance estimates used in the colocalization algorithms.

```
$> cat BMD_RSPO3.h2              
z	ld	h2	varb
BMD_SH.txt,RSPO3_SH.txt	RSPO3.ld,RSPO3.ld	4.86e-03,7.84e-02	2.77e-04,6.76e-03
```
