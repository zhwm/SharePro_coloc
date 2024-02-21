# SharePro for colocalization analysis

SharePro is a command line tool for efficient and accurate colocalization. For analysis conducted in the [SharePro for colocalization analysis paper](https://doi.org/10.1101/2023.07.24.550431), please refer to [SharePro_coloc_analysis](https://github.com/zhwm/sharepro_coloc_analysis).

## Overview 

Colocalization analysis is a commonly used statistical procedure for assessing whether two or more traits share the same genetic signals identified in genome-wide association studies (GWAS). It is important for understanding the interplay between heritable traits.
SharePro takes marginal associations from GWAS summary statistics and linkage disequilibrium (LD) information as inputs, and infers posterior probability of colocalization. Unlike existing methods, SharePro takes an effect group-level approach to integrate LD modelling and colocalization assessment to account for multiple causal variants in colocalization analysis.

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

### Run SharePro with a batch of zscore files and ld files

SharePro can take in a summary file as well as z-scores files and LD files as inputs.

1. The **summary file** contains two mandatory columns: names of z-score files and ld files. Multiple files are allowed and should be separated by comma. An example can be found at [dat/BMD_RSPO3.zld](dat/BMD_RSPO3.zld).

2. The **zscore files** contain two mandatory columns: variant IDs and z-scores. Here are examples for [BMD](dat/BMD_SH.txt) GWAS and [RSPO3](dat/RSPO3_SH.txt) pQTL zscore files.

3. The **LD files** contain Pearson correlation between variants. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/RSPO3.ld](dat/RSPO3.ld) and a working script for matching raw GWAS summary statistics and PLINK bim file is [provided](match_bim_ss.py).

### Run SharePro in one locus with GWAS (BETA/SE/N) and ld files

SharePro can also use the matched GWAS summary statistics including effect size (BETA), standard error (SE) and sample size (N) and ld matrix from one locus. 

1. The **GWAS file** contains three mandatory columns: BETA/SE/N. Here are examples for [BMD](dat/BMD_bse.txt) GWAS and [RSPO3](dat/RSPO3_bse.txt) pQTL zscore files.

2. The **LD file** is the same as above.


## Usage examples

We use the colocalization analysis of a RSPO3 cis-pQTL and eBMD GWAS as an example. The [dat/](dat/) folder contains all files required for this example.
If you want to learn more about this locus, please refer to the [analysis repo](https://github.com/zhwm/SharePro_coloc_analysis/tree/main/dat).

### Run SharePro with a batch of zscore files and ld files

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

### Run SharePro in one locus with GWAS (BETA/SE/N) and ld files

Alternatively, we can use `--z` to provide path to the GWAS summary statistic files with three mandatory columns: BETA/SE/N. With `--ld`, we can provide path to the ld files. Lastly, we use `--save` to provide the filename for the colocalization result.

```
python sharepro_bse.py \
--z dat/BMD_bse.txt dat/RSPO3_bse.txt \
--ld dat/RSPO3.ld dat/RSPO3.ld \
--save res.txt \
--verbose \
--K 10
```

## Output files

Here are the expected outputs from the example above:

```
K = 10
**********************************************************************
Iteration-->0 . Likelihood: 1453.7 . KL_b: -25.1 . KL_c: -104.0 . KL_s: 40.6 . ELBO: 1365.2
**********************************************************************
Iteration-->1 . Likelihood: 1458.1 . KL_b: -23.0 . KL_c: -94.9 . KL_s: 39.2 . ELBO: 1379.5
**********************************************************************
Iteration-->2 . Likelihood: 1459.5 . KL_b: -22.4 . KL_c: -92.4 . KL_s: 39.1 . ELBO: 1383.8
**********************************************************************
Iteration-->3 . Likelihood: 1459.7 . KL_b: -22.4 . KL_c: -92.3 . KL_s: 39.1 . ELBO: 1384.1
Colocalization probability: 1.0
The 0-th effect group contains effective variants:
causal variants: ['rs7741021', 'rs9482773']
variant probabilities for this effect group: [0.5936, 0.4064]
shared probability for this effect group: 1.0

The 1-th effect group contains effective variants:
causal variants: ['rs853974']
variant probabilities for this effect group: [1.0]
shared probability for this effect group: 1.0

The 2-th effect group contains effective variants:
causal variants: ['rs577721086', 'rs72959041']
variant probabilities for this effect group: [0.5313, 0.4686]
shared probability for this effect group: 0.0114

The 3-th effect group contains effective variants:
causal variants: ['rs717796', 'rs7775814', 'rs1512450', 'rs7775090', 'rs1512449', 'rs910536', 'rs6902741', 'rs1569870', 'rs7738255', 'rs1080708', 'rs7756072']
variant probabilities for this effect group: [0.1215, 0.1166, 0.1016, 0.0878, 0.0866, 0.0849, 0.0784, 0.0746, 0.0705, 0.0704, 0.0638]
shared probability for this effect group: 0.0001

The 4-th effect group contains effective variants:
causal variants: ['rs11759578', 'rs77525683', 'rs73593094', 'rs73593068', 'rs73577838', 'rs9482768']
variant probabilities for this effect group: [0.6213, 0.2032, 0.0666, 0.031, 0.0264, 0.0131]
shared probability for this effect group: 0.0003

The 5-th effect group contains effective variants:
causal variants: ['rs2800728', 'rs727330', 'rs727331', 'rs2800720', 'rs2800721', 'rs10456964', 'rs2800733', 'rs2800727', 'rs2745351', 'rs2745356', 'rs2800732', 'rs2745355', 'rs727332', 'rs2800729', 'rs2800730', 'rs719728', 'rs2800719', 'rs2745354', 'rs2800718', 'rs2800723', 'rs2800722']
variant probabilities for this effect group: [0.0687, 0.0664, 0.0658, 0.0646, 0.0608, 0.0532, 0.0521, 0.0514, 0.0464, 0.0411, 0.041, 0.0404, 0.0398, 0.0395, 0.0387, 0.038, 0.034, 0.0328, 0.0324, 0.0323, 0.0322]
shared probability for this effect group: 0.0002
```

The results have been saved into the **colocalization summary** (cs) files with columns:
`cs` for variant representations in effect groups; 
`share` for colocalization probabilities;
`variantProb` for variant representation weights in effect groups;
`k` for the max number of causal signals specified in the model.

```
$> cat res/BMD_SH.txt_RSPO3_SH.txt.cs 
cs      share   variantProb     k
rs7741021/rs9482773     1.0     0.5936/0.4064   10
rs853974        1.0     1.0     10
rs577721086/rs72959041  0.0114  0.5313/0.4686   10
rs717796/rs7775814/rs1512450/rs7775090/rs1512449/rs910536/rs6902741/rs1569870/rs7738255/rs1080708/rs7756072     0.0001  0.1215/0.1166/0.1016/0.0878/0.0866/0.0849/0.0784/0.0746/0.0705/0.0704/0.0638      10
rs11759578/rs77525683/rs73593094/rs73593068/rs73577838/rs9482768        0.0003  0.6213/0.2032/0.0666/0.031/0.0264/0.0131        10
rs2800728/rs727330/rs727331/rs2800720/rs2800721/rs10456964/rs2800733/rs2800727/rs2745351/rs2745356/rs2800732/rs2745355/rs727332/rs2800729/rs2800730/rs719728/rs2800719/rs2745354/rs2800718/rs2800723/rs2800722    0.0002  0.0687/0.0664/0.0658/0.0646/0.0608/0.0532/0.0521/0.0514/0.0464/0.0411/0.041/0.0404/0.0398/0.0395/0.0387/0.038/0.034/0.0328/0.0324/0.0323/0.0322   10
```

## Versions
* Version 1.0.0 (2023/01/25) Initial release
* Version 2.0.0 (2023/08/25) Added K=1 in case of LD mismatch
* Version 3.0.0 (2024/02/21) Added scripts to run SharePro in one locus with GWAS (BETA/SE/N) and ld files

## Citations

If you find SharePro for colocalization analysis useful, please cite:

[Wenmin Zhang, Tianyuan Lu, Robert Sladek, Yue Li, Hamed Najafabadi, Josée Dupuis. SharePro: an accurate and efficient genetic colocalization method accounting for multiple causal signals.](https://doi.org/10.1101/2023.07.24.550431)