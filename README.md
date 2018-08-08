# clustMMRA: cluster MicroRNA Master Regulators Analysis

## Introduction
clustMMRA is a pipeline for the identification of genomically co-clustered microRNAs driving cancer subtypes.

clustMMRA is composed of 3 steps: 
1. Identification of differentially expressed miRNA clusters
2. Target enrichment analysis
3. Network analysis

**ATTENTION:** Due to the computational complexity of step3, the pipeline needs to run on a computational cluster. The scripts provided with the ARACNE distribution, are designed for the Rocks 3.2 distribution of Linux and the Sun Grid Engine (SGE) job scheduler. Thus these scripts use the ‘qsub’ command to submit jobs to the cluster, and this command must be modified based on the user’s job-scheduling software.

## Installation
The following software components are required to run clustMMRA:
* clustMMRA
* [R](http://www.r-project.org/)
* R packages: preprocessCore, plyr, Matching
* Perl
* ARACNE:
  * download file aracne.zip at [ARACNE](http://wiki.c2b2.columbia.edu/califanolab/index.php/Software/ARACNE)
  * save the file in MMRA_pipeline/codes
  * unzip aracne.zip
  * The scripts provided with the ARACNE distribution, are designed for the Rocks 3.2 distribution of Linux and the Sun Grid Engine (SGE) job scheduler. Thus these scripts use the ‘qsub’ command to submit jobs to the cluster, and this command must be modified based on the user’s job-scheduling software

In principle, clustMMRA should run under Linux and Mac.

## Usage
To run clustMMRA:
* Go to `clustMMRA_pipeline/codes`
* `./run_clustMMRA.sh`

## Inputs
While running, the program requires the following inputs:
* The path to the file containing the mRNA expression matrix in log2 with the first column containing genes names and the first row containing sample names.
* The path to the file containing the miRNA expression matrix in log2 with
the first column containing genes names and the first row containing
sample names.

**ATTENTION:** The two matrices have to be composed of the same samples. Moreover the samples have to be ordered according to the classification (from column 1 to N samples of class 1 ,from N+1 to M samples of class 2 …)

* The number of the starting and ending column of each class. If the
mRNA dataset has been previously classified according to a signature
with a classification confidence (as FDR), then give has starting and
ending column of the class only those samples that are significantly
associated to the class according to the confidence parameter.
* The signature:
      * Two files for each class: In fact for each class you have to create a file for up-regulated genes and one for down-regulated genes. The first one has to be called up(number of the class).txt and the second has to be called down(number of the class).txt. Both have to be saved in the `clustMMRA_pipeline/data/signatures/` folder which is automatically created when you install the program. Each file has to be composed of only a column reporting the list of those genes that represent the considered class. The files have to be without any description of the two columns content.

If you don’t have a published signature available for your case of study you can build it. Pay attention to test the capability of your signature to correctly classify the samples under study. Finally, the dimension of the signature influences the dimension of the final output, so pay attention to keep them as large as possible.

* The total number of classes (subtypes)
* Kolmogorov-Smirnov pvalue and log fold change thresholds for the first step.

* The pvalue threshold for the Fisher test executed In the second step of
the second phase of the pipeline. 

## Outputs
The outputs are in the folder `clustMMRA_pipeline/results` and they are:
* step1_results.txt : a report containing the results of the first step of the
pipeline. 
* step2_results.txt : a report containing the results of the second step of
the pipeline. 
* clustMMRA_output.txt : final output of the pipeline.

## Data
The TCGA and Curie datasets used in our work can be dowloaded from the publications cited in the paper.

## Contact
Please feel free to contact us at

<laura.cantini@curie.fr>

Moreover, feel free to change the code according to your needs.

For every use of the original or modified pipeline cite (https://www.biorxiv.org/content/biorxiv/early/2018/03/28/290528.full.pdf).
