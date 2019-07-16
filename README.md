# Single Cell Read Simulator

This Software simulates a single cell read simulator for one chromosom incorporating phylogenetic relationships in order to be able to simulate tumor data sets. It firstly adds random homozygous and heterozygous mutations to the given chromosom to create two starting allels. Then it creates a random  binary tree with according random mutations simulating the process of cell division. At the end it simulates the MDA amplification and illumina sequencing of each cell and each allel.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See "Usage" on how to deploy the software.

### Prerequisites

The workflow of the MDA simulation is handled with snakemake, therefore a working snakemake managment system is necessary. To create the necessary environments anaconda or miniconda is needed.
It's suggested to install anaconda or miniconda first and use it to easily download snakemake.

NOTE: To install snakemake, the miniconda or anaconda python 3 distribution is needed.

For more information on installing anaconda or miniconda for your system, visit their website:
[https://conda.io/projects/conda/en/latest/user-guide/install/index.html](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

After installing anaconda or miniconda, it can be used to install snakemake.
```
conda install -c bioconda -c conda-forge snakemake
```

For more information on installing snakemake, visit their website:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html


### Installing

A step by step manual that tells you how to get a development environment running.

Clone the git repository into the working directory and change into the directory.

```
git clone https://gitlab.ethz.ch/leonors/bachelor-thesis.git path/to/workdir
cd path/to/workdir
```

Now all the necessary files should be installed. To make sure everything works correctly and the necessary environments are installed, run the software with the given test data.

## Test

The necessary test data is given in the Data/TestReferenceData folder and can be used with the original configuration file.

### Running the test

By running the test data, you can make sure that everything was downloaded correctly and the necessary environments get installed. The test data simulates a short genomic sequence and a cell division process leading to 4 cells. You should get 8 .fq output files in the results directory (one for each allel) and the file SingleCellReadSimulator4.txt telling you 'Simulation ran succesfully.'.
There should also be a graph of the generated phylogenetic tree in the Data/Tree folder.

Run the test.
```
snakemake --use-conda "results/SingleCellReadSimulator_Test.txt"
```
To get rid of the test files quickly use the snakemake rule [cleanAll](#running-the-simulator).
More information on this in the usage section.

## Usage

All interaction with the simulation is handled over the configuration file "config.json".
The .fq files containing the chromosoms simulated by the single cell read simulator can be found in the results folder. To have a look at which allel number belongs to which cell in the phylogenetic tree, there is a .pdf file in the Data/Tree with the graph.

This is the structure of the output .fq files in the results folder:
```
simulatedAmplificationAndSequencing_<choosenName>_Allel<allelNumber>_<first/secondRead>.fq
```

### Edit config file

* To adapt the used chromosom, add your fasta file and the vcf metadataframe needed for all files being created to the "setReferenceGenomeAndMetaDataFrame" rule:

 NOTE: The contig header in the vcf files must include an ID (contig=<ID=<n>>) and this ID must be equal to the header given in the fasta file.

 ```
"referenceGenome": "path/to/your/referenceGenome.fa"
"metaDataFrame": "path/to/your/metaDataFrame.vcf"
```

* To change the number of cells (leaves of the tree) being simulated, change "generateBinaryTree":

 ```
"numberOfLeaves": chosenNumberOfCells
```

* The number of introduced germline mutations and the fraction which is heterozygous, can be changed in the rule "createGermlineMutations":

 ```
 "numberOfGermlineMutations": chosenNumberOfGermlineMutations
 "fractionOfHeterozygousGermlineMutations": chosenFractionOfHeterozygousMutations
```

* Adapting the number of mutations occuring from cell division can be done through the "createMutations" rule:

 ```
"numberOfMutations": chosenNumberOfCellDivisionMutations
```

* To change the parameters for the MDA amplification and illumina sequencing, the rule "simulateReads" can be used.
 The MDA amplification is simulated by splitting the chromosomes into different size pieces (bins) with different coverages. With some propability a mutation occurs during this MDA amplification and gets copied to a fraction of the coverage. Then the copies are sequenced and some of the reads get dropped.

 NOTE: the example parameters given below, should be a sensible start for the simulation

 IMPORTANT: Numbers are drawn form a negative binomial distributions with probability p and overdispersion r. Both are calculated from the given means and standart deviations and to not get negative numbers it must hold:
 (standartDeviation)^2 > mean

 Changing the bin size:
```
  "meanBinSize": 500
  "standartDeviationOfBinSize": 100

 ```
Changing the mean coverage and standart deviation of the bins:
  ```
  "meanCoverage": 250
  "standartDeviationOfCoverage": 16
 ```
 Changing the read lengths fed to the illumina sequencer:
  ```
  "readLength": 100
 ```
 Changing the mean fragment length and standart deviation of the illumina sequencer:
  ```
  "meanFragmentSize": 250
  "standartDeviationOfFragmentSize": 16
 ```
 Changing the probability with which an MDA amplification error occurs:
  ```
  "MDAamplificationErrorProbability": 0.0000001
 ```
 Changing the ART Illumina sequencing coverage:
 ```
 "ARTcoverage": 1
 ```
 Changing the probability with which a read is dropped during sequencing:
 ```
 "ARTDropReadProbability": 0.05
 ```

* For every rule including random number drawing there exists a seed for reproducibilaty, which can be changed to get different outcomes:
```
"seed": 479
```


### Running the Simulator

* To run the simulator use snakemake with the conda flag and give it name, to be able to distiguish different simulations:
```
snakemake --use-conda "results/SingleCellReadSimulator_<choosenName>.txt"
```
* It is also possible to let snakemake automatically determine which parts of the workflow can be run in parallel.

  To run with up to N cores:
```
snakemake -j <N> --use-conda "results/SingleCellReadSimulator_<choosenName>.txt"
```
  To run with the number of available CPU's:
```
snakemake -j --use-conda "results/SingleCellReadSimulator_<choosenName>.txt"
```

* To clean up all files and folders created by the simulation use:
```
snakemake cleanAll
```

For further informations on command line options for the execution with snakemake call:
```
snakemake -h
```

## Built With

* [python 3](https://www.python.org/doc/) - Used to generate the scripts
* [numpy](https://www.numpy.org/devdocs/) - Used for data storage in the scripts
* [pyvcf](https://pyvcf.readthedocs.io/en/latest/index.html) - Used to handle vcf files
* [miniconda](https://docs.conda.io/en/latest/miniconda.html#) - Used to compile the environments
* [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) - The workflow framework used
* [biopython](https://biopython.org/wiki/Documentation) - Used handle fasta and fastq files
* [samtools](http://www.htslib.org/doc/#manual-pages) - Used to handle file transformations
* [art](http://www.dropwizard.io/1.0.2/docs/) - Used to simulate illumina sequencing
* [bcftools](http://www.htslib.org/doc/#manual-pages) - Used to handle file transformations
* [htslib](http://www.htslib.org/doc/#manual-pages) - Used to handle file transformations
* [graphiz](https://www.graphviz.org/documentation/) - Used to illustrate the phylogenetic tree

## Contributing

Pull requests are welcome.

## Authors

* **Schubert Leonor** - *Initial work*

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

* Beerenwinkel, Niko, Prof. Dr.
* Singer Jochen
* Marass Francesco, Dr.
