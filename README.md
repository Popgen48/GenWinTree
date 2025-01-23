# Genomic Window Tree Workflow

This repository contains a Snakemake/Python workflow developed to generate phylogenetic trees from genomic data. The workflow processes VCF files, segments chromosomes into windows, and builds phylogenetic trees using [RAxML](https://github.com/amkozlov/raxml-ng).

The first step of the workflow is to produce FASTA alignments for all combinations of the specified taxa, followed by deriving maximum likelihood trees with the RAxML software. Alignments are created directly from the VCF files through the following three steps:
1. Heterozygous SNPs are converted to their corresponding IUPAC codes.
2. Windows with a genotyping rate of less than 95% are excluded.
3. Windows with fewer than 11 total SNPs are also excluded.

## Input Files
#### 1. VCF Files:
The files must be organized by chromosome, with each chromosome having its own VCF file. The file names should align with the pattern specified in the configuration file (vcfPattern). Examples: chromosome1_filtered.vcf.gz, chromosome2_filtered.vcf.gz., etc.
####  2. Configuration File:
A JSON file named _"config.json"_ is required. This file will contain various parameters and tool paths necessary for the workflow. Below is an example of the configuration:

```{
    "tools": {
        "RaxML":"/path/to/raxml-ng"
    },
    "base_dir": {
        "inputDir":"/path/to/input/directory",
        "outputDir": "/path/to/output/directory"
    },
    "params": {
        "sampleList": ["sample1", "sample2", "sample3", "sample4"],
        "windowSize": "10000",
        "stepSize": "50000",
        "threads": "2",
        "vcfPattern": "_filtered.vcf.gz"
    }
}
``` 


|Parameter|Explanation|
|------|------|
| tools.RaxML : | Path to the RAxML-NG executable. |
| base_dir.inputDir : | Directory containing the input VCF files. |
| base_dir.outputDir : | Directory where output files will be saved. |
| params.sampleList : | List of samples to be used. |
| params.windowSize : | The size of each window (in base pairs). |
| params.stepSize : | The step size between windows (in base pairs). |
| params.vcfPattern : | Pattern used to match the VCF files. |

## Output Files
The workflow generates a folder for each chromosome that contains the following files:

 1. General RAxML output files.
 2. A summary file named _"{chrom}.tree.out"_, which contains information for each genomic window. This file includes genotype frequencies, depth statistics, and the corresponding phylogenetic tree. Each line in the file represents a genomic window and includes the following fields::

    - Chromosome Name: Extracted from the header of the VCF file.
    - Window Start: The starting position (in base pairs) of the genomic window on the chromosome.
    - Window End: The ending position (in base pairs) of the genomic window. If the window extends beyond the end of the chromosome, this will correspond to the chromosome length.
    - Total Sites: The number of sites in the VCF file that fall within the genomic window.
    - Genotype Frequencies: A comma-separated list representing the proportion of genotyped sites for each sample in the window.
    - Average Depth: A comma-separated list representing the average sequencing depth for each sample in the window.
    - Phylogenetic Tree: The phylogenetic tree for the window in Newick format.

Additionally, a log file named _"workflow.log"_ is created. This log file records the progress of the workflow and any issues encountered during execution.
 
## Requirements
Make sure the following tools are installed on your system:
- Python: Version 3.x
- Snakemake: For workflow orchestration
- RAxML-NG: For building phylogenetic trees
- pysam: Python module for interacting with VCF files
- glob: For file path manipulation
- json: For reading configuration files

## Installation and usage

To use this workflow, clone the repository to your local machine:

```
git clone clone https://github.com/Popgen48/GenWinTree.git
cd GenWinTree
```

To run the workflow, you need to specify the configuration file (config.json) and ensure that the necessary input files (VCF files) are available in the specified input directory.
Command:
```
snakemake --snakefile GenWinTree.py --cores <number_of_threads>
```
## Credits

The workflow GenWinTree.py was originally written by @BioInf2305 with contributions from @NPogo.

## Citation

## License
MIT
