# FAnnP
Functional Annotation Pipeline

**Current version:** 1.1

The functional annotation pipeline is mainly designed for the annotation of prokaryotic MAGs (Metagenome Assembled Genomes), although it can also be used for functional annotation of any other type of assembly (genomic or metagenomic).

The workflow implemented in this pipeline will process all the sequences in one or more fasta files in order to first identify coding sequences, and then annotate its possible functional role, according to different databases and diverse homology searching methods.

This workflow was designed based on the pipeline implemented for the analysis from (Dombrowski et al., 2020) 1.

More detailed description of the steps and databases can be accessed via the following [repository](https://zenodo.org/record/3839790#.X2r6VLJR2lN)

## Quick start

**Required input files**
Directory with all the genomes, bins or assemblies to annotate.*

_*All of them need to have the same extension (symbolic links are allowed)_

**Download or clone the repository**

> git clone https://github.com/AlejandroAb/FAnnP



**Edit configuration file**
Next, edit the configuration file in order to enter the directory with the sequences to annotate and select your target databases.

Input files and Run name
First, enter the name of your RUN, then the FULL path to your directory containing the bins to annotate, and the bin extension.

```yaml
RUN: "test_run1" 
bin_dir: "/export/lv3/scratch/projects_XX/MyBins/" 
bin_list: ["bin1","bin2","bin3","bin4",...,"binN"]
bin_ext: "fna" 
```

Sometimes the bin_list can contain several bins to be annotated, in this sense, it can be a fastidious task to manually enter all this information. For these cases, we supply a script in order to enter all the bins/genomes given a directory (bin_dir), the extension of the files in such directory (bin_ext) the target configuration file and optionally the name of a new configuration file. If the last argument (new name) is not supplied the bin_list will be replaced on the target config file:

> Scripts/list_Bins.sh test_dir/ fasta config.yaml config_new.yaml

## Run the pipeline
Once that you have all the files in place and your configuration file done:

snakemake --configfile config.yaml -np
