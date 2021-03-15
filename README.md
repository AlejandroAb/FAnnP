# FAnnP
Functional Annotation Pipeline

**Current version:** 1.1

The functional annotation pipeline is mainly designed for the annotation of prokaryotic MAGs (Metagenome Assembled Genomes), although it can also be used for functional annotation of any other type of assembly (genomic or metagenomic).

The workflow implemented in this pipeline will process all the sequences in one or more fasta files in order to first identify coding sequences, and then annotate its possible functional role, according to different databases and diverse homology searching methods.

This workflow was designed based on the pipeline implemented for the analysis from ([Dombrowski et al., 2020](https://www.nature.com/articles/s41467-020-17408-w)).

More detailed description of the steps and databases can be accessed via the following [repository](https://zenodo.org/record/3839790#.X2r6VLJR2lN)

## Quick start

**Required input files**

Directory with all the genomes, bins or assemblies to annotate.*

_*All of them need to have the same extension (symbolic links are allowed)_

**Download or clone the repository**

> git clone https://github.com/AlejandroAb/FAnnP



**Edit configuration file**

Next, edit the configuration file (config.yaml) in order to enter the directory with the sequences to annotate and select your target databases.

Input files and Run name
First, enter the name of your RUN, then the FULL path to your directory containing the bins to annotate, and the bin extension.

```yaml
RUN: "test_run1" 
bin_dir: "/export/lv3/scratch/projects_XX/MyBins/" 
bin_list: ["bin1","bin2","bin3","bin4",...,"binN"]
bin_ext: "fna" 
```

Sometimes the bin_list can contain several bins to be annotated, in this sense, it can be a fastidious task to manually enter all this information. For these cases, we supply a script in order to enter all the bins/genomes given a directory (bin_dir), the extension of the files in such directory (bin_ext) the target configuration file and optionally the name of a new configuration file. If the last argument (new name) is not supplied the bin_list will be replaced on the target config file:

`Scripts/list_Bins.sh test_dir/ fasta config.yaml config_new.yaml`

## Run the pipeline
Once that you have all the files in place and your configuration file done:

`snakemake --configfile config.yaml`

*Gene calling*

You can choose between Prokka and Prodigal. We recommend the first one for small to medium datasets; Prodigal for bigger datasets.

```yaml
GENE_CALLING: "PROKKA" #PRODIGAL OR PROKKA

prodigal:
 procedure: "meta" #Select procedure (single or meta)
 extra_params: "" 

prokka:
 cpus: 15
 extra_params: " --mincontiglen 500" 
 increment: "1" 
Map genome/bin taxonomy
```

You can provide a file with taxonomy annotation from your MAGS or your genomes and this information will end up in the summary file and can be used later on for further cleaning steps.

```yaml
bin_taxonomy:
 add_taxonomy: "T"  # Add taxonomy T or F
 taxonomy_file: ""  # Full path to the tab separated file with the bin/genome name (same as the one supplied at the bin_list) and its taxonomy
 taxonomy_sep: ";"  # Character that separates the different taxonomy levels.
 bin_col: "1"       # Column whith the bin/genome name.  [Def. "1" (Metagenomic pipeline summary file)]
 taxonomy_col: "17" # Column with the taxonomy. [Def. "17" (Metagenomic pipeline summary file)]
```

*Bin cleaning*

It is possible for the users to clean their bins/MAGs (Metagenomic Assembled Genomes) by inspecting different parameters within the contigs of each single bin. This by applying different criteria, that can be fine-tuned according to the user needs.

* *1. GC percentage.* The user can set the maximum tolerated difference between the GC% of the bin vs the GC% of the contig. If this difference is greater than the user threshold (10% by default) the complete contig will be discarded.
* *2. Coverage ratio.* The coverage ratio is obtained by dividing the contig average depth / Bin average depth if this ratio is grater or lower than the user threshold (0.5 by default) the contig is discarded.
* *3. Contig domain identity.* If the taxonomy for the bins/MAGs is supplied, this filter step will compute the percentage of proteins within the same contig belonging to the same taxonomy domain as the bin. If this percentage is sufficient enough the contig is retained, otherwise will be discarded. By default 50%.

```yaml
bin_cleaning:
 clean_bins: "T"          # Clean the bins, T or F
 bin_coverage_gc_file: "" # Full path to the file with GC% and coverage information (Metagenomic pipeline summary file)
 bin_col: "1"             # Column whith the bin/genome name.  [Def. "1" (Metagenomic pipeline summary file)]
 avg_depth_col: "23"      # Column whith the coverage depth data. [Def. "23" (Metagenomic pipeline summary file)]
 avg_gc_col: "24"         # Column whith the GC% data.  [Def. "24" (Metagenomic pipeline summary file)]
 contig_coverage_file: "" # Full path to the file with coverage information by contig. Output file from Metagenomic pipline: /binning/FinalBins/contig_coverage.txt
 rules:
  GC_max_diff: 0.1        # Maximum tolerated GC% difference between bin and contig. Values between 0 and 1. [Def. 0.1 (10%)]
  Coverage_ratio: 0.5     # Maximum tolerated ratio difference between bin and contig. Having:  (Contig_avg_depth/Bin_avg_depth = val) ->  (1-Coverage_ratio) > val < (1+Coverage_ratio). Values between 0 and 1.
  Contig_domain_identity: 50 # Minimun percentage of proteins from a contig assigned to the same domain as the bin/genome. Value between 0 and 100. [Def. 50]
```

*Custom your annoation databases*

The name of the following sections/options is equivalent to the name of the databases or feature to be annotated. In general, the name is followed by the tool that is used for making such annotation (diamond, blast, hmmr). Set the option annotate: "T" or "F" in order to turn on / off such database annotation. The concatenate_bins flag if set to "T" it will perform the annotation for all the proteins from all the bins/genomes at once (concatenating all of them in a single file); otherwise, it will perform this operation individually for all the bins/genomes. This last scenario is more useful when you have a lot of computational power, so then, you can parallelize more jobs, each one with its own number of CPUs.


## Run the pipeline
Once that you have all the files in place and your configuration file done:

`snakemake --configfile config.yaml`

## Output files structure

```
{RUN}
├── FunctionalAnnotation.tsv   Summary file with the annotations by protein
├── clean_bins                 Directory with bins after cleaning steps
│   ├── clean.log
│   ├── clean_bin1.fna
...  ...
    └── clean_binN.fna
├── prokka                       Directory with protein predictions
│   ├── All_bins_clean.faa
│   ├── concatenate_bins.benchmark
│   ├── prodigal.benchmark
│   ├── renamed
│   ├── bin_1                     Directory with protein predictions per bin
│   └── bin_2
├── arCOG                        Results per database
│   └── arcogs_map.tsv
├── cog_diamond
│   └── cog_map.tsv
├── contig_mapping
│   └── Mapping files
├── diamond
│   └──  diamond_map.tsv
├── kfam
│   └── kfam_map.tsv
├── pfam
│   └── pfam_map.tsv
....
└── <DB>
    └── <DB>_map.tsv
```
