"""
FunctionalAnnotation Pipeline
Version: 2.0
Author: Alejandro Abdala, Julia Engelmann, Nina Dombrowski, and Tijn Tensen
Last update: 03/03/2025
"""
from random import randint
import datetime
import sys


"""
Get the run name and create the run ID.
"""
run=config["RUN"]
def run_id_generator(file):
    name_list = file.split("/")
    return f"{name_list[-1][0:3].lower()}{''.join([str(randint(0,9)) for i in range(5)])}"
run_id=run_id_generator(config["RUN"])
args = sys.argv
try:
    config_path = args[args.index("--configfiles") + 1]
except:
    config_path = args[args.index("--configfile") + 1]

"""
Get the current date time for the start of the pipeline.
"""
def log_date_time():
    return datetime.datetime.now()
date_time=log_date_time()


"""
All the other rules will need to run because of this rule. 
This rule will execute when the pipeline is being run without specifying a rule.
"""
rule all:
    input:
        #"{run}/{bin}/report.txt"
        #dynamic("{run}/{bin}/report.txt")        
        #expand("{run}/{bin}/report.txt" if config["concatenate_bins"] == "T" else "{run}/report.txt",bin=config["bin_list"], run=run)
        expand("{run}/report.txt",bin=config["bin_list"], run=run) #if config["concatenate_bins"] == "T" else
        #expand("{run}/{bin}/report.txt",bin=config["bin_list"], run=run) 
        #"{run}/report.txt"

if config["GENE_CALLING"] == "PRODIGAL":
    """
    Runs Prodigal to extract all the genes from each bin out of the contigs.
    The faa file will store the protein sequences from the genes.
    The gbk file will store additional information about the genes.
    The ffn file will store the nucleotide sequences from the genes.
    """
    rule prodigal_bins:
        input:
            config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            gbk=temp("{run}/gene_calling/{bin}/{bin}.gbk"),
            genes="{run}/gene_calling/{bin}/{bin}.ffn",
            prots="{run}/gene_calling/{bin}/{bin}.faa"
        benchmark:
            "{run}/benchmark/prodigal_bins.{bin}.benchmark"
        shell:
            "prodigal -a {output.prots} -d {output.genes} -f gbk  -i {input} "
            "-p "+ config["prodigal"]["procedure"]+"  "+ config["prodigal"]["extra_params"]+"  -o {output.gbk} "
elif config["GENE_CALLING"] == "PROKKA":
    """
    Runs Prokka to extract all the genes from each bin out of the contigs.
    The faa file will store the protein sequences from the genes.
    The gbk file will store additional information about the genes.
    """
    rule prokka_bins:
        input:
            config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            "{run}/gene_calling/{bin}/{bin}.faa",
            "{run}/gene_calling/{bin}/{bin}.gbk"
        params:
            output_dir="{run}/gene_calling/{bin}",
            file_ext=config["bin_ext"],
            bin_prefix="{bin}"
        benchmark:
            "{run}/benchmark/prokka_bins.{bin}.benchmark"
        threads:
            int(config["prokka"]["cpus"])
        shell:
            "prokka --prefix {wildcards.bin} --outdir {params.output_dir}  --addgenes --force --increment  "+str(config["prokka"]["increment"]) + " "
            " --compliant --centre UU --cpus "+ str(config["prokka"]["cpus"]) + " --norrna --notrna "+ config["prokka"]["extra_params"]+" {input}"
else:
    """
    Skips gene calling if Prodigal or Prokka is enabled. Requires faa protein files as input.
    """
    rule skip_bins:
        input:
            config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            "{run}/gene_calling/{bin}/{bin}.faa"
        benchmark:
            "{run}/benchmark/skip_bins.{bin}.benchmark"
        shell:
            "cat {input} > {output}"
            
"""
Always uses prodigal to acquire the gff files.
"""
rule prodigal_for_gff:
    input:
        bins=config["bin_dir"]+"{bin}."+config["bin_ext"]
    output:
        gff=temp("{run}/gene_calling/{bin}/{bin}.gff")
        #genes="{run}/gene_calling/{bin}/{bin}.ffn",
        #prots="{run}/gene_calling/{bin}/{bin}.faa"
    benchmark:
        "{run}/benchmark/prodigal_for_gff.{bin}.benchmark"
    shell:
        "prodigal -f gff  -i {input.bins} "
        "-p "+ config["prodigal"]["procedure"]+"  "+ config["prodigal"]["extra_params"]+"  -o {output.gff} "

"""
Pastes all the bins contigs behind the gff files.
"""
rule add_seq_gff:
    input:
        gff="{run}/gene_calling/{bin}/{bin}.gff",
        genome=config["bin_dir"]+"{bin}."+config["bin_ext"]
    output:
        gff="{run}/gene_calling/{bin}/{bin}.gff3"
    benchmark:
        "{run}/benchmark/add_seq_gff.{bin}.benchmark"
    shell:
        "Scripts/concat_gff_seqs.sh {input.gff} {input.genome} {output}"

"""
Checks if the setting to add the bins contigs to the gff files is enables.
"""
rule gff_check:
    input:
        expand("{run}/gene_calling/{bin}/{bin}.gff", bin=config["bin_list"], run=run)
        if config["ADD_SEQ_TO_GFF"].lower() == "f" else
        expand("{run}/gene_calling/{bin}/{bin}.gff3", bin=config["bin_list"], run=run) 
    output:
        "{run}/contig_mapping/Gff_check.txt"
    benchmark:
        "{run}/benchmark/gff_check.benchmark"
    shell:
        "echo {input} > {output}"

"""
Reconstructs the gene names by pasting the bin names before the gene names if config["rename_bins"] == "T"
"""
rule rename_bins:
    input:
        "{run}/gene_calling/{bin}/{bin}.faa"
    output:
        "{run}/gene_calling/renamed/{bin}.faa"
    benchmark:
        "{run}/benchmark/rename_bins.{bin}.benchmark"
    shell:
        "awk -v run=\"{wildcards.run}\" -v bin=\"{wildcards.bin}\" '/>/{{sub(\">\",\"&\"FILENAME\"|\");sub(/\.faa/,x);sub(run\"/gene_calling/\"bin\"/\",x)}}1' {input} | "
        "cut -f1 -d \" \" > {output}" if config["rename_bins"] == "T" else "ln -s ../{wildcards.bin}/{wildcards.bin}.faa {output}"

"""
Checks if the gene names are bigger than 100 characters. If this is the case we won't be able to store it in the FAnnP database.
Skips the rule if config["generate_sql_database"] == "F"
The output file is necessary for snakemake and only exists for the workflow.
"""
rule check_bins_length:
    input:
        expand("{run}/gene_calling/renamed/{bin}.faa", bin=config["bin_list"], run=run)
    output:
        temp("{run}/gene_calling/renamed_checked_log.txt")
    benchmark:
        "{run}/benchmark/check_bins_length.benchmark"
    shell:
        "bash Scripts/gene_length_checker.sh '{input}' {output}" if config["generate_sql_database"] == "T" else "echo 'Skip length check.' > {output}"

"""
Creates a simple list with all the gene names and their respective bin name.
"""
rule make_protein_bin_list:
    input:
        expand("{run}/gene_calling/renamed/{bin}.faa", bin=config["bin_list"], run=run)
    output:
        "{run}/contig_mapping/1_Bins_to_protein_list.txt"
    benchmark:
        "{run}/benchmark/make_protein_bin_list.benchmark"
    shell:
        "cat {input} | grep '^>' | cut -f1 -d \" \" | sed 's/>//g' | awk -F'\\t' -v OFS='\\t' '{{split($1,a,\"|\"); print $1, a[1]}}' | LC_ALL=C sort > {output}"

"""
Gzips genebank files.
"""
rule gzip_gbk:
    input:
        "{run}/gene_calling/{bin}/{bin}.gbk"
    output:
        "{run}/gene_calling/{bin}/{bin}.gbk.gz"
    benchmark:
        "{run}/benchmark/gzip_gbk.{bin}.benchmark"
    shell:
        "gzip {input}"

if config["GENE_CALLING"] == "PRODIGAL":
    """
    Will run with PRODIGAL setting. Creates the contig protein list using the gene bank files using the next three rules.
    """
    rule seqid_to_prot:
        input:
            "{run}/gene_calling/{bin}/{bin}.gbk.gz"
        output:
            "{run}/gene_calling/{bin}/{bin}.seq2prots"
        benchmark:
            "{run}/benchmark/seqid_to_prot.{bin}.benchmark"
        shell:
            "zcat {input} | grep -A1 --no-group-separator \"CDS\" | "
            "grep \"/note\" | sed 's/.*\\/note=\"//' | "
            "cut -f1 -d\";\" | sed 's/ID=// ; s/_/\\t/' > {output}"
    """
    Create the contig protein lists for Prodigal.
    """
    rule parse_gbk_files_prodigal:
        input:
            gbk="{run}/gene_calling/{bin}/{bin}.gbk.gz",
            map="{run}/gene_calling/{bin}/{bin}.seq2prots"
        output:
            temp("{run}/gene_calling/{bin}/Contig_Protein_list.txt")
        params:
            b="{bin}"
        benchmark:
            "{run}/benchmark/parse_gbk_files_prodigal.{bin}.benchmark"
        shell:
            "zcat {input.gbk} | grep DEFINITION | cut -f3 -d\" \" | "
            "cut -f1-3 -d\";\" | sed 's/;/\\t/g ; s/seq.[a-z]*=//g ; s/\"//g' | "
            "awk -F\"\\t\" -v bin=\"{params}\" 'NR==FNR{{contig[$1]=$3;len[$1]=$2;next}} BEGIN{{OFS=\"\\t\"}} {{print $1,bin,contig[$1],len[$1],contig[$1]\"_\"$2,\"NA\"}}' "
            "- {input.map} > {output}"
    """
    Puts the individual files into one file.
    """
    rule concat_gbk_files_prodigal:
        input:
            expand("{run}/gene_calling/{bin}/Contig_Protein_list.txt",bin=config["bin_list"], run=run)
        output:
            "{run}/contig_mapping/Contig_Protein_list.txt"
        benchmark:
            "{run}/benchmark/concat_gbk_files_prodigal.benchmark"
        shell:
            "cat {input} > {output}"
elif config["GENE_CALLING"] == "PROKKA":
    """
    Will run with PROKKA setting. Creates the contig protein list using the gene bank files using the next two rules.
    """
    rule create_list_gbk:
        input:
            expand("{run}/gene_calling/{bin}/{bin}.gbk.gz", bin=config["bin_list"], run=run)
        output:
            "{run}/gene_calling/list_of_gzips"
        benchmark:
            "{run}/benchmark/create_list_gbk.benchmark"
        shell:
            "echo {input} | sed 's/ /\\n/g' > {run}/gene_calling/list_of_gzips"
#it seems that the previous command generates on extra line which 
#ends up in a error.
    """
    Create the contig protein list for Prokka.
    """
    rule parse_gbk_files:
        input:
            "{run}/gene_calling/list_of_gzips"
        output:
            "{run}/contig_mapping/Contig_Protein_list.txt"
        benchmark:
            "{run}/benchmark/parse_gbk_files.benchmark"
        shell:
            "python Scripts/Parse_prokka_for_MAGs_from_gbk-file.py -i {input} -t {output}"
else:
    """
    GENE_CALLING: "SKIP". Will run with SKIP setting. This will format the contig protein list correctly for the Diamond, HMMER and Blast rules.
    """
    rule parse_faa_files_skip_gbk:
        input:
            expand("{run}/gene_calling/renamed/{bin}.faa", bin=config["bin_list"], run=run)
        output:
            "{run}/contig_mapping/Contig_Protein_list.txt"
        benchmark:
            "{run}/benchmark/parse_faa_files_skip_gbk.benchmark"
        shell:
            "bash Scripts/skip_gene_calling_contig_parser.sh '{input}' {output}"

"""
Replaces the protein ID column with the contig ID.
"""
rule protein_to_binID:
    input:
        "{run}/contig_mapping/1_Bins_to_protein_list.txt"
    output:
        "{run}/contig_mapping/Prokka_to_BinID.txt"
    benchmark:
        "{run}/benchmark/protein_to_binID.benchmark"
    shell:
        "cat {input} | awk 'BEGIN{{OFS=\"\\t\"}}{{split($1,a,\"|\"); print a[2],$2}}' | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{parts=split($1,a,\"_\"); contig_id=a[1]; for(i=2;i<parts;i++){{contig_id=contig_id\"_\"a[i]}}; print contig_id,$2}}' | sort | uniq > {output}"

"""
Tries to link the old and new bin IDs. There might not be a difference between them though.
"""
rule link_old_to_new_binIDS:
    input:
        prokka="{run}/contig_mapping/Prokka_to_BinID.txt",
        prot="{run}/contig_mapping/Contig_Protein_list.txt",
        gff="{run}/contig_mapping/Gff_check.txt"
    output:
        "{run}/contig_mapping/temp_Bin_Contig_Protein_list.txt"
    params:
        idx="1" if config["GENE_CALLING"] == "PROKKA" else "2"
    benchmark:
        "{run}/benchmark/link_old_to_new_binIDS.benchmark"
    shell:
        "awk -v idx={params.idx} 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$idx]}}' {input.prokka} {input.prot} | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,$7,$2,$2,$3,$4,$5}}'  > {output}"
        if config["GENE_CALLING"] == "PROKKA" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,$2,$3,$3,$4,$5,$6}}'  {input.prot} > {output}" 

"""
Add in an extra column for contig number.
"""
rule add_contig_nr:
    input:
        "{run}/contig_mapping/temp_Bin_Contig_Protein_list.txt"
    output:
        "{run}/contig_mapping/Bin_Contig_Protein_list.txt"
    benchmark:
        "{run}/benchmark/add_contig_nr.benchmark"
    shell:
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{split($4,a,\"_\"); print $2\"_contig_\"a[2],$1,$2,$3,a[2],$5,$6,$7}}' "
        "{input} > {output}"  if config["GENE_CALLING"] == "PROKKA" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{ print $2\"_contig_\"$1,$1,$2,$3,$1,$5,$6,$7}}' "
        "{input} > {output}" 

"""
Merges contig and protein information with eachother.
"""
rule merge_contid_ids:
    input:
        contig_prot="{run}/contig_mapping/Bin_Contig_Protein_list.txt",
        contig_old="{run}/contig_mapping/Contig_Old_mapping_for_merging.txt"
    output:
        "{run}/contig_mapping/Bin_Contig_Protein_list_merged.txt"
    shell:#headers: accession, BinID, newContigID, oldContigID, mergeContigID, ContigLengthNew, LengthContigOld, GC, ProteinID, prokka
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.contig_old} {input.contig_prot} | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{split($7, gene_name, \" \");print gene_name[1],$3,$4,$10,$1,$6,$14,$13,$7,$8}}' > {output}" if config["GENE_CALLING"] == "SKIP" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.contig_old} {input.contig_prot} | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $3\"|\"$7,$3,$4,$10,$1,$6,$14,$13,$7,$8}}' > {output}"

"""
Extracts the GC percentage and the length for each contig.
"""
rule computeGC_bins:
    input:
        config["bin_dir"]+"{bin}."+config["bin_ext"]
    output:
        temp("{run}/contig_mapping/{bin}_temp1")
    benchmark:
        "{run}/benchmark/computeGC_bins.{bin}.benchmark"
    shell:
        "perl Scripts/length+GC.pl {input} > {output}"

"""
Adds the bin ID as a column and simplifies the contig ID.
"""
rule compute_num_contigs:
    input:
        "{run}/contig_mapping/{bin}_temp1"
    output:
        temp("{run}/contig_mapping/{bin}_temp2")
        #OUTPUT: Original_ID      BIN    Num_consecutive %GC     Length 
        #Prodigal: k141_10861      low_completion-refined_85       1       0.467   751
        #Prokka:   NODE_12_length_19514_cov_795.143746     unbinned        1       0.272   19514
    params:
        sample="{bin}"
    benchmark:
        "{run}/benchmark/compute_num_contigs.{bin}.benchmark"
    shell:
        "cat {input} | awk -v bin=\"{params.sample}\"  -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}} {{split($1,a, \" \"); print a[1],bin,FNR,$2,$3}}' > {output} "
        #Output:k141_10861      low_completion-refined_85       1       0.467   751
        #add num contigs
        #"cat {input} | awk  '$1=(FNR FS $1 FILENAME)' | " 
        #add in binIDs
        #"awk 'BEGIN{{OFS=\"\\t\"}}{{split($2,a,\"-\"); print a[2],$1,$3,$4}}' | "
        #
        #"awk 'BEGIN{{OFS=\"\\t\"}}{{split($1,a,\"/\"); print a[1],a[2], $2,$3,$4}}' | "
        #
        #"sed 's/contig_maping//g' | sed 's/_temp1//g' > {output}"

"""
Puts all the contig tables in one file.
"""
rule concatenate_contig_data:
    input:
        expand("{run}/contig_mapping/{bin}_temp2", bin=config["bin_list"], run=run)
    output:
        "{run}/contig_mapping/Contig_Old_mapping.txt"
    benchmark:
        "{run}/benchmark/concatenate_contig_data.benchmark"
    shell:
        "cat {input} > {output}" 

"""
Formats the columns so it can be merged with another file easier.
Also creates a new column which connects the bin ID and contig number with eachother.
"""
rule create_merging_file:
    input:
        "{run}/contig_mapping/Contig_Old_mapping.txt"
    output:
        "{run}/contig_mapping/Contig_Old_mapping_for_merging.txt"
    benchmark:
        "{run}/benchmark/create_merging_file.benchmark"
    shell:
        "cat {input} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $2\"_contig_\"$3,$0}}' > {output}"

"""
Puts all gene proteins in one file.
"""
rule concatenate_bins:
    input:
        expand("{run}/gene_calling/renamed/{bin}.faa",  bin=config["bin_list"], run=run)
    output:
        temp("{run}/gene_calling/All_bins.faa")
    benchmark:
        "{run}/benchmark/concatenate_bins.benchmark"
    shell:
        "cat {input} > {output}"

"""
Checks if there are any genes present in the dataset that share the same name.
If this is the case the pipeline will crash on purpose so that the researcher can resolve this issue since gene names must be unique.
The output file is necessary for snakemake and only exists for the workflow. 
"""
rule check_name_duplications:
    input:
        "{run}/gene_calling/All_bins.faa"
    output:
        temp("{run}/gene_calling/duplications_checked_log.txt")
    params:
        rename=config["rename_bins"]
    benchmark:
        "{run}/benchmark/check_name_duplications.benchmark"
    shell:
        "bash Scripts/gene_duplicate_checker.sh {input} {params.rename} {output}"

"""
Creates the template file named B_GenomeInfo.txt for the FunctionalAnnotation.tsv file using several other files.
For every gene, information and annotations will be merged into the B_GenomeInfo.txt file using annotations found by the other tools.
This will create the FunctionalAnnotation.tsv file.
If Meta-Cascabel files are added to the config file this script will also try including this information in the output file.
"""
rule extract_all_bin_contig_protein_info:
    input:
        bins="{run}/gene_calling/All_bins.faa",
        merged="{run}/contig_mapping/Bin_Contig_Protein_list_merged.txt"
    output:
        "{run}/contig_mapping/B_GenomeInfo.txt"
    benchmark:
        "{run}/benchmark/extract_all_bin_contig_protein_info.benchmark"
    shell:
        "bash Scripts/contig_parser_7.sh {input.bins} {input.merged} {output} "+ config["meta_cascabel"]["annotate"] +" "+ config["meta_cascabel"]["location_final_bins"] +" "+ config["meta_cascabel"]["contig_coverage"]

"""
If we should clean the bins if that isn't done already.
Config files might need to be added from Meta-Cascabel to run this properly.
"""
if config["bin_cleaning"]["clean_bins"] == "T":
    """
    Filters out contigs based on if there is too much of a difference between the bin.
    It filters based on GC percentage difference and coverage ratio difference between contig and bin.
    """
    rule mark_NOK_gc_coverage_bins:
        input:
            "{run}/contig_mapping/B_GenomeInfo.txt"
        output:
            "{run}/contig_mapping/bin_clean_1"
        benchmark:
            "{run}/benchmark/mark_NOK_gc_coverage_bins.benchmark"
        shell:
            "cat {input} | tail -n +2 | awk -F \"\\t\" 'function abs(v) {{return v < 0 ? -v : v}} "
            "BEGIN{{OFS=\"\\t\"}} {{pass_gc=\"PASS\";pass_cvg=\"PASS\";if($11<=1){{gc1=$11*100}}else{{gc1=$11}}; if($5<=1){{gc2=$5*100}}else{{gc2=$5}}; gc_diff=abs(gc1-gc2); "  
            "if(gc_diff>"+str(config["bin_cleaning"]["rules"]["GC_max_diff"])+"){{pass_gc=\"FAIL\"}}; "
            "ab_ratio=$12/$4;if(ab_ratio<1-(" +str(config["bin_cleaning"]["rules"]["Coverage_ratio"])+ ") || "
            "ab_ratio>1+("+str(config["bin_cleaning"]["rules"]["Coverage_ratio"])+")){{pass_cvg=\"FAIL\"}}; "
            "if(pass_cvg==\"FAIL\" || pass_gc==\"FAIL\"){{print $2,$6,$7,gc_diff,pass_gc,ab_ratio,pass_cvg}}}}' |"
            "sort | uniq > {output}"

    """
    This rule creates a file with the following columns:
    Bin ID, Contig ID, Bin Taxonomy, Prot Taxonomy, Num Prots with that Taxonomy, Percentage, Num Prots in the contig.
    """
    rule mark_NOK_domain_bins_s1:  
        input:
            "{run}/FunctionalAnnotation.tsv"
        output:
            "{run}/contig_mapping/bin_clean_tax1"
        benchmark:
            "{run}/benchmark/mark_NOK_domain_bins_s1.benchmark"
        shell:
            "cat {input} | sed 's/d__//g' |  awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\"}} "
            "{{if(NR==1){{for(i=1;i<=NF;i++){{if($i==\"TaxString\"){{blastTax=i; }} }}" #Identify column with diamond/blast taxonomy, this can be variable
            "}}else{{split($3,bin_tax,\";\");split($blastTax,prot_tax,\",\");" #split taxonomies 
            "if(!bin[$2]){{bin[$2]=$2}};" #select uniq bins
            "if(!contig[$7]){{contig[$7]=bin_tax[1];bin_contig[$2][$7]=$7}};"#select uniq contigs and link them to its taxonomy and to its bin
            "bin_prot[$7][prot_tax[1]]++;contig_num[$7]++}};}}END{{" #count prot domains per contigs
            "for(b in bin){{" #flush the arrays with the summarized info - first iterate over bins
            "   for(c in bin_contig[b]){{" #iterate over contigs per bin
            "       for(d in bin_prot[c]){{" #iterate over different domains founded by contig
            "print b\"\\t\"c\"\\t\"contig[c]\"\\t\"d\"\\t\"bin_prot[c][d]\"\\t\"(bin_prot[c][d]/contig_num[c])*100\"\\t\"contig_num[c]" # print information
            "}} }} }} }}' | sort -k2,2 -k6g > {output}" #sort by contig and % then output

    """
    This rule take the ouput from previous rule (mark_NOK_domain_bins_s1).
    And then just compare the domain of the bin vs the domains from the contigs.
    If the prot domain is equal or above some threshold, the contig is retained.
    """
    rule mark_NOK_domain_bins_s2:
        input:
            "{run}/contig_mapping/bin_clean_tax1"
        output:
            temp("{run}/contig_mapping/bin_clean_tax2")
        benchmark:
            "{run}/benchmark/mark_NOK_domain_bins_s2.benchmark"
        shell:
            "cat {input} | awk -F \"\\t\" '{{if($3==$4 && $6>="+str(config["bin_cleaning"]["rules"]["Contig_domain_identity"])+" || $3==\"-\")"
            "{{print $1\"\\t\"$2}}}}' | sort | uniq > {output}"

    """
    Finally, we took th contigs passing the tax filter to subset the contigs failing.
    """
    rule mark_NOK_domain_bins_s3:
        input:
            all_contigs="{run}/contig_mapping/bin_clean_tax1",
            passing_contigs="{run}/contig_mapping/bin_clean_tax2"
        output:
            "{run}/contig_mapping/bin_clean_tax3"
        benchmark:
            "{run}/benchmark/mark_NOK_domain_bins_s3.benchmark"
        shell:
            "cat {input.passing_contigs} | cut -f2 | grep -F -v -w -f - {input.all_contigs} "
            "| cut -f1,2 | sort | uniq > {output} || true"

    """
    Remove contigs from bins.
    """
    rule clean_bins:
        input:
            gc_cov_nok="{run}/contig_mapping/bin_clean_1",
            taxa_nok="{run}/contig_mapping/bin_clean_tax3",
            bin=config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            temp("{run}/clean_bins/clean.{bin}.log")
        params:
            "{run}/clean_bins"
        benchmark:
            "{run}/benchmark/clean_bins.{bin}.benchmark"
        shell:
            #"Scripts/cleanBins.sh " +str(config["bin_dir"])+" "+str(config["bin_ext"])+ " {input.gc_cov_nok} {input.taxa_nok} {params}"
            "Scripts/cleanSingleBin.sh {input.bin} "+str(config["bin_ext"])+ "  {input.gc_cov_nok} {input.taxa_nok} {params}"

    """
    Puts the individual files into one file.
    """
    rule concat_clean_bins_log:
            input:
                expand("{run}/clean_bins/clean.{bin}.log",  bin=config["bin_list"], run=run)
            output:
                "{run}/clean_bins/clean.log"
            benchmark:
                "{run}/benchmark/concat_clean_bins_log.benchmark"
            shell:
                "cat {input} > {output}"

else:
    """
    Creates the clean bins file to make sure the pipeline has an output file.
    """
    rule skip_clean_bins:
            output:
                "{run}/clean_bins/clean.log"
            params:
                "{run}/clean_bins"
            benchmark:
                "{run}/benchmark/skip_clean_bins.benchmark"
            shell:
                "echo \"The cleaning WF was not executed.\" > {output}"


if config["generate_sql_database"] == "T":
    """
    Runs if we want to store our annotations into the FAnnP database.
    Collects the run ID, run name, start data time (end data time too but later), work directory, version, and user name.
    This information will be put in a TSV table for inserting it into the database.
    """
    rule run_to_sql:
        output:
            temp("{run}/sql/run.tsv")
        params:
            run_id=run_id,
            run=run,
            date_time=date_time
        benchmark:
            "{run}/benchmark/run_to_sql.benchmark"
        shell:
            "bash Scripts/sql_run_parser.sh {params.run_id} {params.run} '{params.date_time}' {output}"

    """
    Collects the bin ID, run ID, bin path, bin ext, taxonomy, bin length, bin avg GC, and bin avg cov.
    The taxonomy, bin avg GC and bin avg cov are extracted from the Meta-Cascabel files if available.
    Runs slightly different script with SKIP setting that gathers the same information. 
    """
    rule bins_to_sql:
        input:
            expand(config["bin_dir"]+"{bin}."+config["bin_ext"], bin=config["bin_list"])
        output:
            temp("{run}/sql/bin.tsv")
        params:
            bin_ext=config["bin_ext"],
            run_id=run_id,
            final_bins=config["meta_cascabel"]["location_final_bins"] if config["meta_cascabel"]["annotate"] == "T" else [],
            contig_coverage=config["meta_cascabel"]["contig_coverage"] if config["meta_cascabel"]["annotate"] == "T" else []
        benchmark:
            "{run}/benchmark/bins_to_sql.benchmark"
        shell:
            "bash Scripts/sql_bin_skip_parser.sh '{input}' {params.bin_ext} {params.run_id} {output}" if config["GENE_CALLING"] == "SKIP" else
            "bash Scripts/sql_bin_parser.sh '{input}' {params.bin_ext} {params.run_id} {output} {params.final_bins} {params.contig_coverage}"

    if config["meta_data"]["annotate"] == "T":
        """
        Collects the bin ID, run ID, column name, and column value.
        Format Example:
        #BinID  column1    column2  ect.
        G1.2wc-1    20  5.67    ect.
        G1.2wc-2    40  8.41    ect.
        """
        rule metadata_to_sql:
            input:
                config["meta_data"]["bin_meta_data"]
            output:
                temp("{run}/sql/metadata.tsv")
            params:
                run_id=run_id
            benchmark:
                "{run}/benchmark/metadata_to_sql.benchmark"
            shell:
                "bash Scripts/sql_metadata_parser.sh {input} {params.run_id} {output}"

    """
    Collects the contig ID, bin ID, oldname, contig length, contig GC, and contig cov.
    The contig cov is extracted from the Meta-Cascabel file if available.
    Runs slightly different script with SKIP setting that gathers the same information. 
    """
    rule contigs_to_sql:
        input:
            expand(config["bin_dir"]+"{bin}."+config["bin_ext"], bin=config["bin_list"])
        output:
            temp("{run}/sql/contig.tsv")
        params:
            bin_ext=config["bin_ext"],
            rename=config["rename_bins"],
            contig_coverage=config["meta_cascabel"]["contig_coverage"] if config["meta_cascabel"]["annotate"] == "T" else []
        benchmark:
            "{run}/benchmark/contigs_to_sql.benchmark"
        shell:
            "bash Scripts/sql_contig_skip_parser.sh '{input}' {params.bin_ext} {output}" if config["GENE_CALLING"] == "SKIP" else
            "bash Scripts/sql_contig_parser.sh '{input}' {params.bin_ext} {params.rename} {output} {params.contig_coverage}"

    """
    Collects the feature (gene) ID, contig ID, bin ID, name, start loc, end loc, dir, ftype, and definition.
    Runs slightly different script with SKIP setting that gathers the same information. 
    """
    rule gff_to_sql:
        input:
            expand("{run}/gene_calling/{bin}/{bin}.faa", bin=config["bin_list"], run=run) if config["GENE_CALLING"] == "SKIP" else
            expand("{run}/gene_calling/{bin}/{bin}.gff", bin=config["bin_list"], run=run)
        output:
            temp("{run}/sql/feature.tsv")
        params:
            bin_ext=config["bin_ext"],
            rename=config["rename_bins"]
        benchmark:
            "{run}/benchmark/gff_to_sql.benchmark"
        shell:
            "bash Scripts/sql_gene_skip_parser.sh '{input}' {params.bin_ext} {params.rename} {output}" if config["GENE_CALLING"] == "SKIP" else
            "bash Scripts/sql_gene_parser.sh '{input}' {params.rename} {output}"
    
    """
    Collects the run ID, database, db version, def file, and db file.
    Only for databases without their own table. 
    Databases: arCOG_hmmr, tigr_hmmr, cazy_hmmr, signalP, merops_blast, transporterDB_blast, 
                hydDB_blast, db1_hmmr, db2_hmmr, db1_blast, db2_blast
    """
    rule annotation_to_sql:
        output:
            temp("{run}/sql/annotation.tsv")
        params:
            run_id=run_id
        benchmark:
            "{run}/benchmark/annotation_to_sql.benchmark"
        shell:
            "python Scripts/sql_annotation_parser.py -r {params.run_id} -c "+ config_path +" -o {output}"

"""
Uses the whole NCBI non-redundant protein database and matches similar sequences to our proteins with diamond blastp.
"""
rule diamond_prots:
    input:
        bins="{run}/gene_calling/All_bins.faa", #if config["concatenate_bins"] == "T" else
        #"{run}/gene_calling/renamed/{bin}.faa"
        log_length="{run}/gene_calling/renamed_checked_log.txt",
        log_duplicate="{run}/gene_calling/duplications_checked_log.txt"
    output:
        "{run}/diamond/All_bins.tsv"  #if config["concatenate_bins"] == "T" else
        #"{run}/diamond/{bin}.tsv"
    benchmark:
        "{run}/benchmark/diamond_prots.All_bins.benchmark" #if config["concatenate_bins"] == "T" else
        #"{run}/benchmark/diamond.{bin}.benchmark"
    threads:
        int(config["diamond"]["threads"])
    shell:  # Used to have --seq instead of --max-target-seqs, also used --taxonmap in previous version.
        "diamond blastp -q {input.bins} --evalue "+ str(config["diamond"]["evalue"]) + " --threads "+ str(config["diamond"]["threads"]) +" "  
        " --max-target-seqs "+ str(config["diamond"]["seq"]) + " --db "+ str(config["diamond"]["db"]) +" "
        "--outfmt 6 qseqid qtitle qlen sseqid salltitles slen qstart qend sstart send evalue bitscore length pident staxids "
        "-o {output} " + str(config["diamond"]["extra_params"])

"""
Formats the diamond hits in a certain way by only printing specific columns.
"""
rule merge_diamond_results:
    input:
        #expand"({run}/diamond/{bin}.tsv"
        "{run}/diamond/All_bins.tsv"
    output:
        temp("{run}/diamond/tmp")
    benchmark:
        "{run}/benchmark/merge_diamond_results.benchmark"
    shell:
        "cat {input} | awk -F'\\t' -v OFS=\"\\t\" '{{ print $1, $5, $6, $11, $12, $14, $15 }}' > {output}"

"""
Adds an entry for every gene ID, meaning a gene ID without a feature gets an empty entry.
"""
rule parse_diamond_result:
    input:
        dmd="{run}/diamond/tmp",
        prot="{run}/contig_mapping/1_Bins_to_protein_list.txt"
    output:#mark as temp
        "{run}/diamond/tmp_2"
    benchmark:
        "{run}/benchmark/parse_diamond_result.benchmark"
    shell:
        "python Scripts/parse_diamond_blast_results_id_taxid.py -i {input.prot} -d {input.dmd} -o {output}"

"""
Removes spaces from the feature ID/accession and shortens it.
"""
rule format_diamond_result:
    input:
        "{run}/diamond/tmp_2"
    output:#mark as temp
        "{run}/diamond/tmp_3"
    benchmark:
        "{run}/benchmark/format_diamond_result.benchmark"
    shell:
        "cat {input} | sed 1d | awk -F\"\\t\" -v OFS=\"\\t\" '{{for(i=2;i<=NF;i+=2)gsub(/[[:blank:]]/,\"_\",$i)}}1'  | "
#the python script above sometimes leaves an empty 7th column, this gets rid of that issue
        "awk -F'\\t' -v OFS=\"\\t\"  '{{if (!$7) {{print $1,$2, $4 , $6, \"-\"}} else {{print $1, $2, $4, $6, $7}}}}' "
        " | LC_ALL=C sort | "
        #split columns with two tax ids
        "awk -F'\\t' -v OFS='\\t' '{{split($5,a,\";\"); print $1, $2, $3, $4, a[1]}}' | "
        #in column 2 remove everything after < (otherwise the name can get too long)
        "awk -F'\\t' -v OFS='\\t' '{{split($2,a,\"<\"); print $1, a[1], $3, $4, $5}} '"
        "> {output}"

"""
Add the taxonomy name from cross referencing the taxonomy ID in a database file.
"""
rule merge_diamond_taxa:
    input:
        "{run}/diamond/tmp_3"
    output:#mark as temp
        "{run}/diamond/tmp_tax"
    benchmark:
        "{run}/benchmark/merge_diamond_taxa.benchmark"
    shell:
        "LC_ALL=C join -a1 -1 5 -2 1 -e'-' -t $'\\t'  -o1.1,1.2,1.3,1.4,1.5,2.2  "
        "<(LC_ALL=C sort -k5  {input}) <(LC_ALL=C sort -k1 "+config["diamond"]["taxid_to_taxonomy"]+" ) | "
        "LC_ALL=C  sort > {output}"

"""
Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
"""
rule add_header_diamond_taxa:
    input:
        "{run}/diamond/tmp_tax"
    output:
        "{run}/diamond/diamond_map.tsv"
    benchmark:
        "{run}/benchmark/add_header_diamond_taxa.benchmark"
    shell:
        "echo -e \"accession\\tDiamond_TopHit\\tE_value\\tPecID\\tTaxID\\tTaxString\" "
        "| cat - {input} > {output}"

# Unused currently.
rule combine_diamond_to_stats:
    input:
        dmd="{run}/contig_mapping/Diamond_map.tsv",
        contigs="{run}/contig_mapping/B_GenomeInfo.txt"
    output:
        "{run}/contig_mapping/Diamond_map_tmp.tsv"
#    benchmark:
#        "{run}/benchmark/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/benchmark/diamond.{bin}.benchmark"
    shell:
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.dmd} {input.contigs} > {output}"

if config["generate_sql_database"] == "T":
    """
    Collects the feature ID, acc, desc, tax ID, tax, query start, q end, subject start, s end, eval, bitscore, and identity.
    """
    rule diamond_to_sql:
        input:
            sql="{run}/diamond/All_bins.tsv"
        output:
            "{run}/sql/diamond.tsv"
        benchmark:
            "{run}/benchmark/diamond_to_sql.benchmark"
        shell:
            "bash Scripts/sql_diamond_parser.sh {input.sql} "+ str(config["diamond"]["taxid_to_taxonomy"]) +" {output}"


if config["ko_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the KO profile database.
    """
    rule ko_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/kfam/All_bins_ko.out.tmp" if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/kfam/{bin}.out"
        benchmark:
            "{run}/benchmark/ko_hmmr.All_bins.benchmark"  if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/ko_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["ko_hmmr"]["hmmr_tool"]
        threads:
            int(config["ko_hmmr"]["cpus"])
        shell:
            #hmmr header: # target name, accession, query name, accession, E-value, score, bias, E-value, score, bias, exp, reg, clu, ov, env,dom, rep, inc, description of target
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["ko_hmmr"]["cpus"]) +" --notextw "+ str(config["ko_hmmr"]["extra_params"])+" "
            " -E " +str(config["ko_hmmr"]["evalue"]) + " " + str(config["ko_hmmr"]["database"]) + " {input}  > {output}"
    
    """
    Puts the individual files into one file.
    """
    rule concat_ko_hmmr:
        input:
            expand("{run}/kfam/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/kfam/All_bins_ko.out.tmp"
        benchmark:
            "{run}/benchmark/concat_ko_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["ko_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_ko_to_hmmsearch:
            input:
                "{run}/kfam/All_bins_ko.out.tmp"
            output:
                "{run}/kfam/All_bins_ko_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_ko_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_ko_hmmr_domains:
        input:
            "{run}/kfam/All_bins_ko.out.tmp" if config["ko_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/kfam/All_bins_ko_parsed.out.tmp"
        output:
            "{run}/kfam/All_bins_ko_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_ko_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_ko_hmmr_search:
        input:
            "{run}/kfam/All_bins_ko.out.tmp" if config["ko_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/kfam/All_bins_ko_parsed.out.tmp"
        output:
            "{run}/kfam/All_bins_ko_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_ko_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_ko_hmmr_count:
        input:
            sort="{run}/kfam/All_bins_ko_sorted.out.tmp",
            counts="{run}/kfam/All_bins_ko_count.out.tmp"
        output:
            "{run}/kfam/All_bins_ko_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_ko_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_ko_hmmr_all:
        input:
            "{run}/kfam/All_bins_ko_merged.out.tmp"
        output:
            normal="{run}/kfam/All_bins_ko.out",
            sql="{run}/kfam/All_bins_ko_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_ko_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["ko_hmmr"]["max_domain_overlap"]) + " -c " + str(config["ko_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["ko_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"



#   rule format_ko_hmmr:
#       input:
#           "{run}/kfam/All_bins_pfam.out"
#       output:
#           temp("{run}/kfam/kfam_sorted_cols")
#       shell:
#           "cat {input} | cut -f1,3,5,6  > {output}"

    
    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_ko_hmmr:
        input:
            ko= "{run}/kfam/All_bins_ko.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/kfam/kfam_merged"
        benchmark:
            "{run}/benchmark/merge_ko_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.ko}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick 
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
    
    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_ko_hmmr:
        input:
            ko="{run}/kfam/kfam_merged",
        output:#mark as temp
            "{run}/kfam/kfam_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_ko_hmmr.benchmark"
        shell:
            "bash Scripts/awk_ko_parser.sh {input.ko} "+config["ko_hmmr"]["database_names"]+" {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_ko_hmmr:
        input:
            "{run}/kfam/kfam_names_tmp"
        output:
            "{run}/kfam/kfam_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_ko_hmmr.benchmark"
        shell:
            "echo -e \"accession\\tKO_hmm\\te_value\\tbit_score\\tbit_score_cutoff\\tDefinition\\tconfidence\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the feature ID, KO ID, start, end, evalue, and score.
        """
        rule ko_to_sql:
            input:
                sql="{run}/kfam/All_bins_ko_sql.out"
            output:
                "{run}/sql/ko.tsv"
            benchmark:
                "{run}/benchmark/ko_to_sql.benchmark"
            shell:
                "bash Scripts/sql_ko_parser.sh {input.sql} {output}"


if config["arCOG_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the arCOG profile database.
    """
    rule arCOG_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/arCOG/All_bins_arCOG.out.tmp" if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/arCOG/{bin}.out"
        benchmark:
            "{run}/benchmark/arCOG_hmmr.All_bins.benchmark"  if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/arCOG_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["arCOG_hmmr"]["hmmr_tool"]
        threads:
            int(config["arCOG_hmmr"]["cpus"])
        shell:
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["arCOG_hmmr"]["cpus"]) +" --notextw "+ str(config["arCOG_hmmr"]["extra_params"])+" "
            " -E " +str(config["arCOG_hmmr"]["evalue"]) + " " + str(config["arCOG_hmmr"]["database"]) + " {input}  "
            " > {output}" 

    """
    Puts the individual files into one file.
    """
    rule concat_arCOG_hmmr:
        input:
            expand("{run}/arCOG/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/arCOG/All_bins_arCOG.out.tmp"
        benchmark:
            "{run}/benchmark/concat_arCOG_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["arCOG_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_hmmscan_to_hmmsearch:
            input:
                "{run}/arCOG/All_bins_arCOG.out.tmp"
            output:
                "{run}/arCOG/All_bins_arCOG_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_hmmscan_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_arCOG_hmmr_domains:
        input:
            "{run}/arCOG/All_bins_arCOG.out.tmp" if config["arCOG_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/arCOG/All_bins_arCOG_parsed.out.tmp"
        output:
            "{run}/arCOG/All_bins_arCOG_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_arCOG_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_arCOG_hmmr_search:
        input:
            "{run}/arCOG/All_bins_arCOG.out.tmp" if config["arCOG_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/arCOG/All_bins_arCOG_parsed.out.tmp"
        output:
            "{run}/arCOG/All_bins_arCOG_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_arCOG_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_arCOG_hmmr_count:
        input:
            sort="{run}/arCOG/All_bins_arCOG_sorted.out.tmp",
            counts="{run}/arCOG/All_bins_arCOG_count.out.tmp"
        output:
            "{run}/arCOG/All_bins_arCOG_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_arCOG_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_arCOG_hmmr_all:
        input:
            "{run}/arCOG/All_bins_arCOG_merged.out.tmp"
        output:
            normal="{run}/arCOG/All_bins_arCOG.out",
            sql="{run}/arCOG/All_bins_arCOG_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_arCOG_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["arCOG_hmmr"]["max_domain_overlap"]) + " -c " + str(config["arCOG_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["arCOG_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    """
    Formats the arCOG ID differently (for FunctionalAnnotation.tsv).
    """
    rule format_arCOG_hmmr:
        input:
            "{run}/arCOG/All_bins_arCOG.out"        
        output:#mark as temp
            "{run}/arCOG/arcog_tmp"
        benchmark:
            "{run}/benchmark/format_arCOG_hmmr.benchmark"
        shell:
            "cat {input} | awk  -v OFS='\\t' '{{split($2,a,\".\"); print $1, a[1], $3,$4}}' | "
            "LC_ALL=C sort > {output}"
    
    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_arCOG_hmmr:
        input:
            arc="{run}/arCOG/arcog_tmp",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/arCOG/arcog_merged"
        benchmark:
            "{run}/benchmark/merge_arCOG_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.arc}) | LC_ALL=C sort  > {output}"
    
    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_arCOG_hmmr:
        input:
            arc="{run}/arCOG/arcog_merged"
        output:#mark as temp
            temp("{run}/arCOG/arcogs_names_tmp")
        benchmark:
            "{run}/benchmark/add_names_arCOG_hmmr.benchmark"
        shell:
            "bash Scripts/awk_arCOG_parser.sh {input.arc} "+ config["arCOG_hmmr"]["database_names"] +" {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_arCOG_hmmr:
        input:
            "{run}/arCOG/arcogs_names_tmp"
        output:
            "{run}/arCOG/arcogs_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_arCOG_hmmr.benchmark"
        shell:
            "echo -e \"accession\\tarcogs\\tarcogs_geneID\\tarcogs_Description\\tPathway\\tarcogs_evalue\" | "
            "cat - {input} > {output}"

    # Unused currently.
    rule merge_arCOG_hmmr_results:
#    """
#    deprecated, now we map everything together from "{run}/arCOG/arcogs_map.tsv"
#    """
        input:
            arc="{run}/arCOG/arcogs_map.tsv",
            map="{run}/contig_mapping/Diamond_map_tmp.tsv"
        output:#mark as temp
            "{run}/contig_mapping/Arcogs_map_tmp.tsv"
    #    benchmark:
    #        "{run}/benchmark/arCOG.All_bins.benchmark"  if config["concatenate_bins"] == "T" else
    #        "{run}/benchmark/arCOG.{bin}.benchmark"
        shell:
            "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.arc} {input.map} > {output}"  

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity (N/A).
        """
        rule arCOG_to_sql:
            input:
                sql="{run}/arCOG/All_bins_arCOG_sql.out"
            output:
                "{run}/sql/arCOG.tsv"
            params:
                db_id="01"
            benchmark:
                "{run}/benchmark/arCOG_to_sql.benchmark"
            shell:
                "bash Scripts/sql_arCOG_parser.sh {input.sql} "+ str(config["arCOG_hmmr"]["database_names"]) +" {params.db_id} {output}"

if config["cog_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the COG profile database.
    """
    rule cog_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["cog_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/cog/All_bins_cog.out.tmp" if config["cog_hmmr"]["concatenate_bins"] == "T" else
            "{run}/cog/{bin}.out"
        benchmark:
            "{run}/benchmark/cog_hmmr.All_bins.benchmark"  if config["cog_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/cog_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["cog_hmmr"]["hmmr_tool"]
        benchmark:
            "{run}/benchmark/concat_cog_hmmr.benchmark"
        threads:
            int(config["cog_hmmr"]["cpus"])
        shell:
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["cog_hmmr"]["cpus"]) +" --notextw "+ str(config["cog_hmmr"]["extra_params"])+" "
            " -E " +str(config["cog_hmmr"]["evalue"]) + " " + str(config["cog_hmmr"]["database"]) + " {input}  "
            " > {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_cog_hmmr:
        input:
            expand("{run}/cog/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/cog/All_bins_cog.out.tmp"
        benchmark:
            "{run}/benchmark/concat_cog_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["cog_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_cog_to_hmmsearch:
            input:
                "{run}/cog/All_bins_cog.out.tmp"
            output:
                "{run}/cog/All_bins_cog_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_cog_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_cog_hmmr_domains:
        input:
            "{run}/cog/All_bins_cog.out.tmp" if config["cog_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/cog/All_bins_cog_parsed.out.tmp"
        output:
            "{run}/cog/All_bins_cog_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_cog_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_cog_hmmr_search:
        input:
            "{run}/cog/All_bins_cog.out.tmp" if config["cog_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/cog/All_bins_cog_parsed.out.tmp"
        output:
            "{run}/cog/All_bins_cog_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_cog_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_cog_hmmr_count:
        input:
            sort="{run}/cog/All_bins_cog_sorted.out.tmp",
            counts="{run}/cog/All_bins_cog_count.out.tmp"
        output:
            "{run}/cog/All_bins_cog_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_cog_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_cog_hmmr_all:
        input:
            "{run}/cog/All_bins_cog_merged.out.tmp"
        output:
            normal="{run}/cog/All_bins_cog.out",
            sql="{run}/cog/All_bins_cog_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_cog_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["cog_hmmr"]["max_domain_overlap"]) + " -c " + str(config["cog_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["cog_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    """
    Formats the COG ID differently (for FunctionalAnnotation.tsv).
    """
    rule format_cog_hmmr:
        input:
            "{run}/cog/All_bins_cog.out"
        output:#mark as temp
            "{run}/cog/cog_tmp"
        benchmark:
            "{run}/benchmark/format_cog_hmmr.benchmark"
        shell:
            "cat {input} | awk  -v OFS='\\t' '{{split($2,a,\".\"); print $1, a[1], $3,$4}}' | "
            "LC_ALL=C sort > {output}"

    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_cog_hmmr:
        input:
            cog="{run}/cog/cog_tmp",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/cog/cog_merged"
        benchmark:
            "{run}/benchmark/merge_cog_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.cog}) | LC_ALL=C sort  > {output}"

#    rule rm_cog_def_spaces:
#        output:
#            temp("{run}/cog/cog_names.tsv")
#        shell:
#            "cat " + str(config["cog_hmmr"]["database_names"])+" | sed 's/ /_/g' > {output}"

    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_cog_hmmr:
        input:
            cog="{run}/cog/cog_merged"
#            db="{run}/cog/cog_names.tsv"
        output:#mark as temp
            temp("{run}/cog/cogs_names_tmp")
        benchmark:
            "{run}/benchmark/add_names_cog_hmmr.benchmark"
        shell:
            "bash Scripts/awk_cog_parser.sh {input.cog} "+config["cog_hmmr"]["database_names"]+" {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_cog_hmmr:
        input:
            "{run}/cog/cogs_names_tmp"
        output:
            "{run}/cog/cogs_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_cog_hmmr.benchmark"
        shell:
    #add in headers
            "echo -e \"accession\\tNCBI_COG\\tNCBI_COG_Description\\tCOG_PathwayID\\tCOG_Pathway\\tNCBI_COG_evalue\" | "
            "cat - {input} > {output}"
    
    # Unused currently.
    rule merge_cog_hmmr_results:
#    """
#    deprecated, now we map everything together from "{run}/cog/cogs_map.tsv"
#    """
        input:
            arc="{run}/cog/cogs_map.tsv",
            map="{run}/contig_mapping/Diamond_map_tmp.tsv"
        output:#mark as temp
            "{run}/contig_mapping/cogs_map_tmp.tsv"
    #    benchmark:
    #        "{run}/benchmark/cog.All_bins.benchmark"  if config["concatenate_bins"] == "T" else
    #        "{run}/benchmark/cog.{bin}.benchmark"
        shell:
            "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.arc} {input.map} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the feature ID, COG ID, start, end, evalue, and score.
        """
        rule cog_to_sql:
            input:
                sql="{run}/cog/All_bins_cog_sql.out"
            output:
                "{run}/sql/cog.tsv"
            benchmark:
                "{run}/benchmark/cog_to_sql.benchmark"
            shell:
                "bash Scripts/sql_cog_parser.sh {input.sql} {output}"

if config["pfam_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the Pfam profile database.
    """
    rule pfam_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/pfam/All_bins_pfam.out.tmp" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/pfam/{bin}.out"
        benchmark:
            "{run}/benchmark/pfam_hmmr.All_bins.benchmark" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/pfam_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["pfam_hmmr"]["hmmr_tool"]
        threads:
            int(config["pfam_hmmr"]["cpus"])
        shell:
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["pfam_hmmr"]["cpus"]) +" --notextw "+ str(config["pfam_hmmr"]["extra_params"])+" "
            " -E " +str(config["pfam_hmmr"]["evalue"]) + " "  + str(config["pfam_hmmr"]["database"]) + " {input}  > {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_pfam_hmmr:
        input:
            expand("{run}/pfam/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/pfam/All_bins_pfam.out.tmp"
        benchmark:
            "{run}/benchmark/concat_pfam_hmmr.benchmark"
        shell:
            "cat {input} > {output}"
        
#     rule format_pfam_hmmr:
#        input:
#            "{run}/pfam/All_bins_pfam.out"
#        output:
#            temp("{run}/pfam/pfam_sorted_cols")
#        shell:
#            "cat {input} | cut  -f1,3,5,6 | sed 's/ /\t/g' > {output}"

    if config["pfam_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_pfam_to_hmmsearch:
            input:
                "{run}/pfam/All_bins_pfam.out.tmp"
            output:
                "{run}/pfam/All_bins_pfam_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_pfam_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_pfam_hmmr_domains:
        input:
            "{run}/pfam/All_bins_pfam.out.tmp" if config["pfam_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/pfam/All_bins_pfam_parsed.out.tmp"
        output:
            "{run}/pfam/All_bins_pfam_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_pfam_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_pfam_hmmr_search:
        input:
            "{run}/pfam/All_bins_pfam.out.tmp" if config["pfam_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/pfam/All_bins_pfam_parsed.out.tmp"
        output:
            "{run}/pfam/All_bins_pfam_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_pfam_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_pfam_hmmr_count:
        input:
            sort="{run}/pfam/All_bins_pfam_sorted.out.tmp",
            counts="{run}/pfam/All_bins_pfam_count.out.tmp"
        output:
            "{run}/pfam/All_bins_pfam_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_pfam_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_pfam_hmmr_all:
        input:
            "{run}/pfam/All_bins_pfam_merged.out.tmp"
        output:
            normal="{run}/pfam/All_bins_pfam.out",
            sql="{run}/pfam/All_bins_pfam_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_pfam_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["pfam_hmmr"]["max_domain_overlap"]) + " -c " + str(config["pfam_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["pfam_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_pfam_hmmr:
        input:
            pfam="{run}/pfam/All_bins_pfam.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/pfam/pfam_merged"
        benchmark:
            "{run}/benchmark/merge_pfam_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.pfam}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick 
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
    
    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_pfam_hmmr:
        input:
            pfam="{run}/pfam/pfam_merged",
        output:#mark as temp
            "{run}/pfam/pfam_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_pfam_hmmr.benchmark"
        shell:
            "bash Scripts/awk_pfam_parser.sh {input.pfam} "+config["pfam_hmmr"]["database_names"]+" {output}"

#    rule identify_high_confidence_pfam_hmmr:
#       input:
#           "{run}/pfam/pfam_names_tmp"
#       output:#mark as temp
#           temp("{run}/pfam/pfam_names_tmp_2")
#       shell:
#           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_pfam_hmmr:
        input:
            "{run}/pfam/pfam_names_tmp"
        output:
            "{run}/pfam/pfam_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_pfam_hmmr.benchmark"
        shell:
            "echo -e \"accession\\tPFAM_ID\\tPFAM_hmm\\tPFAM_description\\tPfam_Evalue\\tPfam_Score\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the feature ID, Pfam ID, start, end, evalue, and score.
        """
        rule pfam_to_sql:
            input:
                sql="{run}/pfam/All_bins_pfam_sql.out"
            output:
                "{run}/sql/pfam.tsv"
            benchmark:
                "{run}/benchmark/pfam_to_sql.benchmark"
            shell:
                "bash Scripts/sql_pfam_parser.sh {input.sql} "+ str(config["pfam_hmmr"]["database_names"]) +" {output}"

if config["signalP"]["annotate"] == "T":
    rule signalP:
        input:
            "{run}/gene_calling/All_bins.faa" if config["signalP"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/signalP/prediction_results.txt" if config["signalP"]["concatenate_bins"] == "T" else
            "{run}/signalP/{bin}/prediction_results.txt"
        benchmark:
            "{run}/benchmark/signalP.All_bins.benchmark" if config["signalP"]["concatenate_bins"] == "T" else
            "{run}/benchmark/signalP.{bin}.benchmark"
        params:
            "{run}/signalP/signalP.log" if config["signalP"]["concatenate_bins"] == "T" else
            "{run}/signalP/signalP.{bin}.log",
            "{run}/signalP/" if config["signalP"]["concatenate_bins"] == "T" else
            "{run}/signalP/{bin}/"
        threads:
            int(config["signalP"]["cpus"])
        shell:
            "signalp6 --write_procs "+ str(config["signalP"]["cpus"]) +" "+ str(config["signalP"]["extra_params"])+" "
            " -org " +str(config["signalP"]["organism"]) + " --fastafile {input} --output_dir {params[1]} 2> {params[0]}"

    """
    Selects the feature hit with the best score.
    """
    rule select_best_signalP:
        '''
        Output header looks like this
        # SignalP-6.0   Organism: Other Timestamp: 20220912095354
        # ID    Prediction      OTHER   SP(Sec/SPI)     LIPO(Sec/SPII)  TAT(Tat/SPI)    TATLIPO(Sec/SPII)       PILIN(Sec/SPIII)        CS Position
        G4sed-100.C1_1  OTHER   1.000033        0.000000        0.000000        0.000000        0.000000        0.000000
        G4sed-100.C1_2  SP      0.000259        0.999112        0.000186        0.000150        0.000148        0.000135        CS pos: 28-29. Pr: 0.9589
        After format we have this:
        G4sed-100.C1_2  SP      0.999112        CS pos: 28-29. Pr: 0.9589
        G4sed-100.C1_3  LIPO    0.999852        CS pos: 21-22. Pr: 0.9865
        '''
        input:
            "{run}/signalP/prediction_results.txt" if config["signalP"]["concatenate_bins"] == "T" else
            expand("{run}/signalP/{bin}/prediction_results.txt",bin=config["bin_list"], run=run)
        output:
            "{run}/signalP/prediction_summary.txt"
        benchmark:
            "{run}/benchmark/select_best_signalP.benchmark"
        shell:
            "cat {input} | awk -F\"\\t\"  'BEGIN{{OFS=\"\\t\"}}; "
            "$0 !~ \"^#\" && $2 !~ \"OTHER\"{{max=0;for(i=3;i<NF;i++){{if($i>max){{max=$i}}}};split($1,id,\" \"); print id[1],$2,max,$NF}}' > {output}"

    """
    Adds gene entries for every gene. Even if they didn't have a feature hit (for FunctionalAnnotation.tsv).
    """
    rule merge_signalP:
        input:
            signalP="{run}/signalP/prediction_summary.txt",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/signalP/signalP_merged"
        benchmark:
            "{run}/benchmark/merge_signalP.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.signalP}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_signalP:
        input:
            "{run}/signalP/signalP_merged"
        output:
            "{run}/signalP/signalP_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_signalP.benchmark"
        shell:
            "echo -e \"accession\\tsignalP\\tLikelihood\\tCleavage_site_position\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description (N/A), evalue, score (N/A) and identity.
        """
        rule signalP_to_sql:
            input:
                sql="{run}/signalP/prediction_summary.txt"
            output:
                "{run}/sql/signalP.tsv"
            params:
                db_id="04"
            benchmark:
                "{run}/benchmark/signalP_to_sql.benchmark"
            shell:
                "bash Scripts/sql_signalP_parser.sh {input.sql} {params.db_id} {output}"

            
if config["tigr_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the TIGR profile database.
    """
    rule tigr_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/tigr/All_bins_tigr.out.tmp" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/tigr/{bin}.out"
        benchmark:
            "{run}/benchmark/tigr_hmmr.All_bins.benchmark" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/tigr_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["tigr_hmmr"]["hmmr_tool"]
        threads:
            int(config["tigr_hmmr"]["cpus"])
        shell:
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["tigr_hmmr"]["cpus"]) +" --notextw "+ str(config["tigr_hmmr"]["extra_params"])+" "
            " -E " +str(config["tigr_hmmr"]["evalue"]) + " "  + str(config["tigr_hmmr"]["database"]) + " {input} > {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_tigr_hmmr:
        input:
            expand("{run}/tigr/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/tigr/All_bins_tigr.out.tmp"
        benchmark:
            "{run}/benchmark/concat_tigr_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["tigr_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_tigr_to_hmmsearch:
            input:
                "{run}/tigr/All_bins_tigr.out.tmp"
            output:
                "{run}/tigr/All_bins_tigr_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_tigr_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_tigr_hmmr_domains:
        input:
            "{run}/tigr/All_bins_tigr.out.tmp" if config["tigr_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/tigr/All_bins_tigr_parsed.out.tmp"
        output:
            "{run}/tigr/All_bins_tigr_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_tigr_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_tigr_hmmr_search:
        input:
            "{run}/tigr/All_bins_tigr.out.tmp" if config["tigr_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/tigr/All_bins_tigr_parsed.out.tmp"
        output:
            "{run}/tigr/All_bins_tigr_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_tigr_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_tigr_hmmr_count:
        input:
            sort="{run}/tigr/All_bins_tigr_sorted.out.tmp",
            counts="{run}/tigr/All_bins_tigr_count.out.tmp"
        output:
            "{run}/tigr/All_bins_tigr_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_tigr_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_tigr_hmmr_all:
        input:
            "{run}/tigr/All_bins_tigr_merged.out.tmp"
        output:
            normal="{run}/tigr/All_bins_tigr.out",
            sql="{run}/tigr/All_bins_tigr_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_tigr_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["tigr_hmmr"]["max_domain_overlap"]) + " -c " + str(config["tigr_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["tigr_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_tigr_hmmr:
        input:
            tigr="{run}/tigr/All_bins_tigr.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/tigr/tigr_merged"
        benchmark:
            "{run}/benchmark/merge_tigr_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.tigr}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_tigr_hmmr:
        input:
            tigr="{run}/tigr/tigr_merged",
        output:#mark as temp
            "{run}/tigr/tigr_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_tigr_hmmr.benchmark"
        shell:
            "bash Scripts/awk_tigr_parser.sh {input.tigr} "+config["tigr_hmmr"]["database_names"]+" {output}"
#    rule identify_high_confidence_tigr_hmmr:
#       input:
#           "{run}/tigr/tigr_names_tmp"
#       output:#mark as temp
#           temp("{run}/tigr/tigr_names_tmp_2")
#       shell:
#           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_tigr_hmmr:
        input:
            "{run}/tigr/tigr_names_tmp"
        output:
            "{run}/tigr/tigr_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_tigr_hmmr.benchmark"
        shell:
            "echo -e \"accession\\tTIGR_hmm\\tTIGR_name\\tTIGR_description\\tTIGR_EC\\tTIGR_gene\\tTIGR_Evalue\\tTIGR_Score\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity (N/A).
        """
        rule tigr_to_sql:
            input:
                sql="{run}/tigr/All_bins_tigr_sql.out"
            output:
                "{run}/sql/tigr.tsv"
            params:
                db_id="02"
            benchmark:
                "{run}/benchmark/tigr_to_sql.benchmark"
            shell:
                "bash Scripts/sql_tigr_parser.sh {input.sql} "+ str(config["tigr_hmmr"]["database_names"]) +" {params.db_id} {output}"


if config["cazy_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using the CAZy profile database.
    """
    rule cazy_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/cazy/All_bins_cazy.out.tmp" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/cazy/{bin}.out"
        benchmark:
            "{run}/benchmark/cazy_hmmr.All_bins.benchmark" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/cazy_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["cazy_hmmr"]["hmmr_tool"]
        threads:
            int(config["cazy_hmmr"]["cpus"])
        shell:
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["cazy_hmmr"]["cpus"]) +" --notextw "+ str(config["cazy_hmmr"]["extra_params"])+" "
            " -E " +str(config["cazy_hmmr"]["evalue"]) + " "  + str(config["cazy_hmmr"]["database"]) + " {input}  > {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_cazy_hmmr:
        input:
            expand("{run}/cazy/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/cazy/All_bins_cazy.out.tmp"
        benchmark:
            "{run}/benchmark/concat_cazy_hmmr.benchmark"
        shell:
            "cat {input} | sed 's/\.hmm//'  > {output}"
    
    if config["cazy_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_cazy_to_hmmsearch:
            input:
                "{run}/cazy/All_bins_cazy.out.tmp"
            output:
                "{run}/cazy/All_bins_cazy_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_cazy_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_cazy_hmmr_domains:
        input:
            "{run}/cazy/All_bins_cazy.out.tmp" if config["cazy_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/cazy/All_bins_cazy_parsed.out.tmp"
        output:
            "{run}/cazy/All_bins_cazy_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_cazy_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_cazy_hmmr_search:
        input:
            "{run}/cazy/All_bins_cazy.out.tmp" if config["cazy_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/cazy/All_bins_cazy_parsed.out.tmp"
        output:
            "{run}/cazy/All_bins_cazy_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_cazy_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_cazy_hmmr_count:
        input:
            sort="{run}/cazy/All_bins_cazy_sorted.out.tmp",
            counts="{run}/cazy/All_bins_cazy_count.out.tmp"
        output:
            "{run}/cazy/All_bins_cazy_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_cazy_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_cazy_hmmr_all:
        input:
            "{run}/cazy/All_bins_cazy_merged.out.tmp"
        output:
            normal="{run}/cazy/All_bins_cazy.out",
            sql="{run}/cazy/All_bins_cazy_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_cazy_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["cazy_hmmr"]["max_domain_overlap"]) + " -c " + str(config["cazy_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["cazy_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"
    
    """
    Adds gene entries for every gene. Even if they didn't have a domain hit (for FunctionalAnnotation.tsv).
    """
    rule merge_cazy_hmmr:
        input:
            cazy="{run}/cazy/All_bins_cazy.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/cazy/cazy_merged"
        benchmark:
            "{run}/benchmark/merge_cazy_hmmr.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.cazy}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    """
    Add information to the feature/domain hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_cazy_hmmr:
        input:
            cazy="{run}/cazy/cazy_merged",
        output:#mark as temp
            "{run}/cazy/cazy_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_cazy_hmmr.benchmark"
        shell:
            # "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,2.2,1.3,1.4 "
            # "<(LC_ALL=C sort -k2 {input.cazy}) "
            # "<(LC_ALL=C sort -k1 "+str(config["cazy_hmmr"]["database_names"])+") | LC_ALL=C  sort > {output}"
            "bash Scripts/awk_cazy_parser.sh {input.cazy} "+config["cazy_hmmr"]["database_names"]+" {output}"
#    rule identify_high_confidence_cazy_hmmr:
#       input:
#           "{run}/cazy/cazy_names_tmp"
#       output:#mark as temp
#           temp("{run}/cazy/cazy_names_tmp_2")
#       shell:
#           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_cazy_hmmr:
        input:
            "{run}/cazy/cazy_names_tmp"
        output:
            "{run}/cazy/cazy_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_cazy_hmmr.benchmark"
        shell:
            "echo -e \"accession\\tCAZY_hmm\\tCAZY_description\\tCAZY_Evalue\\tCAZY_Score\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity (N/A).
        """
        rule cazy_to_sql:
            input:
                sql="{run}/cazy/All_bins_cazy_sql.out"
            output:
                "{run}/sql/cazy.tsv"
            params:
                db_id="03"
            benchmark:
                "{run}/benchmark/cazy_to_sql.benchmark"
            shell:
                "bash Scripts/sql_cazy_parser.sh {input.sql} "+ str(config["cazy_hmmr"]["database_names"]) +" {params.db_id} {output}"


if config["merops_blast"]["annotate"] == "T":
    """
    Uses the MEROPS database and matches similar sequences to our proteins with blastp.
    """
    rule merops_blast:
        input:
            "{run}/gene_calling/All_bins.faa" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/merops/All_bins_merops.out.tmp" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.out"
        benchmark:
            "{run}/benchmark/merops_blast.All_bins.benchmark" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/merops_blast.{bin}.benchmark"
        threads:
            int(config["merops_blast"]["threads"])
        shell:
            "blastp -num_threads  "+ str(config["merops_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["merops_blast"]["database"]) + " "
            " -out {output} -evalue " +str(config["merops_blast"]["evalue"]) + " "+str(config["merops_blast"]["extra_params"])

    """
    Selects the feature hit with the best score.
    """
    rule select_best_merops_blast:
        input:
            "{run}/merops/All_bins_merops.out.tmp" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.out"
        output:
            "{run}/merops/All_bins_merops.out" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.unq.out"
        benchmark:
            "{run}/benchmark/select_best_merops_blast.All_bins.benchmark" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/select_best_merops_blast.{bin}.benchmark"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_merops_blast:
        input:
            expand("{run}/merops/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/merops/All_bins_merops.out"
        benchmark:
            "{run}/benchmark/concat_merops_blast.benchmark"
        shell:
            "cat {input}  > {output}"
    
    """
    Adds gene entries for every gene. Even if they didn't have a feature hit (for FunctionalAnnotation.tsv).
    """
    rule merge_merops_blast:
        input:
            merops="{run}/merops/All_bins_merops.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/merops/merops_merged"
        benchmark:
            "{run}/benchmark/merge_merops_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.11 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.merops}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    """
    Add information to the feature hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_merops_blast:
        input:
            merops="{run}/merops/merops_merged"
        output:#mark as temp
            "{run}/merops/merops_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_merops_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,2.2,1.3,1.4 "
            "<(LC_ALL=C sort -k2 {input.merops}) "
            "<(LC_ALL=C sort -k1 "+str(config["merops_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_merops_blast:
        input:
            "{run}/merops/merops_names_tmp"
        output:
            "{run}/merops/merops_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_merops_blast.benchmark"
        shell:
            "echo -e \"accession\\tMEROPS_blast\\tMEROPS_description\\tMEROPS_Idenity\\tMEROPS_Evalue\" | "
            "cat - {input} > {output}"    

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity.
        """
        rule merops_to_sql:
            input:
                sql="{run}/merops/All_bins_merops.out"
            output:
                "{run}/sql/merops.tsv"
            params:
                db_id="05"
            benchmark:
                "{run}/benchmark/merops_to_sql.benchmark"
            shell:
                "bash Scripts/sql_merops_parser.sh {input.sql} "+ str(config["merops_blast"]["database_names"]) +" {params.db_id} {output}"


if config["transporterDB_blast"]["annotate"] == "T":
    """
    Uses the Transporter Classification database and matches similar sequences to our proteins with blastp.
    """
    rule transporterDB_blast:
        input:
            "{run}/gene_calling/All_bins.faa" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/transporter/All_bins_transporter.out.tmp" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/{bin}.out"
        benchmark:
            "{run}/benchmark/transporterDB_blast.All_bins.benchmark" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/transporterDB_blast.{bin}.benchmark"
        threads:
            int(config["transporterDB_blast"]["threads"])
        shell:
            "blastp -num_threads  "+ str(config["transporterDB_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["transporterDB_blast"]["database"]) + " "
            " -out {output} -evalue " +str(config["transporterDB_blast"]["evalue"]) + " "+str(config["transporterDB_blast"]["extra_params"])

    """
    Selects the feature hit with the best score.
    """
    rule select_best_transporter_blast:
        input:
            "{run}/transporter/All_bins_transporter.out.tmp" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/{bin}.out"
        output:
            "{run}/transporter/All_bins_transporter.out.tmp2" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/{bin}.unq.out"
        benchmark:
            "{run}/benchmark/select_best_transporter_blast.All_bins.benchmark" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/select_best_transporter_blast.{bin}.benchmark"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_transporter_blast:
        input:
            expand("{run}/transporter/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/transporter/All_bins_transporter.out"
        benchmark:
            "{run}/benchmark/concat_transporter_blast.benchmark"
        shell:
            "cat {input} "
            "| LC_ALL=C sort > {output}"
    
    """
    Sorts the merged file.
    """
    rule concat_transporter_blast_all:
        input:
            "{run}/transporter/All_bins_transporter.out.tmp2"
        output:
            "{run}/transporter/All_bins_transporter.out"
        benchmark:
            "{run}/benchmark/concat_transporter_blast_all.benchmark"
        shell:
            "cat {input} "
            "| LC_ALL=C sort > {output}"
    
    """
    Adds gene entries for every gene. Even if they didn't have a feature hit (for FunctionalAnnotation.tsv).
    """
    rule merge_transporter_blast:
        input:
            transporter="{run}/transporter/All_bins_transporter.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/transporter/transporter_merged"
        benchmark:
            "{run}/benchmark/merge_transporter_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.transporter} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $11}}') | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    """
    Add information to the feature hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_transporter_blast:
        input:
            transporter="{run}/transporter/transporter_merged",
        output:#mark as temp
            "{run}/transporter/transporter_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_transporter_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o 1.1,1.2,2.2,2.3,1.3 "
            "<(LC_ALL=C sort -k2 {input.transporter}) "
            "<(LC_ALL=C sort -k1 "+str(config["transporterDB_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_transporter_blast:
        input:
            "{run}/transporter/transporter_names_tmp"
        output:
            "{run}/transporter/transporter_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_transporter_blast.benchmark"
        shell:
            "echo -e \"accession\\tTransporter_id\\tTransporter_id2\\tTransporter_description\\tTransporter_Evalue\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity.
        """
        rule transporter_to_sql:
            input:
                sql="{run}/transporter/All_bins_transporter.out"
            output:
                "{run}/sql/transporter.tsv"
            params:
                db_id="06"
            benchmark:
                "{run}/benchmark/transporter_to_sql.benchmark"
            shell:
                "bash Scripts/sql_transporter_parser.sh {input.sql} "+ str(config["transporterDB_blast"]["database_names"]) +" {params.db_id} {output}"


if config["hydDB_blast"]["annotate"] == "T":
    """
    Uses the Hydrogenase database and matches similar sequences to our proteins with blastp.
    """
    rule hydDB_blast:
        input:
            "{run}/gene_calling/All_bins.faa" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/hydDB/All_bins_hydDB.out.tmp" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/hydDB/{bin}.out"
        benchmark:
            "{run}/benchmark/hydDB_blast.All_bins.benchmark" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/hydDB_blast.{bin}.benchmark"
        threads:
            int(config["hydDB_blast"]["threads"])
        shell:
            "blastp -num_threads  "+ str(config["hydDB_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["hydDB_blast"]["database"]) + " "
            " -out {output} -evalue " +str(config["hydDB_blast"]["evalue"]) + " "+str(config["hydDB_blast"]["extra_params"])

    """
    Selects the feature hit with the best score.
    """
    rule select_best_hydDB_blast:
        input:
            "{run}/hydDB/All_bins_hydDB.out.tmp" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/hydDB/{bin}.out"
        output:
            "{run}/hydDB/All_bins_hydDB.out" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/hydDB/{bin}.unq.out"
        benchmark:
            "{run}/benchmark/select_best_hydDB_blast.All_bins.benchmark" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/select_best_hydDB_blast.{bin}.benchmark"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_hydDB_blast:
        input:
            expand("{run}/hydDB/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/hydDB/All_bins_hydDB.out"
        benchmark:
            "{run}/benchmark/concat_hydDB_blast.benchmark"
        shell:
            "cat {input}  > {output}"
    
    """
    Adds gene entries for every gene. Even if they didn't have a feature hit (for FunctionalAnnotation.tsv).
    """
    rule merge_hydDB_blast:
        input:
            hydDB="{run}/hydDB/All_bins_hydDB.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/hydDB/hydDB_merged"
        benchmark:
            "{run}/benchmark/merge_hydDB_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.11,2.12 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.hydDB}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
    
    """
    Add information to the feature hits for this database (for FunctionalAnnotation.tsv).
    """
    rule add_names_hydDB_blast:
        input:
            hydDB="{run}/hydDB/hydDB_merged",
        output:#mark as temp
            "{run}/hydDB/hydDB_names_tmp"
        benchmark:
            "{run}/benchmark/add_names_hydDB_blast.benchmark"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,0,2.2,1.4 "
            "<(LC_ALL=C sort -k2 {input.hydDB}) "
            "<(LC_ALL=C sort -k1 "+str(config["hydDB_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

    """
    Add headers to the columns that will be added to the FunctionalAnnotation.tsv table.
    """
    rule add_header_hydDB_blast:
        input:
            "{run}/hydDB/hydDB_names_tmp"
        output:
            "{run}/hydDB/hydDB_map.tsv"
        benchmark:
            "{run}/benchmark/add_header_hydDB_blast.benchmark"
        shell:
            "echo -e \"accession\\tHydDB\\tDescription\\tHydDB_evalue\" | "
            "cat - {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity.
        """
        rule hydDB_to_sql:
            input:
                sql="{run}/hydDB/All_bins_hydDB.out"
            output:
                "{run}/sql/hydDB.tsv"
            params:
                db_id="07"
            benchmark:
                "{run}/benchmark/hydDB_to_sql.benchmark"
            shell:
                "bash Scripts/sql_hydDB_parser.sh {input.sql} "+ str(config["hydDB_blast"]["database_names"]) +" {params.db_id} {output}"


#
# Custom database rules
#
if config["db1_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using a custom profile database.
    """
    rule db1_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["db1_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/db1_hmmr/All_bins_db1_hmmr.out.tmp" if config["db1_hmmr"]["concatenate_bins"] == "T" else
            "{run}/db1_hmmr/{bin}.out"
        benchmark:
            "{run}/benchmark/db1_hmmr.All_bins.benchmark"  if config["db1_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/db1_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["db1_hmmr"]["hmmr_tool"]
        threads:
            int(config["db1_hmmr"]["cpus"])
        shell:
            #hmmr header: # target name, accession, query name, accession, E-value, score, bias, E-value, score, bias, exp, reg, clu, ov, env,dom, rep, inc, description of target
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["db1_hmmr"]["cpus"]) +" --notextw "+ str(config["db1_hmmr"]["extra_params"])+" "
            " -E " +str(config["db1_hmmr"]["evalue"]) + " " + str(config["db1_hmmr"]["database"]) + " {input}  > {output}"
    
    """
    Puts the individual files into one file.
    """
    rule concat_db1_hmmr:
        input:
            expand("{run}/db1_hmmr/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/db1_hmmr/All_bins_db1_hmmr.out.tmp"
        benchmark:
            "{run}/benchmark/concat_db1_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["db1_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_db1_hmmr_to_hmmsearch:
            input:
                "{run}/db1_hmmr/All_bins_db1_hmmr.out.tmp"
            output:
                "{run}/db1_hmmr/All_bins_db1_hmmr_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_db1_hmmr_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_db1_hmmr_domains:
        input:
            "{run}/db1_hmmr/All_bins_db1_hmmr.out.tmp" if config["db1_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/db1_hmmr/All_bins_db1_hmmr_parsed.out.tmp"
        output:
            "{run}/db1_hmmr/All_bins_db1_hmmr_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_db1_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_db1_hmmr_search:
        input:
            "{run}/db1_hmmr/All_bins_db1_hmmr.out.tmp" if config["db1_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/db1_hmmr/All_bins_db1_hmmr_parsed.out.tmp"
        output:
            "{run}/db1_hmmr/All_bins_db1_hmmr_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_db1_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_db1_hmmr_count:
        input:
            sort="{run}/db1_hmmr/All_bins_db1_hmmr_sorted.out.tmp",
            counts="{run}/db1_hmmr/All_bins_db1_hmmr_count.out.tmp"
        output:
            "{run}/db1_hmmr/All_bins_db1_hmmr_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_db1_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_db1_hmmr_all:
        input:
            "{run}/db1_hmmr/All_bins_db1_hmmr_merged.out.tmp"
        output:
            normal="{run}/db1_hmmr/All_bins_db1_hmmr.out",
            sql="{run}/db1_hmmr/All_bins_db1_hmmr_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_db1_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["db1_hmmr"]["max_domain_overlap"]) + " -c " + str(config["db1_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["db1_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity (N/A).
        """
        rule db1_hmmr_to_sql:
            input:
                sql="{run}/db1_hmmr/All_bins_db1_hmmr_sql.out"
            output:
                "{run}/sql/db1_hmmr.tsv"
            params:
                db_id="08"
            benchmark:
                "{run}/benchmark/db1_hmmr_to_sql.benchmark"
            shell:
                "bash Scripts/sql_hmmr_feature_parser.sh {input.sql} {params.db_id} {output}"


if config["db2_hmmr"]["annotate"] == "T":
    """
    Looks for homologs using a HMMER tool (hmmsearch or hmmscan) using a custom profile database.
    """
    rule db2_hmmr:
        input:
            "{run}/gene_calling/All_bins.faa" if config["db2_hmmr"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/db2_hmmr/All_bins_db2_hmmr.out.tmp" if config["db2_hmmr"]["concatenate_bins"] == "T" else
            "{run}/db2_hmmr/{bin}.out"
        benchmark:
            "{run}/benchmark/db2_hmmr.All_bins.benchmark"  if config["db2_hmmr"]["concatenate_bins"] == "T" else
            "{run}/benchmark/db2_hmmr.{bin}.benchmark"
        params:
            hmm_tool=config["db2_hmmr"]["hmmr_tool"]
        threads:
            int(config["db2_hmmr"]["cpus"])
        shell:
            #hmmr header: # target name, accession, query name, accession, E-value, score, bias, E-value, score, bias, exp, reg, clu, ov, env,dom, rep, inc, description of target
            "{params.hmm_tool} --domtblout  /dev/stdout -o  /dev/null --cpu "+ str(config["db2_hmmr"]["cpus"]) +" --notextw "+ str(config["db2_hmmr"]["extra_params"])+" "
            " -E " +str(config["db2_hmmr"]["evalue"]) + " " + str(config["db2_hmmr"]["database"]) + " {input}  > {output}"
    
    """
    Puts the individual files into one file.
    """
    rule concat_db2_hmmr:
        input:
            expand("{run}/db2_hmmr/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/db2_hmmr/All_bins_db2_hmmr.out.tmp"
        benchmark:
            "{run}/benchmark/concat_db2_hmmr.benchmark"
        shell:
            "cat {input} > {output}"

    if config["db2_hmmr"]["hmmr_tool"] == "hmmscan":
        """
        If hmmscan is used, changes the hmmscan format to the hmmsearch format.
        """
        rule parse_db2_hmmr_to_hmmsearch:
            input:
                "{run}/db2_hmmr/All_bins_db2_hmmr.out.tmp"
            output:
                "{run}/db2_hmmr/All_bins_db2_hmmr_parsed.out.tmp"
            benchmark:
                "{run}/benchmark/parse_db2_hmmr_to_hmmsearch.benchmark"
            shell:
                """cat {input} | egrep -v "^#" | awk '{{print $4,$5,$3,$1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output}"""

    """
    Counts how many homolog hits each gene name has.
    """
    rule count_db2_hmmr_domains:
        input:
            "{run}/db2_hmmr/All_bins_db2_hmmr.out.tmp" if config["db2_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/db2_hmmr/All_bins_db2_hmmr_parsed.out.tmp"
        output:
            "{run}/db2_hmmr/All_bins_db2_hmmr_count.out.tmp"
        benchmark:
            "{run}/benchmark/count_db2_hmmr_domains.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 | awk '{{print $1}}' | uniq -c | awk '{{print $2,$1}}' | sed 's/ /\t/g' > {output}"""

    """
    Sorts the homolog hits based on gene name and secondly its overall bit score.
    """
    rule sort_db2_hmmr_search:
        input:
            "{run}/db2_hmmr/All_bins_db2_hmmr.out.tmp" if config["db2_hmmr"]["hmmr_tool"] == "hmmsearch" else
            "{run}/db2_hmmr/All_bins_db2_hmmr_parsed.out.tmp"
        output:
            "{run}/db2_hmmr/All_bins_db2_hmmr_sorted.out.tmp"
        benchmark:
            "{run}/benchmark/sort_db2_hmmr_search.benchmark"
        shell:
            """cat {input} | egrep -v "^#" | sort -k1,1 -k8gr > {output}"""
    
    """
    Paste the counts in front of the gene name so we can use it for the homolog domain parser.
    """
    rule merge_db2_hmmr_count:
        input:
            sort="{run}/db2_hmmr/All_bins_db2_hmmr_sorted.out.tmp",
            counts="{run}/db2_hmmr/All_bins_db2_hmmr_count.out.tmp"
        output:
            "{run}/db2_hmmr/All_bins_db2_hmmr_merged.out.tmp"
        benchmark:
            "{run}/benchmark/merge_db2_hmmr_count.benchmark"
        shell:
            "join {input.counts} {input.sort} | sed 's/ /\t/g' > {output}"

    """
    Stores multiple HMMER domains based on score and feature overlap.
    Can create an output in two different formats, the first one stores all information per gene name (for FunctionalAnnotation.tsv).
    The second one is used for the creation of the TSV table for storing it into the MySQL FAnnP database.
    The format it uses is storing the information per domain hit.
    """
    rule parse_best_db2_hmmr_all:
        input:
            "{run}/db2_hmmr/All_bins_db2_hmmr_merged.out.tmp"
        output:
            normal="{run}/db2_hmmr/All_bins_db2_hmmr.out",
            sql="{run}/db2_hmmr/All_bins_db2_hmmr_sql.out" if config["generate_sql_database"] == "T" else []
        params:
            sql_tag="-m --sql" if config["generate_sql_database"] == "T" else []
        benchmark:
            "{run}/benchmark/parse_best_db2_hmmr_all.benchmark"
        shell:
            "python Scripts/parse_hmmr_domains.py -i {input} -p " + str(config["db2_hmmr"]["max_domain_overlap"]) + " -c " + str(config["db2_hmmr"]["evalue_domain_cutoff"]) + " "
            "-d " + str(config["db2_hmmr"]["domain_col"]) + " -o {output.normal} {params.sql_tag} {output.sql}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity (N/A).
        """
        rule db2_hmmr_to_sql:
            input:
                sql="{run}/db2_hmmr/All_bins_db2_hmmr_sql.out"
            output:
                "{run}/sql/db2_hmmr.tsv"
            params:
                db_id="09"
            benchmark:
                "{run}/benchmark/db2_hmmr_to_sql.benchmark"
            shell:
                "bash Scripts/sql_hmmr_feature_parser.sh {input.sql} {params.db_id} {output}"


if config["db1_blast"]["annotate"] == "T":
    """
    Uses a custom database and matches similar sequences to our proteins with blastp.
    """
    rule db1_blast:
        input:
            "{run}/gene_calling/All_bins.faa" if config["db1_blast"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/db1_blast/All_bins_db1_blast.out.tmp" if config["db1_blast"]["concatenate_bins"] == "T" else
            "{run}/db1_blast/{bin}.out"
        benchmark:
            "{run}/benchmark/db1_blast.All_bins.benchmark" if config["db1_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/db1_blast.{bin}.benchmark"
        threads:
            int(config["db1_blast"]["threads"])
        shell:
            "blastp -num_threads  "+ str(config["db1_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["db1_blast"]["database"]) + " "
            " -out {output} -evalue " +str(config["db1_blast"]["evalue"]) + " "+str(config["db1_blast"]["extra_params"])

    """
    Selects the feature hit with the best score.
    """
    rule select_best_db1_blast:
        input:
            "{run}/db1_blast/All_bins_db1_blast.out.tmp" if config["db1_blast"]["concatenate_bins"] == "T" else
            "{run}/db1_blast/{bin}.out"
        output:
            "{run}/db1_blast/All_bins_db1_blast.out" if config["db1_blast"]["concatenate_bins"] == "T" else
            "{run}/db1_blast/{bin}.unq.out"
        benchmark:
            "{run}/benchmark/select_best_db1_blast.All_bins.benchmark" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/select_best_db1_blast.{bin}.benchmark"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_db1_blast:
        input:
            expand("{run}/db1_blast/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/db1_blast/All_bins_db1_blast.out"
        benchmark:
            "{run}/benchmark/concat_db1_blast.benchmark"
        shell:
            "cat {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity.
        """
        rule db1_blast_to_sql:
            input:
                sql="{run}/db1_blast/All_bins_db1_blast.out"
            output:
                "{run}/sql/db1_blast.tsv"
            params:
                db_id="10"
            benchmark:
                "{run}/benchmark/db1_blast_to_sql.benchmark"
            shell:
                "bash Scripts/sql_blast_feature_parser.sh {input.sql} {params.db_id} {output}"


if config["db2_blast"]["annotate"] == "T":
    """
    Uses a custom database and matches similar sequences to our proteins with blastp.
    """
    rule db2_blast:
        input:
            "{run}/gene_calling/All_bins.faa" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/gene_calling/renamed/{bin}.faa"
        output:
            "{run}/db2_blast/All_bins_db2_blast.out.tmp" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/db2_blast/{bin}.out"
        benchmark:
            "{run}/benchmark/db2_blast.All_bins.benchmark" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/db2_blast.{bin}.benchmark"
        threads:
            int(config["db2_blast"]["threads"])
        shell:
            "blastp -num_threads  "+ str(config["db2_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["db2_blast"]["database"]) + " "
            " -out {output} -evalue " +str(config["db2_blast"]["evalue"]) + " "+str(config["db2_blast"]["extra_params"])

    """
    Selects the feature hit with the best score.
    """
    rule select_best_db2_blast:
        input:
            "{run}/db2_blast/All_bins_db2_blast.out.tmp" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/db2_blast/{bin}.out"
        output:
            "{run}/db2_blast/All_bins_db2_blast.out" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/db2_blast/{bin}.unq.out"
        benchmark:
            "{run}/benchmark/select_best_db2_blast.All_bins.benchmark" if config["db2_blast"]["concatenate_bins"] == "T" else
            "{run}/benchmark/select_best_db2_blast.{bin}.benchmark"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

    """
    Puts the individual files into one file.
    """
    rule concat_db2_blast:
        input:
            expand("{run}/db2_blast/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/db2_blast/All_bins_db2_blast.out"
        benchmark:
            "{run}/benchmark/concat_db2_blast.benchmark"
        shell:
            "cat {input} > {output}"

    if config["generate_sql_database"] == "T":
        """
        Collects the annotation ID, feature ID, start, end, accession, description, evalue, score and identity.
        """
        rule db2_blast_to_sql:
            input:
                sql="{run}/db2_blast/All_bins_db2_blast.out"
            output:
                "{run}/sql/db2_blast.tsv"
            params:
                db_id="11"
            benchmark:
                "{run}/benchmark/db2_blast_to_sql.benchmark"
            shell:
                "bash Scripts/sql_blast_feature_parser.sh {input.sql} {params.db_id} {output}"

"""
For every gene, information and annotations will be merged into the B_GenomeInfo.txt file using annotations found by the other tools.
The annotations found by the other tools are stored in their respective TSV files when their database is enabled in the config file.
The scripts will merge the tables by using their gene names.
This will create the FunctionalAnnotation.tsv file.
"""
rule summarize_annotation:
    input:
        "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/diamond/diamond_map.tsv",
        "{run}/arCOG/arcogs_map.tsv" if config["arCOG_hmmr"]["annotate"] == "T" else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/kfam/kfam_map.tsv" if config["ko_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/pfam/pfam_map.tsv" if config["pfam_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/tigr/tigr_map.tsv" if config["tigr_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/cazy/cazy_map.tsv" if config["cazy_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/cog/cogs_map.tsv" if config["cog_hmmr"]["annotate"] == "T" else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/merops/merops_map.tsv" if config["merops_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/transporter/transporter_map.tsv" if config["transporterDB_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/hydDB/hydDB_map.tsv" if config["hydDB_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/signalP/signalP_map.tsv" if config["signalP"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt"
    output:
        "{run}/FunctionalAnnotation.tsv"
    benchmark:
        "{run}/benchmark/summarize_annotation.benchmark"
    params:
        run="{run}"
    shell:
        "Scripts/summary_annotation.sh {params.run} {input[0]}  {output}"

"""
Stores all TSV files created by the to_sql rules, into the MySQL FAnnP database.
The script uses the config file to locate the TSV files and sends them into a '/tmp/sql/fannp/{run ID}/' map.
From this map the files can be send into the database with the MySQL 'LOAD DATA INFILE' method.
The run TSV file gets stored by using a different method because MySQL needs a custom query to add datetime.
The output file is necessary for snakemake and only exists for the workflow. 
"""
rule all_sql_inserts:
    input:
        "{run}/sql/run.tsv",
        "{run}/sql/bin.tsv",
        "{run}/sql/feature.tsv",
        "{run}/sql/contig.tsv",
        "{run}/sql/annotation.tsv",
        "{run}/sql/diamond.tsv",
        "{run}/sql/metadata.tsv" if config["meta_data"]["annotate"] == "T" else [],
        "{run}/sql/ko.tsv" if config["ko_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/cog.tsv" if config["cog_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/pfam.tsv" if config["pfam_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/arCOG.tsv" if config["arCOG_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/tigr.tsv" if config["tigr_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/cazy.tsv" if config["cazy_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/signalP.tsv" if config["signalP"]["annotate"] == "T" else [],
        "{run}/sql/merops.tsv" if config["merops_blast"]["annotate"] == "T" else [],
        "{run}/sql/transporter.tsv" if config["transporterDB_blast"]["annotate"] == "T" else [],
        "{run}/sql/hydDB.tsv" if config["hydDB_blast"]["annotate"] == "T" else [],
        "{run}/sql/db1_hmmr.tsv" if config["db1_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/db2_hmmr.tsv" if config["db2_hmmr"]["annotate"] == "T" else [],
        "{run}/sql/db1_blast.tsv" if config["db1_blast"]["annotate"] == "T" else [],
        "{run}/sql/db2_blast.tsv" if config["db2_blast"]["annotate"] == "T" else []
    output:
        "{run}/sql/run_to_sql.txt"
    params:
        run_id=run_id,
        a1="--store_metadata" if config["meta_data"]["annotate"] == "T" else [],
        a2="--store_ko" if config["ko_hmmr"]["annotate"] == "T" else [],
        a3="--store_cog" if config["cog_hmmr"]["annotate"] == "T" else [],
        a4="--store_arCOG" if config["arCOG_hmmr"]["annotate"] == "T" else [],
        a5="--store_pfam" if config["pfam_hmmr"]["annotate"] == "T" else [],
        a6="--store_tigr" if config["tigr_hmmr"]["annotate"] == "T" else [],
        a7="--store_cazy" if config["cazy_hmmr"]["annotate"] == "T" else [],
        a8="--store_signalP" if config["signalP"]["annotate"] == "T" else [],
        a9="--store_merops" if config["merops_blast"]["annotate"] == "T" else [],
        a10="--store_transporter" if config["transporterDB_blast"]["annotate"] == "T" else [],
        a11="--store_hydDB" if config["hydDB_blast"]["annotate"] == "T" else [],
        a12="--store_hmmr1" if config["db1_hmmr"]["annotate"] == "T" else [],
        a13="--store_hmmr2" if config["db2_hmmr"]["annotate"] == "T" else [],
        a14="--store_blast1" if config["db1_blast"]["annotate"] == "T" else [],
        a15="--store_blast2" if config["db2_blast"]["annotate"] == "T" else []
    benchmark:
        "{run}/benchmark/all_sql_inserts.benchmark"
    shell:
        "python Scripts/generate_sql_tables.py -c "+ config_path +" -r {params.run_id} -o {output} {params.a1} {params.a2} {params.a3} "
        "{params.a4} {params.a5} {params.a6} {params.a7} {params.a8} {params.a9} {params.a10} {params.a11} {params.a12} {params.a13} {params.a14} {params.a15}"

"""
Copies the Marimo from the scripts using the run ID.
"""
rule generate_marimo_notebook:
    input:
        "{run}/sql/run_to_sql.txt"
    output:
        "{run}/sql/run_to_notebook.txt"
    params:
        run_id=run_id,
        notebook="Scripts/MySQL_FAnnP_Analyzer.py",
        mysql_user=config["mysql_username"],
        mysql_pass=config["mysql_password"],
        mysql_db=config["mysql_database"]
    benchmark:
        "{run}/benchmark/generate_marimo_notebook.benchmark"
    shell:
        "bash Scripts/create_marimo_notebook.sh {input} {params.run_id} {params.notebook} {params.mysql_user} {params.mysql_pass} {params.mysql_db} {output}"

"""
Output for all the final product we want to produce.
"""
rule report:
    input:
        "{run}/FunctionalAnnotation.tsv",
        "{run}/clean_bins/clean.log",
        "{run}/sql/run_to_notebook.txt" if config["generate_sql_database"] == "T" else []
        #"{run}/diamond/All_bins.tsv", # if config["concatenate_bins"] == "T" else
        #"{run}/diamond/{bin}.tsv",
        #"{run}/kfam/All_bins_ko.out" if config["concatenate_bins"] == "T" else
        #"{run}/kfam/{bin}.out",
        #"{run}/pfam/All_bins_pfam.out" if config["concatenate_bins"] == "T" else
        #"{run}/pfam/{bin}.out",
        #"{run}/tigr/All_bins_tigr.out" if config["concatenate_bins"] == "T" else
        #"{run}/tigr/{bin}.out",
        #"{run}/cazy/All_bins_cazy.out" if config["concatenate_bins"] == "T" else
        #"{run}/cazy/{bin}.out",
        #"{run}/merops/All_bins_merops.out" if config["concatenate_bins"] == "T" else
        #"{run}/merops/{bin}.out",
        #"{run}/transporter/All_bins_transporter.out" if config["concatenate_bins"] == "T" else
        #"{run}/transporter/{bin}.out",
        #"{run}/arCOG/All_bins_arCOG.out" if config["concatenate_bins"] == "T" else
        #"{run}/arCOG/{bin}.out",
        #"{run}/hydDB/All_bins_transporter.out" if config["concatenate_bins"] == "T" else
        #"{run}/hydDB/{bin}.out",
        #"{run}/gene_calling/{bin}.gbk.gz",
        #"{run}/contig_mapping/Contig_Old_mapping.txt",#this we can delete later
        #"{run}/contig_mapping/Contig_Old_mapping_for_merging.txt",
        #"{run}/contig_mapping/B_GenomeInfo.txt",
        #"{run}/contig_mapping/Diamond_map_tmp.tsv",
        #"{run}/contig_mapping/Arcogs_map_tmp.tsv",
        #"{run}/diamond/diamond_map.tsv",
        #"{run}/arCOG/arcogs_map.tsv",
        #"{run}/kfam/kfam_map.tsv",
        #"{run}/pfam/pfam_map.tsv",
        #"{run}/tigr/tigr_map.tsv",
        #"{run}/cazy/cazy_map.tsv",
        #"{run}/merops/merops_map.tsv",
        #"{run}/transporter/transporter_map.tsv",
        #"{run}/hydDB/hydDB_map.tsv"

    output:
        "{run}/report.txt" #if config["concatenate_bins"] == "T" else
        #"{run}/{bin}/report.txt"
    benchmark:
        "{run}/benchmark/report.benchmark"
    shell:
        "touch {output}"

