run=config["RUN"]

rule all:
    input:
        #"{run}/{bin}/report.txt"
        #dynamic("{run}/{bin}/report.txt")        
        #expand("{run}/{bin}/report.txt" if config["concatenate_bins"] == "T" else "{run}/report.txt",bin=config["bin_list"], run=run)
        expand("{run}/report.txt",bin=config["bin_list"], run=run) #if config["concatenate_bins"] == "T" else
        #expand("{run}/{bin}/report.txt",bin=config["bin_list"], run=run) 
        #"{run}/report.txt"

if config["GENE_CALLING"] == "PRODIGAL":
    rule prodigal_bins:
        input:
           config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            gbk=temp("{run}/prokka/{bin}/{bin}.gbk"),
            genes="{run}/prokka/{bin}/{bin}.ffn",
            prots=temp("{run}/prokka/{bin}/{bin}.faa")
        benchmark:
            "{run}/prokka/prodigal.benchmark"
        shell:
            "prodigal -a {output.prots} -d {output.genes} -f gbk  -i {input} "
            "-p "+ config["prodigal"]["procedure"]+"  "+ config["prodigal"]["extra_params"]+"  -o {output.gbk} "
        
else: 
    rule prokka_bins:
        input:
            config["bin_dir"]+"{bin}."+config["bin_ext"]
        output:
            "{run}/prokka/{bin}/{bin}.faa",
            "{run}/prokka/{bin}/{bin}.gbk"
        params:
            output_dir="{run}/prokka/{bin}",
            file_ext=config["bin_ext"],
            bin_prefix="{bin}"
        benchmark:
            "{run}/prokka/prokka.benchmark"
        shell:
            "prokka --prefix {wildcards.bin} --outdir {params.output_dir}  --addgenes --force --increment  "+str(config["prokka"]["increment"]) + " "
            " --compliant --centre UU --cpus "+ str(config["prokka"]["cpus"]) + " --norrna --notrna "+ config["prokka"]["extra_params"]+" {input}"

rule rename_bins:
    input:
        "{run}/prokka/{bin}/{bin}.faa"
    output:
        "{run}/prokka/renamed/{bin}.faa"
    shell:
        "awk -v run=\"{wildcards.run}\" -v bin=\"{wildcards.bin}\" '/>/{{sub(\">\",\"&\"FILENAME\"|\");sub(/\.faa/,x);sub(run\"/prokka/\"bin\"/\",x)}}1' {input} | "
        "cut -f1 -d \" \" > {output}"

rule make_protein_list:
    input:
        expand("{run}/prokka/renamed/{bin}.faa", bin=config["bin_list"], run=run)
    output:
        "{run}/contig_mapping/AllProteins_list.txt"
    shell:
        "cat {input} | grep '^>' | cut -f1 -d \" \" | sed 's/>//g' > {output}"

rule bins_to_protList:
    input:
        "{run}/contig_mapping/AllProteins_list.txt"
    output:
        "{run}/contig_mapping/1_Bins_to_protein_list.txt" 
    shell:
        "cat {input} | awk -F'\\t' -v OFS='\\t' '{{split($1,a,\"|\"); print $1, a[1]}}' | LC_ALL=C sort > {output}"

rule gzip_gbk:
    input:
        "{run}/prokka/{bin}/{bin}.gbk"
    output:
        "{run}/prokka/{bin}/{bin}.gbk.gz"
    shell:
        "gzip {input}"

if config["GENE_CALLING"] == "PRODIGAL":
    rule seqid_to_prot:
        input:
            "{run}/prokka/{bin}/{bin}.gbk.gz"
        output:
            "{run}/prokka/{bin}/{bin}.seq2prots"
        shell:
            "zcat {input} | grep -A1 --no-group-separator \"CDS\" | "
            "grep \"/note\" | sed 's/.*\\/note=\"//' | "
            "cut -f1 -d\";\" | sed 's/ID=// ; s/_/\\t/' > {output}"
    rule parse_gbk_files_prodigal:
        input:
            gbk="{run}/prokka/{bin}/{bin}.gbk.gz",
            map="{run}/prokka/{bin}/{bin}.seq2prots"
        output:
            temp("{run}/prokka/{bin}/Contig_Protein_list.txt")
        shell:
            "zcat {input.gbk} | grep DEFINITION | cut -f3 -d\" \" | "
            "cut -f1-3 -d\";\" | sed 's/;/\\t/g ; s/seq.[a-z]*=//g ; s/\"//g' | "
            "awk -F\"\\t\" 'NR==FNR{{contig[$1]=$3;len[$1]=$2;next}} BEGIN{{OFS=\"\\t\"}} {{print $1,contig[$1],len[$1],contig[$1]\"_\"$2,\"NA\"}}' "
            "- {input.map} > {output}"
    rule concat_gbk_files_prodigal:
        input:
            expand("{run}/prokka/{bin}/Contig_Protein_list.txt",bin=config["bin_list"], run=run)
        output:
            "{run}/contig_mapping/Contig_Protein_list.txt"
        shell:
            "cat {input} > {output}"
else:
    rule create_list_gbk:
        input:
            expand("{run}/prokka/{bin}/{bin}.gbk.gz", bin=config["bin_list"], run=run)
        output:
            "{run}/prokka/list_of_gzips"
        shell:
            "echo {input} | sed 's/ /\\n/g' > {run}/prokka/list_of_gzips"
#it seems that the previous command generates on extra line which 
#ends up in a error.
    rule parse_gbk_files:
        input:
            "{run}/prokka/list_of_gzips"
        output:
            "{run}/contig_mapping/Contig_Protein_list.txt"
        shell:
            "python Scripts/Parse_prokka_for_MAGs_from_gbk-file.py -i {input} -t {output}"

rule protein_to_binID:
    input:
        "{run}/contig_mapping/1_Bins_to_protein_list.txt"
    output:
        "{run}/contig_mapping/Prokka_to_BinID.txt"
    shell:
        "cat {input} | awk 'BEGIN{{OFS=\"\\t\"}}{{split($1,a,\"|\"); print a[2],$2}}' | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{parts=split($1,a,\"_\"); contig_id=a[1]; for(i=2;i<parts;i++){{contig_id=contig_id\"_\"a[i]}}; print contig_id,$2}}' | sort | uniq > {output}"

rule link_old_to_new_binIDS:
    input:
        prokka="{run}/contig_mapping/Prokka_to_BinID.txt",
        prot="{run}/contig_mapping/Contig_Protein_list.txt"
    output:
        "{run}/contig_mapping/temp_Bin_Contig_Protein_list.txt"
    params:
        idx="2" if config["GENE_CALLING"] == "PRODIGAL" else "1"
    shell:
        "awk -v idx={params.idx} 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$idx]}}' {input.prokka} {input.prot} | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,$7,$2,$2,$3,$4,$5}}'  > {output}" 

#add in extra column for contig nr
rule add_contig_nr:
    input:
        "{run}/contig_mapping/temp_Bin_Contig_Protein_list.txt"
    output:
        "{run}/contig_mapping/Bin_Contig_Protein_list.txt"
    shell:
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{ print $2\"_contig_\"$1,$1,$2,$3,$1,$5,$6,$7}}' "
        "{input} > {output}"  if config["GENE_CALLING"] == "PRODIGAL" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{split($4,a,\"_\"); print $2\"_contig_\"a[2],$1,$2,$3,a[2],$5,$6,$7}}' "
        "{input} > {output}" 

rule merge_contid_ids:
    input:
        contig_prot="{run}/contig_mapping/Bin_Contig_Protein_list.txt",
        contig_old="{run}/contig_mapping/Contig_Old_mapping_for_merging.txt"
    output:
        "{run}/contig_mapping/Bin_Contig_Protein_list_merged.txt"
    shell:#headers: accession, BinID, newContigID, oldContigID, mergeContigID,ContigLengthNew, LengthContigOld, GC,ProteinID, prokka
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.contig_old} {input.contig_prot} | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $3\"|\"$7,$3,$4,$10,$1,$6,$14,$13,$7,$8}}' > {output}" 


rule computeGC_bins:
    input:
        config["bin_dir"]+"{bin}."+config["bin_ext"]
    output:
        temp("{run}/contig_mapping/{bin}_temp1")
    shell:
        "perl Scripts/length+GC.pl {input} > {output}"

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

rule concatenate_contig_data:
    input:
         expand("{run}/contig_mapping/{bin}_temp2", bin=config["bin_list"], run=run)
    output:
        "{run}/contig_mapping/Contig_Old_mapping.txt"
    shell:
        "cat {input} > {output}" 
 
rule create_merging_file:
    input:
        "{run}/contig_mapping/Contig_Old_mapping.txt"
    output:
        "{run}/contig_mapping/Contig_Old_mapping_for_merging.txt"
    shell:
        "cat {input} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $2\"_contig_\"$3,$0}}' > {output}"


rule concatenate_bins:
    input:
        expand("{run}/prokka/renamed/{bin}.faa",  bin=config["bin_list"], run=run)
    output:
        temp("{run}/prokka/All_bins.faa")
    benchmark:
        "{run}/prokka/concatenate_bins.benchmark"
    shell:
        "cat {input} > {output}"

rule clean_concatenated_bins:
    input:
        "{run}/prokka/All_bins.faa"
    output:
        "{run}/prokka/All_bins_clean.faa"
    shell:
        "cut -f1 -d \" \" {input} > {output}"

rule length_GC_prots:
    input:
        "{run}/prokka/All_bins_clean.faa"
    output:#mark as temp
        temp("{run}/contig_mapping/gc_base")
    shell:
        "perl Scripts/length+GC.pl {input} > {output}"

#prepare a file with protein length and gc to add protein length into our file
#perl ~/../spang_team/Scripts/Others/length+GC.pl ../Prokka/All_Genomes_clean.faa > temp

rule merge_contig_protein_info:
    input:
        base="{run}/contig_mapping/gc_base",
        merged="{run}/contig_mapping/Bin_Contig_Protein_list_merged.txt"
    output:#mark as temp
        "{run}/contig_mapping/gc_base_tmp"
    shell:
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' "
        "{input.base} {input.merged} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$12,$13,$10}}' "
        "> {output}"

#add in taxon string
#Comment: Bin_to_tax.txt info was copied from bin_stats
#1st column = binID, 2nd column = tax string from gtdb (can be exchange by other tax info of course)
rule merge_with_taxonomy:
    input:
        "{run}/contig_mapping/gc_base_tmp"
    output:#mark as tmp
        "{run}/contig_mapping/gc_base_taxo"
    shell:
        "cat " + config["bin_taxonomy"]["taxonomy_file"]+ " | "
        "awk -v id=\""+ config["bin_taxonomy"]["bin_col"]+ "\" -v tax=\""+ config["bin_taxonomy"]["taxonomy_col"]+ "\"  'NR>1{{print $id\"\\t\"$tax}}' | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}'"
        " - {input} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $2,$1,$14,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' > {output}"
        if config["bin_taxonomy"]["add_taxonomy"] == "T" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $2,$1,\"NA\",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' {input} > {output}"
 
#bin GC and covergae       
rule merge_with_coverage:
    input:
        "{run}/contig_mapping/gc_base_taxo"
    output:#mark as tmp
        "{run}/contig_mapping/gc_base_taxo_cov"
    shell:
        "cat " + config["bin_cleaning"]["bin_coverage_gc_file"]+ " | "
        "awk -v id=\""+ config["bin_cleaning"]["bin_col"]+ "\" -v cov=\""+ config["bin_cleaning"]["avg_depth_col"]+ "\""
        " -v gc=\""+ config["bin_cleaning"]["avg_gc_col"]+ "\" ' NR>1{{print $id\"\\t\"$cov\"\\t\"$gc}}' | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$2\"\\t\"$3;next}}BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}"
        "{{print $1,$2,$3,a[$2],$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}'"
        " - {input}  > {output}"
        if config["bin_cleaning"]["clean_bins"] == "T" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print  $1,$2,$3,\"NA\",\"NA\",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' {input} > {output}"

#contig covergae
rule merge_with_contig_cov:
    input:
        "{run}/contig_mapping/gc_base_taxo_cov"
    output:#mark as tmp
        "{run}/contig_mapping/gc_base_taxo_cov_2"
    shell:
        "cat " + config["bin_cleaning"]["contig_coverage_file"]+ " | "
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$2;next}}BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}"
        "{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,a[$7],$12,$13,$14,$15}}'"
        " - {input}  > {output}"
        if config["bin_cleaning"]["clean_bins"] == "T" else
        "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,\"NA\",$12,$13,$14,$15}}' {input} > {output}"


rule add_header_to_file:
    input:
        "{run}/contig_mapping/gc_base_taxo_cov_2"
    output:
        "{run}/contig_mapping/B_GenomeInfo.txt"
    shell:
        "echo -e \"accession\\tBinID\\tTaxString_Bin\\tBin_avg_cov\\tBin_avg_gc\\tNewContigID\\tOldContigId\\tContigIdMerge"
        "\\tContigNewLength\\tContigOldLength\\tContig_GC\\tContig_cov\\tProteinID\\tProteinGC\\tProteinLength\\tProkka\" | "
        "cat - {input} > {output}"

if config["bin_cleaning"]["clean_bins"] == "T":
    rule mark_NOK_gc_coverage_bins:
            input:
                "{run}/contig_mapping/gc_base_taxo_cov_2"
            output:
                "{run}/contig_mapping/bin_clean_1"
            shell:
                "cat {input} |  awk -F \"\\t\" 'function abs(v) {{return v < 0 ? -v : v}} "
                "BEGIN{{OFS=\"\\t\"}} {{pass_gc=\"PASS\";pass_cvg=\"PASS\"; gc_diff=abs($11-$5); "  
                "if(gc_diff>"+str(config["bin_cleaning"]["rules"]["GC_max_diff"])+"){{pass_gc=\"FAIL\"}}; "
                "ab_ratio=$12/$4;if(ab_ratio<1-(" +str(config["bin_cleaning"]["rules"]["Coverage_ratio"])+ ") || "
                "ab_ratio>1+("+str(config["bin_cleaning"]["rules"]["Coverage_ratio"])+")){{pass_cvg=\"FAIL\"}}; "
                "if(pass_cvg==\"FAIL\" || pass_gc==\"FAIL\"){{print $2,$6,$7,gc_diff,pass_gc,ab_ratio,pass_cvg}}}}' |"
                "sort | uniq > {output}"


    rule mark_NOK_domain_bins_s1:
            """ 
            This rule creates a file with the following columns
            Bin_ID, ContigID, BinTaxonomy, ProtTaxonomy, NumProts with that taxonomy, Percentage, NumProts in the contig.
            """  
            input:
                "{run}/FunctionalAnnotation.tsv"
            output:
                "{run}/contig_mapping/bin_clean_tax1"
            shell:
                "cat {input} | sed 's/d__//g' |  awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\"}} "
                "{{if(NR==1){{for(i=1;i<=NF;i++){{if($i==\"TaxString\"){{blastTax=i; }} }}" #Identify column with diamond/blast taxonomy, this can ve variable
                "}}else{{split($3,bin_tax,\";\");split($blastTax,prot_tax,\",\");" #split taxonomies 
                "if(!bin[$2]){{bin[$2]=$2}};" #select uniq bins
                "if(!contig[$7]){{contig[$7]=bin_tax[1];bin_contig[$2][$7]=$7}};"#select uniq contigs and link them to its taxonomy and to its bin
                "bin_prot[$7][prot_tax[1]]++;contig_num[$7]++}};}}END{{" #count prot domains per contigs
                "for(b in bin){{" #flush the arrays with the summarized info - first iterate over bins
                "   for(c in bin_contig[b]){{" #iterate over contigs per bin
                "       for(d in bin_prot[c]){{" #iterate over different domains founded by contig
                "print b\"\\t\"c\"\\t\"contig[c]\"\\t\"d\"\\t\"bin_prot[c][d]\"\\t\"(bin_prot[c][d]/contig_num[c])*100\"\\t\"contig_num[c]" # print information
                "}} }} }} }}' | sort -k2,2 -k6g > {output}" #sort by contig and % then output

    rule mark_NOK_domain_bins_s2:
            """
            This rule take the ouput from previous rule (mark_NOK_domain_bins_s1)
            And then just compair the domain of the bin vs the domains from the contigs
            and, if the prot domain is equal or above some threshold, the contig is retained
            """
            input:
                "{run}/contig_mapping/bin_clean_tax1"
            output:
                temp("{run}/contig_mapping/bin_clean_tax2")
            shell:
                "cat {input} | awk -F \"\\t\" '{{if($3==$4 && $6>="+str(config["bin_cleaning"]["rules"]["Contig_domain_identity"])+" || $3==\"-\")"
                "{{print $1\"\\t\"$2}}}}' | sort | uniq > {output}"

    rule mark_NOK_domain_bins_s3:
            """
            Finally, we took th contigs passing the tax filter to subset the contigs failing
            """
            input:
                all_contigs="{run}/contig_mapping/bin_clean_tax1",
                passing_contigs="{run}/contig_mapping/bin_clean_tax2"
            output:
                "{run}/contig_mapping/bin_clean_tax3"
            shell:
                "cat {input.passing_contigs} | cut -f2 | grep -F -v -w -f - {input.all_contigs} "
                "| cut -f1,2 | sort | uniq > {output}"

    rule clean_bins:
            """
            Remove contigs from bins
            """
            input:
                gc_cov_nok="{run}/contig_mapping/bin_clean_1",
                taxa_nok="{run}/contig_mapping/bin_clean_tax3",
                bin=config["bin_dir"]+"{bin}."+config["bin_ext"]
            output:
                temp("{run}/clean_bins/clean.{bin}.log")
            params:
                "{run}/clean_bins"
            shell:
                #"Scripts/cleanBins.sh " +str(config["bin_dir"])+" "+str(config["bin_ext"])+ " {input.gc_cov_nok} {input.taxa_nok} {params}"
                "Scripts/cleanSingleBin.sh {input.bin} "+str(config["bin_ext"])+ "  {input.gc_cov_nok} {input.taxa_nok} {params}"

    rule concat_clean_bins_log:
            input:
                expand("{run}/clean_bins/clean.{bin}.log",  bin=config["bin_list"], run=run)
            output:
                "{run}/clean_bins/clean.log"
            shell:
                "cat {input} > {output}"

      
else:

    rule skip_clean_bins:
            output:
                "{run}/clean_bins/clean.log"
            params:
                "{run}/clean_bins"
            shell:
                "echo \"The cleaning WF was not executed.\" > {output}"

rule diamond_prots:
    input:
        "{run}/prokka/All_bins.faa" #if config["concatenate_bins"] == "T" else
        #"{run}/prokka/renamed/{bin}.faa"
         
    output:
        "{run}/diamond/All_bins.tsv"  #if config["concatenate_bins"] == "T" else
        #"{run}/diamond/{bin}.tsv"
    benchmark:
        "{run}/diamond/diamond.All_bins.benchmark" #if config["concatenate_bins"] == "T" else
        #"{run}/diamond/diamond.{bin}.benchmark"
    shell:
        "diamond blastp -q {input} --evalue "+ str(config["diamond"]["evalue"]) + " --threads "+ str(config["diamond"]["threads"]) +" "  
        " --seq "+ str(config["diamond"]["seq"]) + " --db "+ str(config["diamond"]["db"]) + " --taxonmap "+ str(config["diamond"]["taxonmap"]) +" "
        "--outfmt 6 qseqid qtitle qlen sseqid salltitles slen qstart qend sstart send evalue bitscore length pident staxids "
        "-o {output} " + str(config["diamond"]["extra_params"]) 

rule merge_diamond_results:
    input:
        #expand"({run}/diamond/{bin}.tsv"
        "{run}/diamond/All_bins.tsv"
    output:
        temp("{run}/diamond/tmp")
    shell:
        "cat {input} | awk -F'\\t' -v OFS=\"\\t\" '{{ print $1, $5, $6, $11, $12, $14, $15 }}' > {output}"

rule parse_diamond_result:
    input:
        dmd="{run}/diamond/tmp",
        prot="{run}/contig_mapping/AllProteins_list.txt"

    output:#mark as temp
        "{run}/diamond/tmp_2"
#    benchmark:
#        "{run}/diamond/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/diamond/diamond.{bin}.benchmark"
    shell:
         "python Scripts/parse_diamond_blast_results_id_taxid.py -i {input.prot} -d {input.dmd} -o {output}"
rule format_diamond_result:
    input:
        "{run}/diamond/tmp_2"
    output:#mark as temp
        "{run}/diamond/tmp_3"
#    benchmark:
#        "{run}/diamond/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/diamond/diamond.{bin}.benchmark"
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

rule merge_diamond_taxa:
    input:
        "{run}/diamond/tmp_3"
    output:#mark as temp
        "{run}/diamond/tmp_tax"
#    benchmark:
#        "{run}/diamond/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/diamond/diamond.{bin}.benchmark"
    shell:
         "LC_ALL=C join -a1 -1 5 -2 1 -e'-' -t $'\\t'  -o1.1,1.2,1.3,1.4,1.5,2.2  "
         "<(LC_ALL=C sort -k5  {input}) <(LC_ALL=C sort -k1 "+config["diamond"]["taxid_to_taxonomy"]+" ) | "
         "LC_ALL=C  sort > {output}"
rule add_header_diamond_taxa:
    input:
        "{run}/diamond/tmp_tax"
    output:
        "{run}/diamond/diamond_map.tsv"
#    benchmark:
#        "{run}/diamond/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/diamond/diamond.{bin}.benchmark"
    shell:
         "echo -e \"accession\\tDiamond_TopHit\\tE_value\\tPecID\\tTaxID\\tTaxString\" "
         "| cat - {input} > {output}"
rule combine_diamond_to_stats:
    input:
        dmd="{run}/contig_mapping/Diamond_map.tsv",
        contigs="{run}/contig_mapping/B_GenomeInfo.txt"
    output:
        "{run}/contig_mapping/Diamond_map_tmp.tsv"
#    benchmark:
#        "{run}/diamond/diamond.All_bins.benchmark" if config["concatenate_bins"] == "T" else
#        "{run}/diamond/diamond.{bin}.benchmark"
    shell:
         "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.dmd} {input.contigs} > {output}"

if config["ko_hmmr"]["annotate"] == "T":
    rule ko_hmmr:
        input:
            "{run}/prokka/All_bins.faa" if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/kfam/All_bins_ko.out.tmp" if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/kfam/{bin}.out"
        benchmark:
            "{run}/kfam/kfam.All_bins.benchmark"  if config["ko_hmmr"]["concatenate_bins"] == "T" else
            "{run}/kfam/kfam.{bin}.benchmark"
        shell:
            #hmmr header: # target name, accession, query name, accession, E-value, score, bias, E-value, score, bias, exp, reg, clu, ov, env,dom, rep, inc, description of target
            "hmmsearch --tblout  /dev/stdout -o  /dev/null --cpu "+ str(config["ko_hmmr"]["cpus"]) +" --notextw "+ str(config["ko_hmmr"]["extra_params"])+" "
            " -E " +str(config["ko_hmmr"]["evalue"]) + " " + str(config["ko_hmmr"]["hmmr_database"]) + " {input}  > {output}"
    rule select_best_ko_hmmr:
        input:
            "{run}/kfam/{bin}.out"
        output:
            "{run}/kfam/{bin}.unq.out" 
        shell:
            "cat {input} |  awk '$0 !~ \"^#\" {{print $0}}'  | "
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' | "
            " sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1 | sort -t$'\\t' -k6,6gr  > {output}"
 
    rule concat_ko_hmmr:
        input:
            expand("{run}/kfam/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/kfam/All_bins_ko.out"
        shell:
            "cat {input} | cut -f1,3,5,6 > {output}"

    rule select_best_ko_hmmr_all:
        input:
            "{run}/kfam/All_bins_ko.out.tmp"
        output:
            "{run}/kfam/All_bins_ko.out"
        shell:
            "cat {input} |  awk '$0 !~ \"^#\" {{print $0}}'  | " 
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' |  "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1  | "
            "sort -t$'\\t' -k6,6gr |  cut -f1,3,5,6  > {output}"

 

 #   rule format_ko_hmmr:
 #       input:
 #           "{run}/kfam/All_bins_pfam.out"
 #       output:
 #           temp("{run}/kfam/kfam_sorted_cols")
 #       shell:
 #           "cat {input} | cut -f1,3,5,6  > {output}"

            
    rule merge_ko_hmmr:
        input:
            ko= "{run}/kfam/All_bins_ko.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/kfam/kfam_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.ko}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick 
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
            
    rule add_names_ko_hmmr:
        input:
            ko="{run}/kfam/kfam_merged",
        output:#mark as temp
            "{run}/kfam/kfam_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,1.3,1.4,2.2,2.12 "
            "<(LC_ALL=C sort -k2 {input.ko}) "
            "<(LC_ALL=C sort -k1 "+str(config["ko_hmmr"]["database_names"])+") | LC_ALL=C  sort > {output}"
    rule identify_high_confidence_ko_hmmr:
        input:
            "{run}/kfam/kfam_names_tmp"
        output:#mark as temp
            temp("{run}/kfam/kfam_names_tmp_2")
        shell:
            "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"
               
    rule add_header_ko_hmmr:
        input:
            "{run}/kfam/kfam_names_tmp_2"
        output:
            "{run}/kfam/kfam_map.tsv"
        shell:
            "echo -e \"accession\\tKO_hmm\\te_value\\tbit_score\\tbit_score_cutoff\\tDefinition\\tconfidence\" | "
            "cat - {input} > {output}"
                    
if config["arCOG_hmmr"]["annotate"] == "T":
    rule arCOG_hmmr:
        input:
            "{run}/prokka/All_bins.faa" if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/arCOG/All_bins_arCOG.out.tmp" if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/arCOG/{bin}.out"
        benchmark:
            "{run}/arCOG/arCOG.All_bins.benchmark"  if config["arCOG_hmmr"]["concatenate_bins"] == "T" else
            "{run}/arCOG/arCOG.{bin}.benchmark"
        shell:
            "hmmsearch --tblout  /dev/stdout -o  /dev/null --cpu "+ str(config["arCOG_hmmr"]["cpus"]) +" --notextw "+ str(config["arCOG_hmmr"]["extra_params"])+" "
            " -E " +str(config["arCOG_hmmr"]["evalue"]) + " " + str(config["arCOG_hmmr"]["hmmr_database"]) + " {input}  "
            " > {output}" 
    rule select_best_arcCOG_hmmr:
        input:
            "{run}/arCOG/{bin}.out"
        output:
            "{run}/arCOG/{bin}.unq.out"
        shell:
            "cat {input} |  awk '$0 !~ \"^#\" {{print $0}}'  | "
            "perl -lane 'print join \"\t\",@F[0..17],join \" \",@F[18..$#F]' | "
            " sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1 | sort -t$'\\t' -k6,6gr > {output}"
     
    rule concat_arCOG_hmmr:
        input:
            expand("{run}/arCOG/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/arCOG/All_bins_arCOG.out"
        shell:
            "cat {input} > {output}"

    rule select_best_arCOG_hmmr_all:
        input:
            "{run}/arCOG/All_bins_arCOG.out.tmp"
        output:
            "{run}/arCOG/All_bins_arCOG.out"
        shell:
            "cat {input} |   awk '$0 !~ \"^#\" {{print $0}}' | "
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' |  "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1  | "
            "sort -t$'\\t' -k6,6gr |  cut -f1,3,5,6  > {output}"


    rule format_arCOG_hmmr:
        input:
            "{run}/arCOG/All_bins_arCOG.out"        
        output:#mark as temp
            "{run}/arCOG/arcog_tmp"
    #    benchmark:
    #        "{run}/arCOG/arCOG.All_bins.benchmark"  if config["concatenate_bins"] == "T" else
    #        "{run}/arCOG/arCOG.{bin}.benchmark"
        shell:
            "cat {input} | awk  -v OFS='\\t' '{{split($2,a,\".\"); print $1, a[1], $3,$4}}' | "
            "LC_ALL=C sort > {output}"
    
    rule merge_arCOG_hmmr:
        input:
            arc="{run}/arCOG/arcog_tmp",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/arCOG/arcog_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.arc}) | LC_ALL=C sort  > {output}"
    #this is done to format the databse names        
    rule rm_arCOG_def_spaces:
        output:
            temp("{run}/arCOG/arcog_names.tsv")
        shell:
            "cat " + str(config["arCOG_hmmr"]["database_names"])+" | sed 's/ /_/g' > {output}"
    
    rule add_names_arCOG_hmmr:
        input:
            arc="{run}/arCOG/arcog_merged",
            db="{run}/arCOG/arcog_names.tsv"
        output:#mark as temp
            temp("{run}/arCOG/arcogs_names_tmp")
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,0,2.3,2.4,2.2,1.3 "
            "<(LC_ALL=C sort -k2 {input.arc}) "
            "<(LC_ALL=C sort -k1 {input.db}) | LC_ALL=C  sort > {output}"

    rule add_header_arCOG_hmmr:
        input:
            "{run}/arCOG/arcogs_names_tmp"
        output:
            "{run}/arCOG/arcogs_map.tsv"
        shell:
            "echo -e \"accession\\tarcogs\\tarcogs_geneID\\tarcogs_Description\\tPathway\\tarcogs_evalue\" | "
            "cat - {input} > {output}"

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
    #        "{run}/arCOG/arCOG.All_bins.benchmark"  if config["concatenate_bins"] == "T" else
    #        "{run}/arCOG/arCOG.{bin}.benchmark"
        shell:
            "awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}FNR==NR{{a[$1]=$0;next}}{{print $0,a[$1]}}' {input.arc} {input.map} > {output}"  

if config["cog_diamond"]["annotate"] == "T":
     rule cog_diamond:
        input:
            "{run}/prokka/All_bins.faa" if config["cog_diamond"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"

        output:
            "{run}/cog_diamond/All_bins_cog.out.tmp"  if config["cog_diamond"]["concatenate_bins"] == "T" else
            "{run}/cog_diamond/{bin}.out"
        benchmark:
            "{run}/cog_diamond/cog_diamond.All_bins.benchmark" if config["cog_diamond"]["concatenate_bins"] == "T" else
            "{run}/cog_diamond/cog_diamond.{bin}.benchmark"
        shell:
            "diamond blastp -q {input} --evalue "+ str(config["cog_diamond"]["evalue"]) + " --threads "+ str(config["cog_diamond"]["threads"]) +" "
            " --seq "+ str(config["cog_diamond"]["seq"]) + " --db "+ str(config["cog_diamond"]["db"]) + " --taxonmap "+ str(config["cog_diamond"]["taxonmap"]) +" "
            "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
            "-o {output} " + str(config["cog_diamond"]["extra_params"])

     rule concat_cog_diamond:
        input:
            expand("{run}/cog_diamond/{bin}.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/cog_diamond/All_bins_cog.out.tmp"
        shell:
            "cat {input} > {output}"

     rule select_best_cog_diamond:
        input:
            "{run}/cog_diamond/All_bins_cog.out.tmp"
        output:
            "{run}/cog_diamond/All_bins_cog.out.tmp2"
        shell:
            "perl Scripts/best_blast.pl {input} {output}"

#     rule concat_transporter_blast_all:
#        input:
#            "{run}/cog_diamond/All_bins_cog.out.tmp2"
#        output:
#            temp("{run}/cog_diamond/All_bins_cog.out.tmp3")
#        shell:
#            "cat {input} | awk -F'\\t' -v OFS='\\t' '{{split($2,a,\"|\"); print $1, a[1], a[2], $11}}' "
#            "| LC_ALL=C sort > {output}"

     rule merge_cog_diamond:
        input:
            cogs="{run}/cog_diamond/All_bins_cog.out.tmp2",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            temp("{run}/cog_diamond/cog_merged")
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.11 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.cogs}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

     rule add_names_cog_diamond:
        input:
            "{run}/cog_diamond/cog_merged"
        output:#mark as temp
            temp("{run}/cog_diamond/cog_names_tmp")
        shell:
            "LC_ALL=C join -a1 -j2 -e'-' -t $'\\t'  -o1.1,2.3,2.5,1.3 "
            "<(LC_ALL=C sort -k2 {input}) "
            "<(LC_ALL=C sort -k2 "+str(config["cog_diamond"]["database_names"])+") | LC_ALL=C  sort > {output}"

     rule add_header_cog_diamond:
        input:
            "{run}/cog_diamond/cog_names_tmp"
        output:
            "{run}/cog_diamond/cog_map.tsv"
        shell:
            "echo -e \"accession\\tcogs\\tcogs_Description\\tcogs_evalue\" | "
            "cat - {input} > {output}"

if config["pfam_hmmr"]["annotate"] == "T":
     rule pfam_hmmr:
        input:
            "{run}/prokka/All_bins.faa" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/pfam/All_bins_pfam.out.tmp" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/pfam/{bin}.out"
        benchmark:
            "{run}/pfam/pfam.All_bins.benchmark" if config["pfam_hmmr"]["concatenate_bins"] == "T" else
            "{run}/pfam/pfam.{bin}.benchmark"
        shell:
            "hmmsearch --tblout  /dev/stdout -o  /dev/null --cpu "+ str(config["pfam_hmmr"]["cpus"]) +" --notextw "+ str(config["pfam_hmmr"]["extra_params"])+" "
            " -E " +str(config["pfam_hmmr"]["evalue"]) + " "  + str(config["pfam_hmmr"]["hmmr_database"]) + " {input}  > {output}"

     rule select_best_pfam_hmmr:
        input:
            "{run}/pfam/{bin}.out" 
        output:
            "{run}/pfam/{bin}.unq.out"
        shell:
            "cat {input} |   awk '$0 !~ \"^#\" {{print $0}}' | perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' | "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1 | sort -t$'\\t' -k6,6gr > {output}"

     rule select_best_pfam_hmmr_all:
        input:
            "{run}/pfam/All_bins_pfam.out.tmp"
        output:
            "{run}/pfam/All_bins_pfam.out"
        shell:
            "cat {input} |   awk '$0 !~ \"^#\" {{print $0}}' | "
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' |  "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1  | "
            "sort -t$'\\t' -k6,6gr |  cut -f1,3,5,6  > {output}"

     rule concat_pfam_hmmr:
        input:
            expand("{run}/pfam/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/pfam/All_bins_pfam.out"
        shell:
            "cat {input} |  cut  -f1,3,5,6  > {output}"
        
#     rule format_pfam_hmmr:
#        input:
#            "{run}/pfam/All_bins_pfam.out"
#        output:
#            temp("{run}/pfam/pfam_sorted_cols")
#        shell:
#            "cat {input} | cut  -f1,3,5,6 | sed 's/ /\t/g' > {output}"
    
     rule merge_pfam_hmmr:
        input:
            pfam="{run}/pfam/All_bins_pfam.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/pfam/pfam_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.pfam}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick 
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
            
     rule add_names_pfam_hmmr:
        input:
            pfam="{run}/pfam/pfam_merged",
        output:#mark as temp
            "{run}/pfam/pfam_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 4 -e'-' -t $'\\t' -o1.1,2.1,0,2.5,1.3,1.4 "
            "<(LC_ALL=C sort -k2 {input.pfam}) "
            "<(LC_ALL=C sort -k4 "+str(config["pfam_hmmr"]["database_names"])+") | LC_ALL=C  sort > {output}"
 #    rule identify_high_confidence_pfam_hmmr:
 #       input:
 #           "{run}/pfam/pfam_names_tmp"
 #       output:#mark as temp
 #           temp("{run}/pfam/pfam_names_tmp_2")
 #       shell:
 #           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"
               
     rule add_header_pfam_hmmr:
        input:
            "{run}/pfam/pfam_names_tmp"
        output:
            "{run}/pfam/pfam_map.tsv"
        shell:
            "echo -e \"accession\\tPFAM_ID\\tPFAM_hmm\\tPFAM_description\\tPfam_Evalue\\tPfam_Score\" | "
            "cat - {input} > {output}"
            
if config["tigr_hmmr"]["annotate"] == "T":
    rule tigr_hmmr:
        input:
            "{run}/prokka/All_bins.faa" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/tigr/All_bins_tigr.out.tmp" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/tigr/{bin}.out"
        benchmark:
            "{run}/tigr/tigr.All_bins.benchmark" if config["tigr_hmmr"]["concatenate_bins"] == "T" else
            "{run}/tigr/tigr.{bin}.benchmark"
        shell:
            "hmmsearch --tblout  /dev/stdout -o  /dev/null --cpu "+ str(config["tigr_hmmr"]["cpus"]) +" --notextw "+ str(config["tigr_hmmr"]["extra_params"])+" "
            " -E " +str(config["tigr_hmmr"]["evalue"]) + " "  + str(config["tigr_hmmr"]["hmmr_database"]) + " {input} > {output}"

    rule select_best_tigr_hmmr:
        input:
            "{run}/tigr/{bin}.out"
        output:
            "{run}/tigr/{bin}.unq.out"
        shell:
            "cat {input} |   awk '$0 !~ \"^#\" {{print $0}}' | perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' | "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1 | sort -t$'\\t' -k6,6gr > {output}"

    rule concat_tigr_hmmr:
        input:
            expand("{run}/tigr/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/tigr/All_bins_tigr.out"
        shell:
            "cat {input} |  cut  -f1,3,5,6  > {output}"

    rule select_best_tigr_hmmr_all:
        input:
            "{run}/tigr/All_bins_tigr.out.tmp"
        output:
            "{run}/tigr/All_bins_tigr.out"
        shell:
            "cat {input} |   awk '$0 !~ \"^#\" {{print $0}}' | "
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' |  "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1  | "
            "sort -t$'\\t' -k6,6gr |  cut -f1,3,5,6  > {output}"
    rule merge_tigr_hmmr:
        input:
            tigr="{run}/tigr/All_bins_tigr.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/tigr/tigr_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.tigr}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

    rule add_names_tigr_hmmr:
        input:
            tigr="{run}/tigr/tigr_merged",
        output:#mark as temp
            "{run}/tigr/tigr_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,2.2,2.3,2.4,1.3,1.4 "
            "<(LC_ALL=C sort -k2 {input.tigr}) "
            "<(LC_ALL=C sort -k1 "+str(config["tigr_hmmr"]["database_names"])+") | LC_ALL=C  sort > {output}"
 #    rule identify_high_confidence_tigr_hmmr:
 #       input:
 #           "{run}/tigr/tigr_names_tmp"
 #       output:#mark as temp
 #           temp("{run}/tigr/tigr_names_tmp_2")
 #       shell:
 #           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"

    rule add_header_tigr_hmmr:
        input:
            "{run}/tigr/tigr_names_tmp"
        output:
            "{run}/tigr/tigr_map.tsv"
        shell:
            "echo -e \"accession\\tTIGR_hmm\\tTIGR_name\\tTIGR_description\\tTIGR_EC\\tTIGR_Evalue\\tTIGR_Score\" | "
            "cat - {input} > {output}"


if config["cazy_hmmr"]["annotate"] == "T":
     rule cazy_hmmr:
        input:
            "{run}/prokka/All_bins.faa" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/cazy/All_bins_cazy.out.tmp" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/cazy/{bin}.out"
        benchmark:
            "{run}/cazy/cazy.All_bins.benchmark" if config["cazy_hmmr"]["concatenate_bins"] == "T" else
            "{run}/cazy/cazy.{bin}.benchmark"
        shell:
            "hmmsearch --tblout  /dev/stdout -o  /dev/null --cpu "+ str(config["cazy_hmmr"]["cpus"]) +" --notextw "+ str(config["cazy_hmmr"]["extra_params"])+" "
            " -E " +str(config["cazy_hmmr"]["evalue"]) + " "  + str(config["cazy_hmmr"]["hmmr_database"]) + " {input}  > {output}"

     rule select_best_cazy_hmmr:
        input:
            "{run}/cazy/{bin}.out"
        output:
            "{run}/cazy/{bin}.unq.out"
        shell:
            "cat {input} |  awk '$0 !~ \"^#\" {{print $0}}'  | perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' | "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1 | sort -t$'\\t' -k6,6gr > {output}"

     rule concat_cazy_hmmr:
        input:
            expand("{run}/cazy/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/cazy/All_bins_cazy.out"
        shell:
            "cat {input} |  cut  -f1,3,5,6 | sed 's/\.hmm//'  > {output}"

     rule select_best_cazy_hmmr_all:
        input:
            "{run}/cazy/All_bins_cazy.out.tmp"
        output:
            "{run}/cazy/All_bins_cazy.out"
        shell:
            "cat {input} |  awk '$0 !~ \"^#\" {{print $0}}'  | "
            "perl -lane 'print join \"\\t\",@F[0..17],join \" \",@F[18..$#F]' |  "
            "sort -t$'\\t' -k6,6gr  | sort -t$'\\t' --stable -u -k1,1  | "
            "sort -t$'\\t' -k6,6gr |  cut -f1,3,5,6 | sed 's/\.hmm//'  > {output}"
     rule merge_cazy_hmmr:
        input:
            cazy="{run}/cazy/All_bins_cazy.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/cazy/cazy_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.cazy}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
     rule add_names_cazy_hmmr:
        input:
            cazy="{run}/cazy/cazy_merged",
        output:#mark as temp
            "{run}/cazy/cazy_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,2.2,1.3,1.4 "
            "<(LC_ALL=C sort -k2 {input.cazy}) "
            "<(LC_ALL=C sort -k1 "+str(config["cazy_hmmr"]["database_names"])+") | LC_ALL=C  sort > {output}"
 #    rule identify_high_confidence_cazy_hmmr:
 #       input:
 #           "{run}/cazy/cazy_names_tmp"
 #       output:#mark as temp
 #           temp("{run}/cazy/cazy_names_tmp_2")
 #       shell:
 #           "cat {input} | awk  -v OFS='\\t' '{{ if ($4 > $5){{ $7=\"high_score\" }}else{{ $7=\"-\" }} print }} ' > {output}"

     rule add_header_cazy_hmmr:
        input:
            "{run}/cazy/cazy_names_tmp"
        output:
            "{run}/cazy/cazy_map.tsv"
        shell:
            "echo -e \"accession\\tCAZY_hmm\\tCAZY_description\\tCAZY_Evalue\\tCAZY_Score\" | "
            "cat - {input} > {output}"

if config["merops_blast"]["annotate"] == "T":
     rule merops_blast:
        input:
            "{run}/prokka/All_bins.faa" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/merops/All_bins_merops.out.tmp" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.out"
        benchmark:
            "{run}/merops/merops.All_bins.benchmark" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/merops.{bin}.benchmark"
        shell:
            "blastp -num_threads  "+ str(config["merops_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["merops_blast"]["blast_database"]) + " "
            " -out {output} -evalue " +str(config["merops_blast"]["evalue"]) + " "+str(config["merops_blast"]["extra_params"])

     rule select_best_merops_blast:
        input:
            "{run}/merops/All_bins_merops.out.tmp" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.out"
        output:
            "{run}/merops/All_bins_merops.out" if config["merops_blast"]["concatenate_bins"] == "T" else
            "{run}/merops/{bin}.unq.out"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

     rule concat_merops_blast:
        input:
            expand("{run}/merops/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/merops/All_bins_merops.out"
        shell:
            "cat {input}  > {output}"
            
     rule merge_merops_blast:
        input:
            merops="{run}/merops/All_bins_merops.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/merops/merops_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.11 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.merops}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

     rule add_names_merops_blast:
        input:
            merops="{run}/merops/merops_merged"
        output:#mark as temp
            "{run}/merops/merops_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,1.2,2.2,1.3,1.4 "
            "<(LC_ALL=C sort -k2 {input.merops}) "
            "<(LC_ALL=C sort -k1 "+str(config["merops_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

     rule add_header_merops_blast:
        input:
            "{run}/merops/merops_names_tmp"
        output:
            "{run}/merops/merops_map.tsv"
        shell:
            "echo -e \"accession\\tMEROPS_blast\\tMEROPS_description\\tMEROPS_Idenity\\tMEROPS_Evalue\" | "
            "cat - {input} > {output}"    

if config["transporterDB_blast"]["annotate"] == "T":
     rule transporterDB_blast:
        input:
            "{run}/prokka/All_bins.faa" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/transporter/All_bins_transporter.out.tmp" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/{bin}.out"
        benchmark:
            "{run}/transporter/transporter.All_bins.benchmark" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/transporter.{bin}.benchmark"
        shell:
            "blastp -num_threads  "+ str(config["transporterDB_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["transporterDB_blast"]["blast_database"]) + " "
            " -out {output} -evalue " +str(config["transporterDB_blast"]["evalue"]) + " "+str(config["transporterDB_blast"]["extra_params"])

     rule select_best_transporter_blast:
        input:
            "{run}/transporter/All_bins_transporter.out.tmp" if config["transporterDB_blast"]["concatenate_bins"] == "T" else 
            "{run}/transporter/{bin}.out"
        output:
            "{run}/transporter/All_bins_transporter.out.tmp2" if config["transporterDB_blast"]["concatenate_bins"] == "T" else
            "{run}/transporter/{bin}.unq.out"
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

     rule concat_transporter_blast:
        input:
            expand("{run}/transporter/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/transporter/All_bins_transporter.out"
        shell:
            "cat {input} | awk -F'\\t' -v OFS='\\t' '{{split($2,a,\"|\"); print $1, a[1], a[2], $11}}' "
            "| LC_ALL=C sort > {output}"
     rule concat_transporter_blast_all:
        input:
            "{run}/transporter/All_bins_transporter.out.tmp2"
        output:
            "{run}/transporter/All_bins_transporter.out"
        shell:
            "cat {input} | awk -F'\\t' -v OFS='\\t' '{{split($2,a,\"|\"); print $1, a[1], a[2], $11}}' "
            "| LC_ALL=C sort > {output}"
            
     rule merge_transporter_blast:
        input:
            transporter="{run}/transporter/All_bins_transporter.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/transporter/transporter_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.3,2.4 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.transporter}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"

     rule add_names_transporter_blast:
        input:
            transporter="{run}/transporter/transporter_merged",
        output:#mark as temp
            "{run}/transporter/transporter_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 2 -e'-' -t $'\\t' -o1.1,0,2.1,2.3,1.3 "
            "<(LC_ALL=C sort -k2 {input.transporter}) "
            "<(LC_ALL=C sort -k2 "+str(config["transporterDB_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

     rule add_header_transporter_blast:
        input:
            "{run}/transporter/transporter_names_tmp"
        output:
            "{run}/transporter/transporter_map.tsv"
        shell:
            "echo -e \"accession\\tTransporter_id\\tTransporter_id2\\tTransporter_description\\tTransporter_Evalue\" | "
            "cat - {input} > {output}"

if config["hydDB_blast"]["annotate"] == "T":
     rule hydDB_blast:
        input:
            "{run}/prokka/All_bins.faa" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/prokka/renamed/{bin}.faa"
        output:
            "{run}/hydDB/All_bins_hydDB.out.tmp" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/hydDB/{bin}.out"
        benchmark:
            "{run}/hydDB/hydDB.All_bins.benchmark" if config["hydDB_blast"]["concatenate_bins"] == "T" else
            "{run}/hydDB/hydDB.{bin}.benchmark"
        shell:
            "blastp -num_threads  "+ str(config["hydDB_blast"]["threads"]) +" -outfmt 6 -query {input} -db "  + str(config["hydDB_blast"]["blast_database"]) + " "
            " -out {output} -evalue " +str(config["hydDB_blast"]["evalue"]) + " "+str(config["hydDB_blast"]["extra_params"])

     rule select_best_hydDB_blast:
        input:
            "{run}/hydDB/All_bins_hydDB.out.tmp" if config["hydDB_blast"]["concatenate_bins"] == "T" else 
            "{run}/hydDB/{bin}.out"
        output:
            "{run}/hydDB/All_bins_hydDB.out" if config["hydDB_blast"]["concatenate_bins"] == "T" else 
            "{run}/hydDB/{bin}.unq.out" 
        shell:   #try to replace with sort
            "perl Scripts/best_blast.pl {input} {output}"

     rule concat_hydDB_blast:
        input:
            expand("{run}/hydDB/{bin}.unq.out",  bin=config["bin_list"], run=run)
        output:
            "{run}/hydDB/All_bins_hydDB.out"
        shell:
            "cat {input}  > {output}"
            
     rule merge_hydDB_blast:
        input:
            hydDB="{run}/hydDB/All_bins_hydDB.out",
            prots="{run}/contig_mapping/1_Bins_to_protein_list.txt"
        output:#mark as temp
            "{run}/hydDB/hydDB_merged"
        shell:
            "LC_ALL=C join -a1 -j1 -e'-' -t $'\\t' -o 0,2.2,2.3,2.11,2.12 "
            "<(LC_ALL=C sort {input.prots}) "
            "<(LC_ALL=C sort {input.hydDB}) | LC_ALL=C sort  "
            #get rid of empty space this was in the past the -e "-" was missing and now it makes the trick
            #"awk 'BEGIN {{FS = OFS = \"\\t\"}} {{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"-\" }}; 1' "
            "> {output}"
    
     rule add_names_hydDB_blast:
        input:
            hydDB="{run}/hydDB/hydDB_merged",
        output:#mark as temp
            "{run}/hydDB/hydDB_names_tmp"
        shell:
            "LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\\t' -o1.1,0,2.2,1.4 "
            "<(LC_ALL=C sort -k2 {input.hydDB}) "
            "<(LC_ALL=C sort -k1 "+str(config["hydDB_blast"]["database_names"])+") | LC_ALL=C  sort > {output}"

     rule add_header_hydDB_blast:
        input:
            "{run}/hydDB/hydDB_names_tmp"
        output:
            "{run}/hydDB/hydDB_map.tsv"
        shell:
            "echo -e \"accession\\tHydDB\\tDescription\\tHydDB_evalue\" | "
            "cat - {input} > {output}"

rule summarize_annotation:
    input:
        "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/diamond/diamond_map.tsv",
        "{run}/arCOG/arcogs_map.tsv" if config["arCOG_hmmr"]["annotate"] == "T" else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/kfam/kfam_map.tsv" if config["ko_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/pfam/pfam_map.tsv" if config["pfam_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/tigr/tigr_map.tsv" if config["tigr_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/cazy/cazy_map.tsv" if config["cazy_hmmr"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/cog_diamond/cog_map.tsv" if config["cog_diamond"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/merops/merops_map.tsv" if config["merops_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/transporter/transporter_map.tsv" if config["transporterDB_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt",
        "{run}/hydDB/hydDB_map.tsv" if config["hydDB_blast"]["annotate"] == "T"  else "{run}/contig_mapping/B_GenomeInfo.txt"
    output:
        "{run}/FunctionalAnnotation.tsv"
    params:
        run="{run}"
    shell:
        "Scripts/summary_annotation.sh {params.run} {input[0]}  {output}"

rule report:
    input:
        "{run}/FunctionalAnnotation.tsv",
        "{run}/clean_bins/clean.log"
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
        #"{run}/prokka/{bin}.gbk.gz",
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
    shell:
        "touch {output}"

