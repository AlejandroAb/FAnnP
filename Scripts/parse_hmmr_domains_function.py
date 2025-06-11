#!/usr/bin/python
# Tijn 2024-2025 :)
# Last recorded update: 02-2025
#
#

# Libraries
import re


# Put in a function for pytest, called by parse_hmmr_domains.py
def parse_hmmr_domains(input_file, percentage_overlap_allowed, cutoff_value, domain_index_value, 
                        output_file_normal, create_output_file_sql, output_file_sql):
    domain_list = []
    feature_list = []
    approved_domains = []
    no_overlap = True
    unique_domain = True
    low_evalue = True
    count = 0

    # Open sql specific output file if we want to make a database.
    if create_output_file_sql:
        output_handle2 = open(output_file_sql, "w")

    with open(input_file) as input_handle1, open(output_file_normal, "w") as output_handle1:
        for line in input_handle1:
            # Debug gene entry limit.
            #count += 1
            #if count == 2000:
            #	break
            # Add domain with features in the gene to a list with domains.
            feature_list.append(line.strip("\n"))
            feature_list = re.sub(r" +", r"\t", feature_list[0]).split("\t")
            domain_list.append(feature_list)
            # If it is the last found domain in the gene. Go through the domain list.
            if int(feature_list[1]) == len(domain_list):
                # Sorts by independant domain score (14), overall domain score = (8).
                domain_list.sort(key = lambda row: float(row[14]), reverse=True)
                for domain_entry in domain_list:
                    # First domain will have the highest score within the gene and will always get accepted.
                    if not approved_domains and float(cutoff_value) > float(domain_entry[13]):
                        approved_domains.append(domain_entry)
                    else:
                        # Checks if the domain has a low independent evalue.
                        # Doesn't currently check if the first domain has a low independent evalue.
                        # Default evalue cutoff will be 1e-4.
                        if float(cutoff_value) < float(domain_entry[13]):
                            low_evalue = False
                        # Compare every approved domain to the incoming candidates.
                        for approved_entry in approved_domains:
                            # Here it is calculated how much overlap we will allow domains to have.
                            # Calculates the amount of overlap which is by default 25% of smallest domain.
                            # (20) = begin domain position, (21) = end domain position
                            percentage_list = sorted([float(domain_entry[21]) - float(domain_entry[20]), float(approved_entry[21]) - float(approved_entry[20])])
                            overlap_amount = int(round(percentage_list[0] / (100 / int(percentage_overlap_allowed))))
                            # Checks if the domain wasn't found already in another domain entry.
                            #if approved_entry[4] == domain_entry[4]:
                            #	unique_domain = False
                            # Checks wether or not there is any overlap between the domains.
                            if float(approved_entry[20]) + overlap_amount < float(domain_entry[21]) and float(approved_entry[21]) - overlap_amount > float(domain_entry[20]):
                                no_overlap = False
                        # Add the domain to the approved domain list within the gene.
                        if no_overlap and unique_domain and low_evalue:
                            approved_domains.append(domain_entry)
                        low_evalue = True
                        unique_domain = True
                        no_overlap = True
                # Write out all approved domain entries of a gene.
                # Using --domtblout (0) = gene name, (4) = domain name, (7) = overall domain evalue, (8) = overall domain score,
                # (13) = independant domain e-value, (14) independant domain score
                # Write to sql specific output file if we want to make a database.
                if approved_domains:
                    if create_output_file_sql:
                        for approved_entry in approved_domains:
                            # Debug output
                            #output_handle2.write(approved_entry[0] + "\t" + approved_entry[int(domain_index_value)] + "\ti-E-val: " + approved_entry[13] + "\tScore: " + approved_entry[14] + "\tBegin_pos: " + approved_entry[20] + "\tEnd_pos: " + approved_entry[21] + "\tTot: " + approved_entry[1] + "\t" + str(len(approved_domains)) + "\n")
                            # Regular output (Gene - Domain - Domain Begin - Domain End - I-Evalue - I-Score)
                            output_handle2.write(approved_entry[0] + "\t" + approved_entry[int(domain_index_value)] + "\t" + approved_entry[20] + "\t" + approved_entry[21] + "\t" + approved_entry[13] + "\t" + approved_entry[14] + "\n")
                    # Merged output
                    domain_name = ','.join([row[int(domain_index_value)] for row in approved_domains])
                    domain_evalue = ','.join([row[13] for row in approved_domains])
                    domain_score = ','.join([row[14] for row in approved_domains])
                    domain_coverage = ','.join([row[22] for row in approved_domains])
                    output_handle1.write(approved_domains[0][0] + "\t" + domain_name + "\t" + domain_evalue + "\t" + domain_score + "\t" + domain_coverage + "\n")
                # Empty the lists containing domains for the next gene.
                approved_domains = []
                domain_list = []
            # Empty domain entry list for next found domain.
            feature_list = []
    # Close sql specific output file if we want to make a database.
    if create_output_file_sql:
        output_handle2.close()

