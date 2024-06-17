#!/usr/bin/python3
# -*-coding:Latin-1 -*

'''This script takes as input a blast output formatted as the exemplar (csv) dataset given :
it uses blast alignment score and compare it between all given alignment and find 
the Lowest Common Ancestor; one per query sequence.

usage : python3 from from_csv_to_lca.py input.csv > output.tsv

@author: Tirera Sourakhata <stirera(at)pasteur-cayenne.fr>
'''

import sys
from operator import itemgetter
from collections import defaultdict
from itertools import groupby
from ete3 import NCBITaxa, Tree
import csv

def transform_tabulated_file(file_path, delimiter='\t'):
    results2 = {}
    with open(file_path, 'r') as file:
        for line in file:
            values = line.strip().split(delimiter)
            sequence_id = values[0]
            other_cols = values[1:10]
            alignment_score = float(values[11])  # Assuming alignment score is an integer
            other_cols.append(alignment_score)
            other_cols.extend(values[12:])
            #other_cols.append(values[12])

            entry = {'sequence_id': sequence_id, 'other_cols': other_cols}

            if entry['sequence_id'] in results2:
                results2[sequence_id].append(entry['other_cols'])
            else:
                results2[sequence_id] = [other_cols]
    
    # Convert the dictionary to a list of dictionaries

    results_list = [{'sequence_id': key, 'other_cols': value} for key, value, in results2.items()]
    
    #print ("results_list2", results_list)

    return results_list

def choose_top_scoring_hits(results):
    '''Try to keep the strucuture results at the end
        if key in bests, just append hit into the 'other_columns'
        else create and append entry

    '''
    bests = []

    for res in results:
        best_scores = []
        
        for hit in res['other_cols']:
            if len(best_scores) > 0:
                if hit[9] > max(best_scores):
                    best_scores = [hit[9]]
                    entry = {'sequence_id': res['sequence_id'], 'other_cols': hit}
                    bests = [entry]
                elif hit[9] == max(best_scores):
                    best_scores.append(hit[9])
                    bests.append({'sequence_id': res['sequence_id'], 'other_cols': hit})
            else:
                entry = {'sequence_id': res['sequence_id'], 'other_cols': hit}
                best_scores.append(hit[9])
                bests.append(entry)

    #print("This is the top hits: ", bests)
    #print("End of section\n\n")
    # Return the top-scoring hits
    #return top_hits[:]
    #bests2 = [{'sequence_id': key, 'other_cols': value} 
    bests2=[]
    sorted_results = sorted(bests, key=lambda x: x['sequence_id'])
    '''
    for hit in bests:
        if str(hit['sequence_id']) in bests2.keys():
            #print("####before",bests2['other_cols'])
            bests2['other_cols'].append(hit['other_cols'])
            #print("####after",bests2['other_cols'])

        else:
            entry = {'sequence_id' : hit['sequence_id'] , 'other_cols': hit['other_cols']}
            bests2.append(entry)                    
            #print("#####afterelse",bests2)
    '''
    grouped_results = {}
    for key, group in groupby(sorted_results, key=lambda x: x['sequence_id']):
        grouped_results[key] = list(group)
    
    #print("grouped_results=====>",grouped_results)

    return grouped_results

def read_blast_output(blast_output_file):
    blast_results = defaultdict(list)
    with open(blast_output_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            query_id = fields[0]
            tax_id = fields[12]
            blast_results[query_id].append(tax_id)
    return blast_results

#def find_lca(tax_ids):
    # Implement your LCA algorithm here
    # This function should find the Lowest Common Ancestor (LCA) from the taxonomic IDs
    # You can use a database or pre-built tree structure to determine the LCA
    
    # Placeholder implementation: returning the first taxonomic ID as the LCA
    # return tax_ids[0]

def read_taxonomic_database():
    taxonomic_db = {}
    nodes_file = './new_taxdump/nodes.dmp' 
    names_file = './new_taxdump/names.dmp'

    # Read the taxonomic nodes file
    with open(nodes_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t|\t')
            tax_id = fields[0]
            parent_tax_id = fields[1]
            rank = fields[2]
            taxonomic_db[tax_id] = {
                'parent': parent_tax_id,
                'rank': rank
            }

    # Read the taxonomic names file
    with open(names_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t|\t')
            tax_id = fields[0]
            name = fields[1]
            name_class = fields[3]
            if name_class == 'scientific name':
                taxonomic_db[tax_id]['name'] = name

    return taxonomic_db

def find_lca_0(tax_ids, taxonomic_db):
    # Find the LCA from the taxonomic IDs using the taxonomic database

    # Step 1: Get the lineage for each taxonomic ID
    lineages = []
    for tax_id in tax_ids:
        lineage = []
        while tax_id != '1':  # '1' represents the root of the taxonomic tree
            if tax_id in taxonomic_db:
                lineage.append(tax_id)
                tax_id = taxonomic_db[tax_id]['parent']
            else:
                break
        lineages.append(lineage)

    # Step 2: Find the LCA by comparing the lineages
    lca_tax_id = '1'  # Assume the root as the initial LCA
    for i in range(len(lineages[0])):
        tax_id = lineages[0][i]
        found_in_all = all(tax_id in lineage for lineage in lineages)
        if found_in_all:
            lca_tax_id = tax_id
        else:
            break

    return lca_tax_id

def store_lineages_in_tree(taxids):
    ncbi = NCBITaxa()
    root_taxid = find_lca(taxids)
    
    tree = Tree(name=str(root_taxid))
    
    for taxid in taxids:
        lineage = ncbi.get_lineage(taxid)
        lineage_names = ncbi.get_taxid_translator(lineage)
        current_node = tree
        for taxid in lineage[1:]:
            name = lineage_names[taxid]
            child_node = current_node.add_child(name=name, dist=1.0)
            current_node = child_node
    
    return tree


def find_lca_1(taxids):
    
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxids[0])
    for taxid in taxids[1:]:
        lineage = ncbi.get_common_ancestor(lineage[-1], taxid)
    return lineage[-1]

def find_lca(taxids):
    ncbi = NCBITaxa()
    topology = ncbi.get_topology(taxids)
    lca = topology.get_tree_root().name
    return int(lca)

def assign_taxonomic_labels(blast_results):
    taxonomic_labels = {}
    tax_db=read_taxonomic_database()
    for query_id, hits in blast_results.items():
        tax_ids=[]
        for hit in hits:
            tax_id=hit['other_cols'][10]
            tax_ids.append(tax_id)

        lca = find_lca(tax_ids)
        #print("Lowest Common Ancestor (LCA) TaxID:", lca)
        taxonomic_labels[query_id] = lca

    return taxonomic_labels

def taxid_to_lineage(lca_taxid) :   
    #for key, value in taxonomic_labels_top.items():
    #    outline=""
    #    print("key_label_items = ",key,"\tvalue_label_items:",value)
    taxon = NCBITaxa()
    lineage = taxon.get_lineage(str(lca_taxid))
    names = taxon.get_taxid_translator(lineage)
    lineage_names = [names[taxid] for taxid in lineage]
        #print(" lineage_names = ", lineage_names)

    #    outline+=key+"\t"+str(value)+"\t"
        #print("key;",key,";value;",value)
    outline=""
    for ln in lineage_names : 
            #print(ln,";")
    
        outline+=ln+";"
        #print(outline)
        #print('childs===>',top_hits[key])
        #childs=top_hits[key]
        #tab_out=[]
        """
        for c in childs:
            tab_out.append(key)
            tab_out.extend(c['other_cols'])
            
            outline+="|\t"+c['other_cols'][-1]
        """
            #outline+="|\t"+c['other_cols'][]
        #print_table_to_csv(tab_out)
        #print("tab_out===>", tab_out)
    #print(outline)
    return outline
        #print("other_cols",c['other_cols'])  

def print_table_to_csv(table_data):
    writer = csv.writer(sys.stdout, delimiter=';', quoting=csv.QUOTE_MINIMAL)
    for row in table_data:
        writer.writerow(row)


if __name__ == '__main__':
    
    blast_output=sys.argv[1]
    #data=transform_tabulated_file(blast_output)
    
    data=transform_tabulated_file(blast_output)
    #for result in data2:
    #    #print("#", result)
    #    #print("##", result['other_cols'])

    top_hits=choose_top_scoring_hits(data)
    #print("**top_hits==>", top_hits) 
    taxonomic_labels_top=assign_taxonomic_labels(top_hits)


    for key, value in top_hits.items():
        hitline=""
        k_test=taxonomic_labels_top[key]
        #print("k_test==", k_test)
        tax_lineage=taxid_to_lineage(k_test)
        for v in value :
            hitline+=v['sequence_id']+"\t"
            c=0
            for x in v['other_cols']:
                c+=1
                hitline+=str(x)+"\t"
                if c==15:
                    hitline+=tax_lineage+"\n"
        print(hitline,"\n")
        #print("hitline===>", hitline)
        #print("tax_lineage===>", tax_lineage)


    #taxonomic_labels_all=assign_taxonomic_labels(top_hits)
    #taxonomic_labels_top=assign_taxonomic_labels(top_hits)
    #print("**labels==>", taxonomic_labels_top)
    #|cut -f1-14,17-
