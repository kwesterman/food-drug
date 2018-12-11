#!/usr/bin/python
"""Search database A for all compounds with a Tanimoto cutoff of 0.85 compared to compounds in database B"""

import urllib
import xml.etree.cElementTree as ET
import sys

def sim_search (url,dbB_id): 
	"""Take in a search url and the ID of the origin compound (from dbB), return a dict of similar compound info (from dbA)"""
	dbA_id, listA_structures, listA_names = ([] for i in xrange(3))
	searchResults = urllib.urlopen(url)
 	try: #To avoid the few problematic XML files returned
 		tree = ET.parse(searchResults) #Parse XML file returned from database A API
 	except:
 		pass
	else: #Append all similar list A compound info to growing info lists for input list B compound
		for elem in tree.findall('.//molecule/molecule_chembl_id'):
			dbA_id.append(elem.text)
		for elem in tree.findall('.//canonical_smiles'):
			listA_structures.append(elem.text)
		for elem in tree.findall('.//pref_name'):
			listA_names.append(elem.text)	
	return {'dbA_id': dbA_id, 'listA_structures': listA_structures, 'listA_names': listA_names, 'dbB_id': [dbB_id]*len(dbA_id)}

if __name__ == '__main__':
	if len(sys.argv) != 2:
		exit('Usage: ./foodcmpds_to_drugs.py thread_number')
	task_id = int(sys.argv[1])

	with open('foodcmpds_name.txt') as file: #Create lists of compound names and structures
		listB_names = file.read().split('\n')
	with open('foodcmpds_smiles.txt') as file:
		listB_structures = file.read().split('\n')

	results = {'dbA_id': [], 'listA_structures': [], 'listA_names': [], 'dbB_id': []}

	for i in xrange(500*(task_id-1), 500*task_id): #Perform sim. searches on a 500 cmpd. "chunk" of much larger input compound list
		if i < len(listB_names):
			url = 'https://www.ebi.ac.uk/chembl/api/data/similarity/%s/85?limit=1000' % listB_structures[i] #Create query URL
			for key, value in sim_search(url, listB_names[i]).items(): #Append results of most recent cmpd. search
				results[key].extend(value)

	with open('out%s' % sys.argv[1], 'w') as outfile: #Write results to a 4-column tab-sep. file
		for i in xrange(len(results['dbA_id'])):
			outfile.write('%s\t%s\t%s\t%s\n' % (results['dbB_id'][i],results['dbA_id'][i],results['listA_names'][i],results['listA_structures'][i]))
				

#SHELL FILE USED TO CALL THIS ONE:
# #!/bin/bash
# #SBATCH -o %A_%a.out
# #SBATCH -e %A_%a.err
# #SBATCH -c 1
# #SBATCH -p batch

# ./compound_similarity.py $SLURM_ARRAY_TASK_ID


#TO CALL ENTIRE JOB USING SLURM:
#sbatch -N 3 --array=1-50 ./run_sim_search.sh