This project aims to discover novel bioactivity of foods and food compounds through chemical similarity tests of these food compounds with large-scale biological and pharmaceutical databases. It uses the basic QSAR (quantitative structure activity relationship) model in which chemical similarity between two compounds suggests similarity in bioactivity. Thus, drugs of interest, or those targeting proteins of interest, can be compared to a database of food compounds to suggest novel small molecule-based dietary interventions for specific phenotypes.

Components:
- ChemmineR_manipulation.R/compound_similarity.py: These R and Python scripts are slightly different incarnations of the same basic process: taking in compound structure data and performing chemical similarity searches.
- find_foods.txt/get_foodcmpds_info.txt/get_foodcmpds_smiles.txt: SQL queries for retrieving data from FooDB (food- and food compound-related database). 
- chembl_getDrugs.txt: SQL query for retrieving drugs with structures and metadata from the ChEMBL database.
- OpenBabel_commands: Miscellaneous commands for use with OpenBabel command-line tool.
- ppi_target_prioritization.R: R script used to prioritize gene targets. A PPI network is constructed based on a set of genes involved in triglyceride metabolism. Clustering is then performed to detect functional pathways to be targeted by nutritional interventions. 
