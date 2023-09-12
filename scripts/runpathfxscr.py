
##################################################

print('##################################################')
print('')

'''
Run PathFX for drugs included:
- in a '.txt' file
- or in a list
'''

# run PathFX

import os
import pandas as pd

main_path = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
os.chdir(main_path)

#drugs_ls = pd.read_csv('all_890_drugs.txt', sep='\t')
#dl = drugs_ls['Drugs'].tolist()
#print('Made the list of drugs!')
#print('')

dl = ['alteplase', 'atropine']

os.chdir(os.path.join(main_path, 'PathFX/scripts/'))
for drug_name in dl:
  cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s' %(drug_name, drug_name)
  #cmd = 'python phenotype_enrichment_pathway_Pfxevalall4SEs1.py -d %s -a %s' %(drug_name, drug_name)
  os.system(cmd)
print('Ran PathFX for all drugs!')

print('')
print('##################################################')

##################################################
