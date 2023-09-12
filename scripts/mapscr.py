
##################################################

'''
Import libraries in Python
'''

# libraries

import os
import pickle
import pandas as pd

##################################################

print('##################################################')
print('')

'''
Map the drug names in the drug toxicity data to the PathFX datbase
'''

print('The drug mapping process:')
print('')

# change the directory

os.chdir('/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/PathFX/rscs/')

# read the drugbankid_to_name file and collect the list of drug names

db2n = pickle.load(open(('drugbankid_to_name.pkl'),'rb'))
drugbank_names = []
for name, value in db2n.items():
  drugbank_names.append(value.lower())
print('# of drugs in the PathFX drugbank:', len(drugbank_names))
print('')

# change the directory

os.chdir('/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/')

# read the drug toxicity dataset and collect the list of drug names in the drug toxicity data

df = pd.read_csv('Drugs_labeled_for_AEs.txt', sep='\t', low_memory=False);
dfse = df.fillna('nodrug')
dfls = dfse.values.tolist()
df_tox = []
for i in range(len(dfls)):
  df_tox.extend(dfls[i])
drug_tox = list(set(df_tox))
drug_tox.remove('nodrug')
drug_tox_ls = []
for drugname in drug_tox:
  drug_tox_ls.append(drugname.lower())
print('# of drugs in the drug toxicity data w/o duplicates:', len(drug_tox_ls))
print('')

# find the similar drug names in both drug lists

intersection = set(drug_tox_ls).intersection(drugbank_names)
drug_intersection_list = list(intersection)
print('# of drugs in both the toxicity data & the PathFX drugbank:', len(drug_intersection_list))
drug_diff_list = [x for x in drug_tox_ls if x not in drug_intersection_list]
print('# of drugs in the toxicity data but not in the PathFX drugbank:', len(drug_diff_list))
print('')

##################################################

print('##################################################')
print('')

'''
Map the side effects in the drug toxicity data to the PathFX datbase
'''

print('The side effect mapping process:')
print('')

# change the directory

os.chdir('/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/PathFX/rscs/')

# read the file including phenotypes and genes in PathFX

phen = pickle.load(open(('Pfx050120_all_assoc_to_nodes.pkl'),'rb'))
print('# of phenotype-genes pairs in the PathFX dictionary:', len(phen))
print('')
phen_ls = []
for name, value in phen.items():
  phen_ls.append(name.lower())
gen_ls = []
for name, value in phen.items():
  gen_ls.append(value)

# read the side effects in the drug toxicity data

tox_phen = df.columns.tolist()
for i in range(len(tox_phen)):
    tox_phen[i] = tox_phen[i].lower()

print('# of side effects in the drug toxicity data:', len(tox_phen))
print('')

# find the similar phenotype names in both lists

intersection = set(tox_phen).intersection(phen_ls)
phen_intersection_list = list(intersection)
print('# of side effects in the drug toxicity data & the phenotype-genes dictionary:', len(phen_intersection_list))
print('Side effect names in the drug toxicity data & the phenotype-genes dictionary:', phen_intersection_list)
phen_diff_list = [x for x in tox_phen if x not in phen_intersection_list]
print('# of side effects not in the the phenotype-genes dictionary:', len(phen_diff_list))
print('Side effect names not in the phenotype-genes dictionary:', phen_diff_list)

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

'''
Find the PathFX phenotypes relevant to labeled side effects
'''

print('The phenotype-matching process:')
print('')

def find_closet_match(string_name, list_name):

  '''
  Install this library every time you restart runtime on google colab:
  !pip install strsimpy

  # An ensemble function of three distinct techniques capable of detecting
  similar matches between two lists of words.
  - The first component of the string-matching search checks similar letters
  with the same positions in both strings and sees how many similar letters
  they have. Eventually, the algorithm takes the one with the highest number
  of similar letters and returns one match. This algorithm can be changed to
  return more similar matches too.
  - The second method uses a Python library called Jaro-Winkler that computes
  the similarity between 2 strings and the returned value lies in the
  interval of 0.0 and 1.0. When we have no similarity, the score will be 0.0.
  - The last component applies a Python module, difflib.get_close_matches,
  by giving back a list of good enough matches.
  '''

  import numpy as np
  import difflib
  from strsimpy.jaro_winkler import JaroWinkler

  !pip install strsimpy

  scores = {}
  for ii in list_name:
    cnt = 0
    if len(string_name)<=len(ii):
      str1, str2 = string_name, ii
    else:
      str1, str2 = ii, string_name
    for jj in range(len(str1)):
      cnt += 1 if str1[jj]==str2[jj] else 0
    scores[ii] = cnt
  scores_values = np.array(list(scores.values()))
  closest_match_idx = np.argsort(scores_values, axis=0, kind='quicksort')[-1]
  match1 = np.array(list(scores.keys()))[closest_match_idx]

  sims = []
  jw = JaroWinkler()
  for word in list_name:
    sims.append(jw.similarity(string_name, word))
    match2 = list_name[np.argmax(sims)]

  match3 = difflib.get_close_matches(string_name, list_name)

  similar_matches = [match1, match2]
  similar_matches.extend(match3)
  similar_matches_ls = set(similar_matches)
  print(similar_matches_ls)

  return

print('We developed a two step phonotype-matching approach, an ensemble function & validating the relevance by manual review of the literature!')
print('')
print('##################################################')

##################################################
