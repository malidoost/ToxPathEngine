
##################################################

'''
Define novel pathways
'''

##################################################

print('##################################################')
print('')

# function 1:

'''
A function to define new pathways using distinct TP genes
'''

def pfxg_input_tp (sideeffect):

  import os
  import pandas as pd
  import numpy as np
  import pickle

  # make the raw data with TP results.

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/Eval_SE2D_Res_Done/'))

  dfmet1 = pd.DataFrame()
  df_g = pd.read_csv('all_GenesResults.csv')
  df_genes = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'Distinct_TP_Genes'].tolist()

  g1 = df_genes[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  g4 = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    g4.append(gg)

  cui_dummy = ["cui_tp_" + ''.join(sideeffect) for i in range(len(g4))]
  Phene_name = ["phen_tp_" + ''.join(sideeffect) for i in range(len(g4))]

  res = list(zip(g4, cui_dummy, Phene_name))
  dfmet1 = dfmet1.append(res)

  dfmet1.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet1['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet1.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_DistinctTPGenes/', SE))
  dfmet1.to_csv('./' + '%s' %sideeffect.capitalize() + '_TPGenes.txt', index=False, sep='\t')

  ####################

  # randomization 1
  # gene knockout (15%)

  np.random.seed(7)

  genes_ls = g4
  randomsize = round(0.85 * len(genes_ls))
  gene_subset = np.random.choice(genes_ls, randomsize, replace=False)

  cui_dummy = ["cui_knockout_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_knockout_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet2 = pd.DataFrame()
  dfmet2 = dfmet2.append(res)

  dfmet2.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet2['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet2.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_DistinctTPGenes/', SE))
  dfmet2.to_csv('./' + '%s' %sideeffect.capitalize() + '_KnockoutTPGenes.txt', index=False, sep='\t')

  ####################

  # randomization 2
  # add noise (15%) from the shared genes list

  os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/Eval_SE2D_Res_Done/'))

  np.random.seed(7)

  dfmet3 = pd.DataFrame()

  df_g = pd.read_csv('all_GenesResults.csv')
  df_sh_genes = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'Shared_Genes'].tolist()
  g1 = df_sh_genes[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  g4 = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    g4.append(gg)

  genes_ls = g4
  randomsize = round(0.15 * len(genes_ls))
  gene_subset = np.random.choice(g4, randomsize, replace=False)

  cui_dummy = ["cui_noise_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_noise_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet3 = dfmet3.append(res)

  dfmet3.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet3['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet3.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_DistinctTPGenes/', SE))
  dfmet3.to_csv('./' + '%s' %sideeffect.capitalize() + '_NoiseSharedGenes.txt', index=False, sep='\t')

  ####################

  # combine all data sources

  dfmet = [dfmet1, dfmet2, dfmet3]
  df_data = pd.concat(dfmet)
  df_data.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_DistinctTPGenes/', SE))
  df_data.to_csv('./' + '%s' %sideeffect.capitalize() + '_TPRandoms.txt', index=False, sep='\t')

  ####################

  # TP & FP drugs

  os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/Eval_SE2D_Res_Done/'))

  dfmet11 = pd.DataFrame()
  dfmet12 = pd.DataFrame()

  df_g = pd.read_csv('all_GenesResults.csv')

  TP_drugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'TP_Drugs'].tolist()
  g1 = TP_drugs[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  TP_drugs_ls = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    TP_drugs_ls.append(gg)

  FP_drugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'FP_Drugs'].tolist()
  g1 = FP_drugs[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  FP_drugs_ls = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    FP_drugs_ls.append(gg)

  res11 = list(zip(TP_drugs_ls))
  dfmet11 = dfmet11.append(res11)

  res12 = list(zip(FP_drugs_ls))
  dfmet12 = dfmet12.append(res12)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_DistinctTPGenes/', SE))
  dfmet11.columns = ['TP_Drugs']
  dfmet11.to_csv('./' + '%s' %sideeffect.capitalize() + '_TPDrugs.txt', index=False, sep='\t')
  dfmet12.columns = ['FP_Drugs']
  dfmet12.to_csv('./' + '%s' %sideeffect.capitalize() + '_FPDrugs.txt', index=False, sep='\t')

  return df_data, dfmet11, dfmet12

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# function 2:

'''
A function to define new pathways using baseline & signature shared TP genes.
'''

def pfxg_input_dist_sig (sideeffect):

  import os
  import pandas as pd
  import numpy as np
  import pickle

  # distinct tp_drugs & fp_drugs

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/Eval_SE2D_Res_Done/'))

  df_g = pd.read_csv('all_GenesResults.csv')

  df_tpdrugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'TP_Drugs'].tolist()
  d1 = df_tpdrugs[0]
  d2 = d1.strip("[]")
  d3 = d2.split(',')
  d4 = []
  for i in range(len(d3)):
    d = d3[i].replace("'", "")
    dd = d.replace(" ", "")
    d4.append(dd)
  tp_drugs = d4

  df_fpdrugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'FP_Drugs'].tolist()
  d1 = df_fpdrugs[0]
  d2 = d1.strip("[]")
  d3 = d2.split(',')
  d4 = []
  for i in range(len(d3)):
    d = d3[i].replace("'", "")
    dd = d.replace(" ", "")
    d4.append(dd)
  fp_drugs = d4

  # signatures

  data_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/Data/'
  os.chdir(data_dir)
  df_drug_sig = pd.read_csv('PharmOmics_drug_signature_database.txt', sep='\t')
  df_hs = df_drug_sig.loc[(df_drug_sig['species'] == 'Homo sapiens')]
  df_dr = df_hs['drug'].tolist()
  df_dr_ls = [[x for x in word.lower().split()] for word in df_dr]
  df_drug_ls = [''.join(x) for x in df_dr_ls]

  # DrugBank

  dir1 = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/PathFX/rscs/'
  os.chdir(dir1)
  db2n = pickle.load(open(('drugbankid_to_name.pkl'),'rb'))
  drugbank_names = []
  for name, value in db2n.items():
    drugbank_names.append(value.lower())

  # map

  map_dr_ls = [x for x in df_drug_ls if x in drugbank_names]

  # map dataframe to analyze

  map_drs_ls = [i.capitalize() for i in map_dr_ls]
  df_data = df_hs.loc[df_hs['drug'].isin(map_drs_ls)]
  df_a = df_data[['drug', 'genes_top']]

  # similar drugs between tp & fp side effect drugs & sig drugs
  # and their corresponding genes in the sig data

  int_tp_dr = [x for x in tp_drugs if x in df_dr]
  df_b = df_a.loc[df_a['drug'].isin(int_tp_dr)]
  all_genes = df_b['genes_top'].tolist()
  cleaned_all_genes = [x for x in all_genes if pd.notnull(x)] # nan is not a string
  all_genes_ls1 = []
  for i in range(len(cleaned_all_genes)):
    agl = (cleaned_all_genes[i].split(','))
    all_genes_ls1.extend(agl)
  all_genes_list1 = list(set(all_genes_ls1))

  int_fp_dr = [x for x in fp_drugs if x in df_dr]
  df_b = df_a.loc[df_a['drug'].isin(int_fp_dr)]
  all_genes = df_b['genes_top'].tolist()
  cleaned_all_genes = [x for x in all_genes if pd.notnull(x)] # nan is not a string
  all_genes_ls2 = []
  for i in range(len(cleaned_all_genes)):
    agl = (cleaned_all_genes[i].split(','))
    all_genes_ls2.extend(agl)
  all_genes_list2 = list(set(all_genes_ls2))

  tp_genes_sig = all_genes_list1
  fp_genes_sig = all_genes_list2

  # distinct tp_genes & fp_genes

  df_genes = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'Sorted_TP_Genes'].tolist()
  g1 = df_genes[0]
  # convert string to list of tuples
  g1_ls = eval(g1)
  # extract strings from tuples
  tp_genes_dist = [m[0] for m in g1_ls]

  df_genes = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'Sorted_FP_Genes'].tolist()
  g1 = df_genes[0]
  # convert string to list of tuples
  g1_ls = eval(g1)
  # extract strings from tuples
  fp_genes_dist = [m[0] for m in g1_ls]

  # shared genes

  genTP_DisSig = [x for x in tp_genes_dist if x in tp_genes_sig]
  genFP_DisSig = [x for x in fp_genes_dist if x in fp_genes_sig]

  genTP_DisSig_noFP_Dis = [x for x in genTP_DisSig if x not in fp_genes_dist]
  genTP_DisSig_noFP_Sig = [x for x in genTP_DisSig if x not in fp_genes_sig]
  genTP_DisSig_noFP_DisSig = [x for x in genTP_DisSig_noFP_Dis if x not in fp_genes_sig]

  ###

  tp_genes = genTP_DisSig_noFP_Dis

  dfmet1a = pd.DataFrame()
  cui_dummy = ["cui_sig_tp_noFP_Dis_" + ''.join(sideeffect) for i in range(len(tp_genes))]
  Phene_name = ["phen_sig_tp_noFP_Dis_" + ''.join(sideeffect) for i in range(len(tp_genes))]
  res = list(zip(tp_genes, cui_dummy, Phene_name))
  dfmet1a = dfmet1a.append(res)

  datas = [dfmet1a]
  dfmet1 = pd.concat(datas)
  dfmet1.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet1['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet1.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet1.to_csv('./' + '%s' %sideeffect.capitalize() + '_Sig1TPGenes.txt', index=False, sep='\t')

  #

  tp_genes = genTP_DisSig_noFP_Sig

  dfmet4a = pd.DataFrame()
  cui_dummy = ["cui_sig_tp_noFP_Sig_" + ''.join(sideeffect) for i in range(len(tp_genes))]
  Phene_name = ["phen_sig_tp_noFP_Sig_" + ''.join(sideeffect) for i in range(len(tp_genes))]
  res = list(zip(tp_genes, cui_dummy, Phene_name))
  dfmet4a = dfmet4a.append(res)

  datas = [dfmet4a]
  dfmet4 = pd.concat(datas)
  dfmet4.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet4['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet4.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet4.to_csv('./' + '%s' %sideeffect.capitalize() + '_Sig2TPGenes.txt', index=False, sep='\t')

  ####################

  # randomization 1
  # gene knockout (15%)

  np.random.seed(7)

  genes_ls = genTP_DisSig_noFP_Dis
  randomsize = round(0.85 * len(genes_ls))
  gene_subset = np.random.choice(genes_ls, randomsize, replace=False)

  cui_dummy = ["cui_sig_knockout_noFP_Dis_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_sig_knockout_noFP_Dis_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet2 = pd.DataFrame()
  dfmet2 = dfmet2.append(res)

  dfmet2.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet2['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet2.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet2.to_csv('./' + '%s' %sideeffect.capitalize() + '_KnockoutSig1TPGenes.txt', index=False, sep='\t')

  #

  np.random.seed(7)

  genes_ls = genTP_DisSig_noFP_Sig
  randomsize = round(0.85 * len(genes_ls))
  gene_subset = np.random.choice(genes_ls, randomsize, replace=False)

  cui_dummy = ["cui_sig_knockout_noFP_Sig_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_sig_knockout_noFP_Sig_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet5 = pd.DataFrame()
  dfmet5 = dfmet5.append(res)

  dfmet5.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet5['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet5.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet5.to_csv('./' + '%s' %sideeffect.capitalize() + '_KnockoutSig2TPGenes.txt', index=False, sep='\t')

  ####################

  # randomization 2
  # add noise (15%) from the shared genes list

  np.random.seed(7)

  genes_ls = genFP_DisSig
  randomsize = round(0.15 * len(genes_ls))
  gene_subset = np.random.choice(genes_ls, randomsize, replace=False)

  cui_dummy = ["cui_sig_noise_noFP_Dis_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_sig_noise_noFP_Dis_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  dfmet3 = pd.DataFrame()
  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet3 = dfmet3.append(res)

  dfmet3.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet3['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet3.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet3.to_csv('./' + '%s' %sideeffect.capitalize() + '_NoiseSig1SharedGenes.txt', index=False, sep='\t')

  #

  np.random.seed(7)

  genes_ls = genFP_DisSig
  randomsize = round(0.15 * len(genes_ls))
  gene_subset = np.random.choice(genes_ls, randomsize, replace=False)

  cui_dummy = ["cui_sig_noise_noFP_Sig_" + ''.join(sideeffect) for i in range(len(gene_subset))]
  Phene_name = ["phen_sig_noise_noFP_Sig_" + ''.join(sideeffect) for i in range(len(gene_subset))]

  dfmet6 = pd.DataFrame()
  res = list(zip(gene_subset, cui_dummy, Phene_name))
  dfmet6 = dfmet6.append(res)

  dfmet6.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  dfmet6['Gene Symbol'].replace('', np.nan, inplace=True)
  dfmet6.dropna(subset=['Gene Symbol'], inplace=True)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet6.to_csv('./' + '%s' %sideeffect.capitalize() + '_NoiseSig2SharedGenes.txt', index=False, sep='\t')

  ####################

  # combine all data sources

  dfmet = [dfmet1, dfmet2, dfmet3, dfmet4, dfmet5, dfmet6]
  df_data = pd.concat(dfmet)
  df_data.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']
  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  df_data.to_csv('./' + '%s' %sideeffect.capitalize() + '_SigTPRandoms.txt', index=False, sep='\t')

  ####################

  # TP & FP drugs

  os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/Eval_SE2D_Res_Done/'))

  dfmet11 = pd.DataFrame()
  dfmet12 = pd.DataFrame()

  df_g = pd.read_csv('all_GenesResults.csv')

  TP_drugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'TP_Drugs'].tolist()
  g1 = TP_drugs[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  TP_drugs_ls = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    TP_drugs_ls.append(gg)

  FP_drugs = df_g.loc[(df_g['Side-Effect'] == sideeffect), 'FP_Drugs'].tolist()
  g1 = FP_drugs[0]
  g2 = g1.strip("[]")
  g3 = g2.split(',')
  FP_drugs_ls = []
  for i in range(len(g3)):
    g = g3[i].replace("'", "")
    gg = g.replace(" ", "")
    FP_drugs_ls.append(gg)

  res11 = list(zip(TP_drugs_ls))
  dfmet11 = dfmet11.append(res11)

  res12 = list(zip(FP_drugs_ls))
  dfmet12 = dfmet12.append(res12)

  SE = str(sideeffect.capitalize())
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_Omics/', SE))
  dfmet11.columns = ['TP_Drugs']
  dfmet11.to_csv('./' + '%s' %sideeffect.capitalize() + '_TPDrugs.txt', index=False, sep='\t')
  dfmet12.columns = ['FP_Drugs']
  dfmet12.to_csv('./' + '%s' %sideeffect.capitalize() + '_FPDrugs.txt', index=False, sep='\t')
  return df_data, dfmet11, dfmet12

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# function 3:

'''
A function to find the phenotype-gene associations for mapped 121 phenotypes.
'''

def pfxg_input_phen ():

  import os
  import pickle
  import pandas as pd

  # make the raw data with all 121 phenotypes and their corresponding genes

  phene_ls_all = [# 4 phenotypes
                  'neuroleptic malignant syndrome', 'delirium', 'completed suicide',
                  'hepatic necrosis',

              #####
              # 101 phenotypes

              'gastric ulcer', 'peptic ulcer',

              'hyperlipidemia', 'hypertriglyceridemia', 'familial hypercholesterolemia 1',
              'familial hypercholesterolemia', 'hypercholesterolemia',

              'edema', 'brain edema',

              'tardive dyskinesia', 'oral dyskinesia', 'drug-induced tardive dyskinesia',

              'myocardial infarction', 'myocardial infarction 1',
              'myocardial infarction susceptibility to, 1 (finding)', 'old myocardial infarction',
              'myocardial ischemia', 'acute myocardial infarction',
              'myocardial failure',

              'thrombocytopenia', 'thrombocytopenia 5', 'thrombocytopenia 6',
              'autoimmune thrombocytopenia', 'macrothrombocytopenia',
              'idiopathic thrombocytopenia', 'thrombocythemia, essential',
              'thrombocytopenia due to platelet alloimmunization',

              'cerebral infarction', 'brain infarction',
              'cerebrovascular accident', 'left middle cerebral artery infarction',
              'embolic infarction, middle cerebral artery',
              'right middle cerebral artery infarction',
              'brain ischemia', 'acute cerebrovascular accidents',
              'middle cerebral artery occlusion', 'cerebral hemorrhage',

              'pancreatitis', 'acute pancreatitis', 'pancreatitis idiopathic',
              'carcinoma of pancreas', 'adenocarcinoma of pancreas',
              'pancreatitis, chronic', 'pancreatitis, alcoholic',

              'peripheral neuropathy', 'peripheral motor neuropathy',
              'peripheral axonal neuropathy',
              'sciatic neuropathy', 'neuropathy',

              'hemorrhage', 'skin hemorrhages', 'brain hemorrhages',

              'proteinuria', 'mild proteinuria',

              'hypertension', 'prehypertension',
              'hypertensive disease', 'idiopathic pulmonary arterial hypertension',
              'genetic hypertension', 'pulmonary hypertension',
              'essential hypertension', 'hypertension, renovascular',
              'idiopathic pulmonary hypertension', 'renal hypertension',
              'familial primary, pulmonary hypertension', 'ocular hypertension',
              'hypertension, portal',

              'pulmonary edema', 'fluid retention in lung',

              'deep vein thrombosis', 'thrombosis', 'thrombus', 'thrombosis of inferior vena cava',
              'thrombocytosis',
              'middle cerebral artery thrombosis', 'venous thrombosis',

              'sepsis', 'sepsis of the newborn',
              'septicemia',

              'cardiac arrest', 'sudden cardiac arrest', 'cardiac arrest in children',

              'myopathy', 'gne myopathy', 'myopathy with exercise intolerance, swedish type',
              'inclusion body myopathy, sporadic', 'myofibrillar myopathy',
              'myotubular myopathy', 'tubular aggregate myopathy', 'cardiomyopathy',

              'pneumonia', 'pleuropneumonia',
              'lobar pneumonia', 'streptococcal pneumonia',
              'pneumocystis jiroveci pneumonia',

              'stevens-johnson syndrome toxic epidermal necrolysis spectrum',
              'stevens-johnson syndrome',
              'drug-induced stevens johnson syndrome',
              'mycoplasma-induced stevens-johnson syndrome',

              'agranulocytosis', 'granulocytosis',

              #####
              # 21 phenotypes

              'lung cyst',

              'anemia',
              'anemia of chronic disease',

              'tachycardia', 'tachyarrhythmia',

              'seizures', 'tonic seizures',
              'tonic - clonic seizures', 'seizures, clonic', 'convulsive seizures',
              'generalized seizures',

              'anaphylaxis (non medication)', 'anaphylaxis',

              'prolonged qt interval',

              'sleep disorders', 'sleep wake disorders',
              'sleep disturbances', 'sleeplessness',

              'sleep apnea syndromes', 'sleep apnea central', 'sleep apnea obstructive']

  os.chdir(os.path.join(main_dir, 'PathFX/rscs/'))
  p2c = pickle.load(open('Pfx050120_all_phens_to_cuis.pkl','rb'))
  name_ls_all_r = [name.lower() for name, value in p2c.items() if name.lower() in phene_ls_all]
  cui_ls_all_r = [value for name, value in p2c.items() if name.lower() in phene_ls_all]

  name_ls_all = list(set(name_ls_all_r))
  cui_ls_all = list(set(cui_ls_all_r))
  # venous thrombosis : C3160733
  # venous thrombosis : C0042487
  cui_ls_all.remove('C0042487')
  cui_ls = cui_ls_all

  c2g = pickle.load(open('Pfx050120_merged_unique_cuis2genes.pkl','rb'))

  genes_ls = []
  for na in cui_ls:
    for name, value in c2g.items():
      if na==name:
        genes_ls.append(value)

  res_ls = []
  for i in range(len(genes_ls)):
    gls = list(genes_ls[i])
    for j in range(len(gls)):
      ll = [gls[j], cui_ls_all[i], name_ls_all[i]]
      res_ls.append(ll)
  dfmet = pd.DataFrame()
  dfmet = dfmet.append(res_ls)
  dfmet.columns = ['Gene Symbol', 'CUI Term', 'Phenotype Name']

  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Pathways_121Phenotypes/'))
  dfmet.to_csv('./121PhenotypesGenes.txt', index=False, sep='\t')

  return dfmet

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# include all functions

import os
import pandas as pd

main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'

se_ls = ['hypertension', 'pancreatitis', 'thrombocytopenia',
         'myocardial infarction']

appended_df_tp = pd.DataFrame()
appended_df_tpdr = pd.DataFrame()
appended_df_fpdr = pd.DataFrame()
appended_df_tpsig = pd.DataFrame()
appended_df_tpdrsig = pd.DataFrame()
appended_df_fpdrsig = pd.DataFrame()

for sideeffect in se_ls:
  df_tp, df_tpdr, df_fpdr = pfxg_input_tp (sideeffect)
  df_tpsig, df_tpdrsig, df_fpdrsig = pfxg_input_dist_sig (sideeffect)

  appended_df_tp=appended_df_tp.append(df_tp, ignore_index=True)
  appended_df_tpsig=appended_df_tpsig.append(df_tpsig, ignore_index=True)

df_phenegene = pfxg_input_phen ()

os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/'))

dfall = [appended_df_tp, appended_df_tpsig, df_phenegene]
df_data = pd.concat(dfall)
phengene = df_data.drop_duplicates()
phengene.to_csv('./pfxg_input_phengene.txt', index=False, sep='\t')

#tp_drug = appended_df_tpdr.drop_duplicates()
#tp_drug.columns = ['TP_Drugs']
#tp_drug.to_csv('./pfxg_tpdrug.txt', index=False, sep='\t')
#fp_drug = appended_df_fpdr.drop_duplicates()
#fp_drug.columns = ['FP_Drugs']
#fp_drug.to_csv('./pfxg_fpdrug.txt', index=False, sep='\t')

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# useful functions for post-processing

# p-values statistics

def pval_stat (pval_ls):

  from statistics import mean, median, stdev

  min_pval_ls = min(pval_ls)
  print('min p_value:', min_pval_ls)
  max_pval_ls = max(pval_ls)
  print('max p_value:', max_pval_ls)
  mean_pval_ls = mean(pval_ls)
  print('mean p_value:', mean_pval_ls)
  median_pval_ls = median(pval_ls)
  print('median p_value:', median_pval_ls)
  stdev_pval_ls = stdev(pval_ls)
  print('stdev p_value:', stdev_pval_ls)

  return

####################

# pvalues_tp & pvalues_fp

def pval_tp_fp (file_name):

  import os
  import pandas as pd

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/P-values/'))

  pval_ls_tp = []
  df1 = pd.read_csv('%s.csv' %file_name)
  bh1 = df1['Benjamini-Hochberg'].tolist()
  for j in range(0, df1.shape[0]):
    g1 = bh1[j]
    g2 = g1.strip("[]")
    g3 = g2.split(',')
    for i in range(len(g3)):
      g = g3[i].replace("'", "")
      gg = g.replace(" ", "")
      pval_ls_tp.append(gg)

  pval_ls_tp_f = []
  for i in pval_ls_tp:
    new_i = float(i)
    pval_ls_tp_f.append(new_i)

  print('TP-values statistics:')
  pval_stat(pval_ls_tp_f)

  return pval_ls_tp_f

####################

# p_value plots

def pval_plot (sideeffect, pval_ls_tp, pval_ls_fp):

  import os
  import matplotlib.pyplot as plt
  import numpy as np

  plt.figure(1)
  plt.figure(figsize=(10, 8))
  plt.hist([pval_ls_tp, pval_ls_fp], color=['b', 'y'], label=['TP P-values', 'FP P-values'])
  plt.legend(loc='upper right')
  plt.title('P-values for %s phenotypes' %sideeffect)
  plt.xlabel('P-values Range')
  plt.ylabel('Frequency')

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/P-values/'))
  plt.savefig('P-values_' + '%s' %sideeffect + '.png')

  return

####################

# count the number of pathways

def pval_inf (sideeffect):

  import os
  import pandas as pd

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/P-values'))

  df_tp = pd.read_csv('pvalues_tpdrugs_' + '%s' %sideeffect + '.csv')
  df_fp = pd.read_csv('pvalues_fpdrugs_' + '%s' %sideeffect + '.csv')

  df_tp_tp = df_tp.loc[(df_tp.Phenotype == 'phen_tp_' + '%s' %sideeffect)]
  print('')
  print(df_tp_tp.shape)
  print('df_tp_tp:\n', df_tp_tp)
  df_fp_tp = df_fp.loc[(df_fp.Phenotype == 'phen_tp_' + '%s' %sideeffect)]
  print('')
  print(df_fp_tp.shape)
  print('df_fp_tp:\n', df_fp_tp)
  df_tp_tp = df_tp.loc[(df_tp.Phenotype == 'phen_sig_tp_' + '%s' %sideeffect)]
  print('')
  print(df_tp_tp.shape)
  print('df_tp_tp_sig:\n', df_tp_tp)
  df_fp_tp = df_fp.loc[(df_fp.Phenotype == 'phen_sig_tp_' + '%s' %sideeffect)]
  print('')
  print(df_fp_tp.shape)
  print('df_fp_tp_sig:\n', df_fp_tp)

  df_tp_ko = df_tp.loc[(df_tp.Phenotype == 'phen_knockout_' + '%s' %sideeffect)]
  print('')
  print(df_tp_ko.shape)
  print('df_tp_knockout:\n', df_tp_ko)
  df_fp_ko = df_fp.loc[(df_fp.Phenotype == 'phen_knockout_' + '%s' %sideeffect)]
  print('')
  print(df_fp_ko.shape)
  print('df_fp_knockout:\n', df_fp_ko)
  df_tp_ko = df_tp.loc[(df_tp.Phenotype == 'phen_sig_knockout_' + '%s' %sideeffect)]
  print('')
  print(df_tp_ko.shape)
  print('df_tp_knockout_sig:\n', df_tp_ko)
  df_fp_ko = df_fp.loc[(df_fp.Phenotype == 'phen_sig_knockout_' + '%s' %sideeffect)]
  print('')
  print(df_fp_ko.shape)
  print('df_fp_knockout_sig:\n', df_fp_ko)

  df_tp_n = df_tp.loc[(df_tp.Phenotype == 'phen_noise_' + '%s' %sideeffect)]
  print('')
  print(df_tp_n.shape)
  print('df_tp_noise:\n', df_tp_n)
  df_fp_n = df_fp.loc[(df_fp.Phenotype == 'phen_noise_' + '%s' %sideeffect)]
  print('')
  print(df_fp_n.shape)
  print('df_fp_noise:\n', df_fp_n)
  df_tp_n = df_tp.loc[(df_tp.Phenotype == 'phen_sig_noise_' + '%s' %sideeffect)]
  print('')
  print(df_tp_n.shape)
  print('df_tp_noise_sig:\n', df_tp_n)
  df_fp_n = df_fp.loc[(df_fp.Phenotype == 'phen_sig_noise_' + '%s' %sideeffect)]
  print('')
  print(df_fp_n.shape)
  print('df_fp_noise_sig:\n', df_fp_n)

  return

print('')
print('##################################################')

##################################################
