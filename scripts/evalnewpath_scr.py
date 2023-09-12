
##################################################

print('##################################################')
print('')

# Genes associated with TP & FP drugs (baseline & signature)

'''
A function to find genes associated with TP & FP drugs (baseline & signature)
'''

def genes_dist_sig (sideeffect):

  import os
  import pandas as pd
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

  # mapped dataframe to analyze

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

  return tp_genes_sig, fp_genes_sig, tp_genes_dist, fp_genes_dist, genTP_DisSig, genFP_DisSig, genTP_DisSig_noFP_Dis, genTP_DisSig_noFP_Sig, genTP_DisSig_noFP_DisSig

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# Sensitivity & Specificity of drugs with new defined pathways (per side effect)

'''
- name the side effects
- change the directory
'''

se_ls = ['hypertension', 'pancreatitis', 'thrombocytopenia',
         'myocardial infarction']

rows = []
for sideeffect in se_ls:
  print('%s:' %sideeffect)
  print('')
  print('truth:')

  # truth list

  import os
  import pandas as pd

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'

  # all 890 drugs

  path1 = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/Files/Results2/Eval_SE2D/SE2D_analysis_files/'
  files_and_folders = os.listdir(path1)
  folders = [f for f in files_and_folders if os.path.isdir(os.path.join(path1, f))]

  num_folders = len(folders)
  #print(f'There are {num_folders} folders in {path1}:')

  folder_ls = []
  for folder in folders:
    folder_ls.append(folder)

  # drugs for each side effect using the drug toxicity data (truth)

  os.chdir(main_dir)
  df = pd.read_csv('Drugs_labeled_for_AEs.txt', sep='\t', low_memory=False)

  # map all phenotypes with no direct matches to the side effects
  orig_se_ls = ['seizures', 'lung cyst', 'anemia', 'tachycardia', 'sleep disorders',
                'sleep apnea syndromes', 'anaphylaxis', 'Prolonged QT interval']
  map_to_orig_se_dict = {
    'lung cyst': 'Interstitial lung disease',
    'anemia': 'Hemolytic anemia',
    'tachycardia': 'Ventricular tachycardia',
    'seizures': 'Generalized tonic-clonic seizure',
    'anaphylaxis': 'Anaphylactic reaction',
    'Prolonged QT interval': 'QT Prolongation',
    'sleep disorders': 'Sleep disorder',
    'sleep apnea syndromes': 'Sleep apnea syndrome'
    }

  if sideeffect.capitalize() in df:
    dfwr = df[sideeffect.capitalize()]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  elif sideeffect.title() in df:
    dfwr = df[sideeffect.title()]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  elif sideeffect in orig_se_ls:
    orig_se = map_to_orig_se_dict[sideeffect]
    dfwr = df[orig_se]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  else:
    dfse_ls = []
    print('sideeffect not in the drug toxicity dataset')

  all_drugs_ls = [string.capitalize() for string in folder_ls]
  truth_ls_ = [x for x in dfse_ls if x in all_drugs_ls]
  truth_ls = list(set(truth_ls_))
  print('%s truth drugs list to model:' %sideeffect, len(truth_ls))

  print('')
  print('prediction:')

  # prediction list

  # change the directory
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Evals/Eval4/'))

  df = pd.read_csv('pathways_' + '%s' %sideeffect + '.csv')

  df_tp = df.loc[(df.Phenotype == 'phen_tp_' + '%s' %sideeffect)]
  print('df_tp:', df_tp.shape)

  dfs_tp = df.loc[(df.Phenotype == 'phen_sig_tp_' + '%s' %sideeffect)]
  print('df_tp_sig:', dfs_tp.shape)

  df1 = df.loc[(df.Phenotype != 'phen_knockout_' + '%s' %sideeffect)]
  df2 = df1.loc[(df1.Phenotype != 'phen_sig_knockout_' + '%s' %sideeffect)]
  df3 = df2.loc[(df2.Phenotype != 'phen_noise_' + '%s' %sideeffect)]
  df4 = df3.loc[(df3.Phenotype != 'phen_sig_noise_' + '%s' %sideeffect)]

  print('df without knockout & noise:', df4.shape)

  dr_ls = df['Drug_name'].tolist()
  dr_pred = list(set(dr_ls))
  print('%s prediction drugs list to model:' %sideeffect, len(dr_pred))
  print('')

  # evaluation

  '''
  results: sensitivity & specificity & precision
  '''

  print('results: confusion matrix')
  print('')

  pred_ls = dr_pred
  print('prediction list:', len(pred_ls))
  truth_ls = truth_ls
  print('truth list:', len(truth_ls))
  print('')

  TP_ls = [x for x in pred_ls if x in truth_ls]
  TP = len(TP_ls)
  print('TPs:', TP)
  print('TP drugs:', TP_ls)
  print('')

  FP_ls = [x for x in pred_ls if x not in truth_ls]
  FP = len(FP_ls)
  print('FPs:', FP)
  print('FP drugs:', FP_ls)
  print('')

  FN_ls = [x for x in truth_ls if x not in pred_ls]
  FN = len(FN_ls)
  print('FNs:', FN)
  print('FN drugs:', FN_ls)
  print('')

  ls_predtruth = list(set(pred_ls + truth_ls))

  all_drugs_ls = all_drugs_ls
  TN_ls = [x for x in all_drugs_ls if x not in ls_predtruth]
  TN = len(TN_ls)
  print('TNs:', TN)
  print('TN drugs:', TN_ls)
  print('')

  print('results: metrics')
  print('')

  conf_sensitivity = (TP / float(TP + FN + 0.00000001))
  conf_specificity = (TN / float(TN + FP + 0.00000001))
  conf_precision = (TP / float(TP + FP + 0.00000001))

  print(f'Sensitivity: {round(conf_sensitivity,2)}')
  print(f'Specificity: {round(conf_specificity,2)}')
  print(f'Precision: {round(conf_precision,2)}')
  print('')
  print('')
  print('')

  header = ('Side-Effect', 'TP', 'TN', 'FP', 'FN', 'Sensitivity', 'Specificity',
  'Precision', 'Distnct TP Pathways', 'Signature TP Pathways')
  met = [sideeffect, TP, TN, FP, FN, conf_sensitivity, conf_specificity,
  conf_precision, df_tp.shape, dfs_tp.shape]
  rows.append(met)
  dfmet = pd.DataFrame.from_records(rows, columns=header)
  dfmet.to_csv('./SideeffectMetrics_NewPathways.csv', index=False)

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# Sensitivity & Specificity of drugs with new defined pathways (per phenotype)

'''
- name the side effects
- change the directory: line 90
'''

se_ls = ['hypertension', 'pancreatitis', 'thrombocytopenia',
         'myocardial infarction']

for sideeffect in se_ls:
  print('%s:' %sideeffect)
  print('')
  print('truth:')

  # truth list

  import os
  import pandas as pd

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  #os.chdir(os.path.join(main_dir, 'Files/Results2/Eval_SE2D/SE2D_analysis_files/'))

  # all 890 drugs

  path1 = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/Files/Results2/Eval_SE2D/SE2D_analysis_files/'
  files_and_folders = os.listdir(path1)
  folders = [f for f in files_and_folders if os.path.isdir(os.path.join(path1, f))]

  num_folders = len(folders)
  #print(f'There are {num_folders} folders in {path1}:')

  folder_ls = []
  for folder in folders:
    folder_ls.append(folder)

  # drugs for each side effect using the drug toxicity data (truth)

  os.chdir(main_dir)
  df = pd.read_csv('Drugs_labeled_for_AEs.txt', sep='\t', low_memory=False)

  # map all phenotypes with no direct matches to the side effects
  orig_se_ls = ['seizures', 'lung cyst', 'anemia', 'tachycardia', 'sleep disorders',
                'sleep apnea syndromes', 'anaphylaxis', 'Prolonged QT interval']
  map_to_orig_se_dict = {
    'lung cyst': 'Interstitial lung disease',
    'anemia': 'Hemolytic anemia',
    'tachycardia': 'Ventricular tachycardia',
    'seizures': 'Generalized tonic-clonic seizure',
    'anaphylaxis': 'Anaphylactic reaction',
    'Prolonged QT interval': 'QT Prolongation',
    'sleep disorders': 'Sleep disorder',
    'sleep apnea syndromes': 'Sleep apnea syndrome'
    }

  if sideeffect.capitalize() in df:
    dfwr = df[sideeffect.capitalize()]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  elif sideeffect.title() in df:
    dfwr = df[sideeffect.title()]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  elif sideeffect in orig_se_ls:
    orig_se = map_to_orig_se_dict[sideeffect]
    dfwr = df[orig_se]
    dfse = dfwr.dropna()
    dfse_ls = dfse.tolist()
    print('drugs associated with %s on the drug toxicity data:' %sideeffect, len(dfse_ls))
  else:
    dfse_ls = []
    print('sideeffect not in the drug toxicity dataset')

  all_drugs_ls = [string.capitalize() for string in folder_ls]
  truth_ls_ = [x for x in dfse_ls if x in all_drugs_ls]
  truth_ls = list(set(truth_ls_))
  print('%s truth drugs list to model:' %sideeffect, len(truth_ls))

  print('')
  print('prediction:')
  print('')
  print('')
  print('')

  # prediction list

  # change the directory
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Evals/EvalAll1/'))

  # phenotype lists
  SE_name = '%s' %sideeffect

  if 'hypertension' == SE_name:
    SE_ls = ['hypertension', 'prehypertension', 'hypertensive disease',
             'idiopathic pulmonary arterial hypertension',
             'genetic hypertension', 'pulmonary hypertension',
             'essential hypertension', 'hypertension, renovascular',
             'idiopathic pulmonary hypertension', 'renal hypertension',
             'ocular hypertension', 'hypertension, portal',
             #
             'phen_tp_hypertension',
             #
             'phen_sig_tp_noFP_Dis_hypertension',
             'phen_sig_tp_noFP_Sig_hypertension']

  if 'pancreatitis' == SE_name:
    SE_ls = ['pancreatitis', 'acute pancreatitis', 'pancreatitis idiopathic',
             'carcinoma of pancreas', 'adenocarcinoma of pancreas',
             'pancreatitis, chronic', 'pancreatitis, alcoholic',
             #
             'phen_tp_pancreatitis',
             #
             'phen_sig_tp_noFP_Dis_pancreatitis',
             'phen_sig_tp_noFP_Sig_pancreatitis']

  if 'thrombocytopenia' == SE_name:
    SE_ls = ['thrombocytopenia', 'thrombocytopenia 5', 'thrombocytopenia 6',
             'autoimmune thrombocytopenia', 'macrothrombocytopenia',
             'idiopathic thrombocytopenia', 'thrombocythemia, essential',
             'thrombocytopenia due to platelet alloimmunization',
             #
             'phen_tp_thrombocytopenia',
             #
             'phen_sig_tp_noFP_Dis_thrombocytopenia',
             'phen_sig_tp_noFP_Sig_thrombocytopenia']

  if 'myocardial infarction' == SE_name:
    SE_ls = ['myocardial infarction', 'myocardial infarction 1',
             'myocardial infarction susceptibility to, 1 (finding)',
             'old myocardial infarction', 'myocardial ischemia',
             'acute myocardial infarction', 'myocardial failure',
             #
             'phen_tp_myocardial infarction',
             #
             'phen_sig_tp_noFP_Dis_myocardial infarction',
             'phen_sig_tp_noFP_Sig_myocardial infarction']

  # reading the SE file
  df = pd.read_csv('pathways_' + '%s' %sideeffect + '.csv')

  # analyzing the phenotypes
  rows = []
  for phenotype in SE_ls:
    df_phene = df.loc[(df.Phenotype == phenotype)]
    print(df_phene)
    dr_ls = df_phene['Drug_name'].tolist()
    dr_pred = list(set(dr_ls))
    print('%s prediction drugs list to model:' %phenotype, len(dr_pred))

    if len(dr_pred) != 0:
      gene_ls = df_phene['Genes'].tolist()
      fgenels = []
      # Iterate over each element in the input list
      for item in gene_ls:
        # Remove the outer brackets and single quotes
        cleaned_item = item.strip("[]").replace("'", "")
        # Split the string by commas and extend the final list
        fgenels.extend(cleaned_item.split(","))
        gene_pred = list(set(fgenels))
      print('%s prediction genes list:' %phenotype, len(gene_pred))
    else:
      gene_pred = []
      print('%s prediction genes list:' %phenotype, len(gene_pred))


    df_input = pd.read_csv('pfxg_input_phengene_evalall_4SEs1.txt', sep='\t', low_memory=False)
    gene_counts = df_input['Phenotype Name'].value_counts()
    for side_effect, count in gene_counts.items():
      if side_effect == phenotype:
        inputgene_count = count
    print('%s input genes list:' %phenotype, inputgene_count)
    print('')

    # evaluation

    '''
    results: sensitivity & specificity & precision
    '''

    print('results: confusion matrix')
    print('')

    pred_ls = dr_pred
    print('prediction list:', len(pred_ls))
    truth_ls = truth_ls
    print('truth list:', len(truth_ls))
    print('')

    TP_ls = [x for x in pred_ls if x in truth_ls]
    TP = len(TP_ls)
    print('TPs:', TP)
    print('TP drugs:', TP_ls)
    print('')

    FP_ls = [x for x in pred_ls if x not in truth_ls]
    FP = len(FP_ls)
    print('FPs:', FP)
    print('FP drugs:', FP_ls)
    print('')

    FN_ls = [x for x in truth_ls if x not in pred_ls]
    FN = len(FN_ls)
    print('FNs:', FN)
    print('FN drugs:', FN_ls)
    print('')

    ls_predtruth = list(set(pred_ls + truth_ls))

    all_drugs_ls = all_drugs_ls
    TN_ls = [x for x in all_drugs_ls if x not in ls_predtruth]
    TN = len(TN_ls)
    print('TNs:', TN)
    print('TN drugs:', TN_ls)
    print('')

    print('results: metrics')
    print('')

    conf_sensitivity = (TP / float(TP + FN + 0.00000001))
    conf_specificity = (TN / float(TN + FP + 0.00000001))
    conf_precision = (TP / float(TP + FP + 0.00000001))

    print(f'Sensitivity: {round(conf_sensitivity,2)}')
    print(f'Specificity: {round(conf_specificity,2)}')
    print(f'Precision: {round(conf_precision,2)}')
    print('')
    print('')
    print('')

    header = ('Side-Effect', 'Phenotype', 'TP', 'TN', 'FP', 'FN',
              'Sensitivity', 'Specificity', 'Precision',
              '# of Predicted Associated Genes', '# of Input Associated Genes' )
    met = [sideeffect, phenotype, TP, TN, FP, FN,
           conf_sensitivity, conf_specificity, conf_precision, len(gene_pred),
           inputgene_count]
    rows.append(met)
    dfmet = pd.DataFrame.from_records(rows, columns=header)
    sorted_dfmet = dfmet.sort_values(by='Specificity', ascending=False)
  SE = str(sideeffect.capitalize())
  sorted_dfmet.to_csv('./' + '%s' %SE + '_PhenotypesMetrics_NewPathways.csv', index=False)

print('')
print('##################################################')

##################################################

print('##################################################')
print('')

# Evaluation plots for 4 SEs

import os
import pandas as pd
import matplotlib.pyplot as plt

se_ls = ['Hypertension', 'Thrombocytopenia', 'Myocardial Infarction', 'Pancreatitis']

for sideeffect in se_ls:

  main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'
  os.chdir(os.path.join(main_dir, 'Files/PathFX_Gen_Res/Evals/Metrics2/'))

  df_hyp = pd.read_csv('MetricsPlot_' + '%s' %sideeffect + '.csv')
  x = df_hyp['Specificity'].tolist()
  y = df_hyp['Sensitivity'].tolist()

  markers = ['o', 's', 'D', '^', 'P', '*', 'X']
  colors = ['r', 'g', 'b', 'm', 'y', 'k', 'c']
  labels = df_hyp['Pathways'].tolist()

  fig, ax = plt.subplots()

  for i in range(len(x)):
    ax.scatter(x[i], y[i], marker=markers[i], color=colors[i], label=f'{labels[i]}')

  ax.set_xlabel('Specificity')
  ax.set_ylabel('Sensitivity')
  ax.set_xlim([-0.1, 1.1])
  ax.set_ylim([-0.1, 1.1])
  ax.legend(title=None, loc='upper left')
  if sideeffect == 'Pancreatitis':
    organ_sys = 'Gastrointestinal'
  else:
    organ_sys = 'Cardiovascular'
  ax.set_title('Evaluation metrics for ' + '%s' %sideeffect + ' (' + '%s' %organ_sys + ')')

  plt.savefig('Evaluation_Plot_' + '%s' %sideeffect + '.png', dpi=300)

  plt.show()

print('')
print('##################################################')

##################################################
