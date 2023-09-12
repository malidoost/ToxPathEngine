
##################################################

print('##################################################')
print('')

'''
Map phenotypes cuis to side effects cuis to present the evaluations per side effect.
'''

# map phenotypes cuis to side effects cuis

map_to_orig_cui_dict = {

  'C0040034': 'C0040034',
  'C4310789': 'C0040034',
  'C4015537': 'C0040034',
  'C0040028': 'C0040034',
  'C0242584': 'C0040034',
  'C0920163': 'C0040034',
  'C2751260': 'C0040034',
  'C0272286': 'C0040034',
  'C0038358': 'C0038358',
  'C0030920': 'C0038358',
  'C0013604': 'C0013604',
  'C1527311': 'C0013604',
  'C0686347': 'C0686347',
  'C3714760': 'C0686347',
  'C0454606': 'C0686347',
  'C0020473': 'C0020473',
  'C0020445': 'C0020473',
  'C0020557': 'C0020473',
  'C0020443': 'C0020473',
  'C0745103': 'C0020473',
  'C0027051': 'C0027051',
  'C0155668': 'C0027051',
  'C1832662': 'C0027051',
  'C0151744': 'C0027051',
  'C0155626': 'C0027051',
  'C1959583': 'C0027051',
  'C0032285': 'C0032285',
  'C0032241': 'C0032285',
  'C0032300': 'C0032285',
  'C1535939': 'C0032285',
  'C0155862': 'C0032285',
  'C0001824': 'C0001824',
  'C1282609': 'C0001824',
  'C0038325': 'C0038325',
  'C3658302': 'C0038325',
  'C1274933': 'C0038325',
  'C3658301': 'C0038325',
  'C0034063': 'C0034063',
  'C0848538': 'C0034063',
  'C0243026': 'C0243026',
  'C0456103': 'C0243026',
  'C0036690': 'C0243026',
  'C0018790': 'C0018790',
  'C3826614': 'C0018790',
  'C1720824': 'C0018790',
  'C0149871': 'C0149871',
  'C0040053': 'C0149871',
  'C0087086': 'C0149871',
  'C0836924': 'C0149871',
  'C2712843': 'C0149871',
  'C0042487': 'C0149871',
  'C0740376': 'C0149871',
  'C3278737': 'C3278737',
  'C1696708': 'C3278737',
  'C3203102': 'C3278737',
  'C0085580': 'C3278737',
  'C0020542': 'C3278737',
  'C0020545': 'C3278737',
  'C0152171': 'C3278737',
  'C0028840': 'C3278737',
  'C0020541': 'C3278737',
  'C0598428': 'C3278737',
  'C0020544': 'C3278737',
  'C0020538': 'C3278737',
  'C0026848': 'C0026848',
  'C1853926': 'C0026848',
  'C1850718': 'C0026848',
  'C2678065': 'C0026848',
  'C0175709': 'C0026848',
  'C0410207': 'C0026848',
  'C0878544': 'C0026848',
  'C0751713': 'C0026848',
  'C0033687': 'C0033687',
  'C4022832': 'C0033687',
  'C0019080': 'C0019080',
  'C0852361': 'C0019080',
  'C0031117': 'C0031117',
  'C1263857': 'C0031117',
  'C0235025': 'C0031117',
  'C0149940': 'C0031117',
  'C0442874': 'C0031117',
  'C0030305': 'C0030305',
  'C0747198': 'C0030305',
  'C0001339': 'C0030305',
  'C0376670': 'C0030305',
  'C0149521': 'C0030305',
  'C0235974': 'C0030305',
  'C0279176': 'C0030305',
  'C0007785': 'C0007785',
  'C0751955': 'C0007785',
  'C0751956': 'C0007785',
  'C0038454': 'C0007785',
  'C0751846': 'C0007785',
  'C0751847': 'C0007785',
  'C0751849': 'C0007785',
  'C0740391': 'C0007785',
  'C2937358': 'C0007785',
  'C0007786': 'C0007785',

  'C0027849': 'C0027849', #1
  'C0011206': 'C0011206', #1
  'C0852733': 'C0852733', #1
  'C0151798': 'C0151798', #1

  'C0039231': 'C0039231', #0
  'C0080203': 'C0039231', #0
  'C0151878': 'C0151878', #0
  'C0002792': 'C0002792', #0
  'C0850803': 'C0002792', #0
  'C0036572': 'C0036572', #0
  'C3809174': 'C0036572', #0
  'C0751494': 'C0036572', #0
  'C0234535': 'C0036572', #0
  'C0494475': 'C0036572', #0
  'C0234533': 'C0036572', #0
  'C0546483': 'C0546483', #0
  'C0002871': 'C0002871', #0
  'C0002873': 'C0002871', #0
  'C0037315': 'C0037315', #0
  'C0851578': 'C0851578', #0
  'C4042891': 'C0851578', #0
  'C0037317': 'C0851578', #0
  'C0917801': 'C0851578'  #0

}

##################################################

def cui_to_genes (cui_sideeffect):

    '''
    - This function reads all PathFX result tables and stores the important
    information in a table for the future analysis.
    - The input for this function is a CUI side effect.
    - The output will be a table of phenotypes, CUIs, drug names, and genes.
    '''

    import os
    import pandas as pd

    # dir1 is the folder including PathFX results for all 890 drugs
    dir1 = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/SE2D_analysis_files/'
    # dir2 is the folder to store the desired results
    dir2 = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/Eval_SE2D_Res/'

    file_ext = '_merged_neighborhood__assoc_table_.txt'

    genes_ls = []
    res_table = []

    for root, subdirectories, files in os.walk(dir1):

      for file_name in files:
        drugname = os.path.split(root)[1]

        if file_ext in file_name:
          file_name_to_read = file_name
          os.chdir(os.path.join(root))
          df_drug = pd.read_csv(file_name_to_read, sep='\t')
          cui_phene = df_drug['cui'].tolist()
          if cui_sideeffect in cui_phene:
            phen_row = df_drug.loc[df_drug['cui'] == cui_sideeffect]
            genes_ls = phen_row['genes'].tolist()
            phen = phen_row['phenotype'].tolist()
            err = 'NoErr'
          elif cui_sideeffect not in cui_phene:
            genes_ls = ['nogenes']
            phen = ['nophen']
            err = 'NoPhenErr'
          else:
            genes_ls = ['nogenes']
            phen = ['nophen']
            err = 'GeneralErr'

          res_table.append([phen, cui_sideeffect, drugname, err, genes_ls])

    header = ('Phenotype', 'CUI', 'Drug_name', 'Error', 'Genes')
    dfres = pd.DataFrame.from_records(res_table, columns=header)
    os.chdir(dir2)
    dfres.to_csv('./phene_gene_table.csv', index=False)

    return

##################################################

def eval_se2d (sideeffect):

    '''
    - This function outputs the evaluations per side effect.
    - The input for this function is a side effect.
    - The outputs will be the associated genes to the side effect,
    the true pathways, drugs, and the evaluation metrics.
    '''

    print('%s analysis results:' %sideeffect)
    print('')
    print('#####')
    print('')

    import os
    import pickle
    import pandas as pd
    from itertools import chain

    main_dir = '/content/gdrive/MyDrive/PhD_Lab/Project_Drug_Toxicity_Network_Predictions/'

    # convert the side effect to the CUI term

    os.chdir(os.path.join(main_dir, 'PathFX/rscs/'))
    p2c = pickle.load(open('Pfx050120_all_phens_to_cuis.pkl','rb'))

    se_name_pr = ['anaphylaxis', 'Prolonged QT interval']
    if sideeffect in se_name_pr:
      cui_label = [value for name, value in p2c.items() if name==sideeffect]
    else:
      cui_label = [value for name, value in p2c.items() if name==sideeffect.capitalize()]
      cui_label_extra = [value for name, value in p2c.items() if name==sideeffect.title()]
      cui_label.extend(cui_label_extra)

    cui_ls = list(set(cui_label))

    # collect all phenotype matches of a side effect in a list

    keys_cui_ls = []
    keys = [k for k, vs in map_to_orig_cui_dict.items() if ''.join(cui_ls) in vs]
    keys_cui_ls.append(keys)
    #for i in range(len(keys_cui_ls[0])):
    #  print(keys_cui_ls[0][i])

    # calculate the associated genes to the phenotypes

    rows1 = []
    genes_all_ls = []
    pred_drugs_ls = []
    for i in range(len(keys_cui_ls[0])):
      cui_term = keys_cui_ls[0][i]
      cui_to_genes (cui_term)

      # read the 'genes associated to the side effects' file

      os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))
      phg = pd.read_csv('phene_gene_table.csv')

      phen = [k for k, v in p2c.items() if v==cui_term]
      print('phenotype:', phen)

      phg.to_csv('./phene_gene_table_' + ''.join(phen) + '.csv', index=False)

      gene_ls1 = phg.loc[phg.CUI == cui_term, 'Genes'].tolist()

      # convert the list of strings to the list of lists
      # convert the list of lists to the list of strings
      # remove 'nogenes' strings
      # split a list by ','
      # convert the list of lists to the list of strings
      # count and sort the genes

      gene_ls2 = []
      for item in gene_ls1:
        gene2 = [x.strip() for x in eval(item)]
        gene_ls2.append(gene2)

      gene_ls3 = []
      for item in gene_ls2:
        gene3 = ''.join(item)
        gene_ls3.append(gene3)

      gene_ls4 = gene_ls3

      gene_ls5 = []
      for item in gene_ls4:
        gene5 = item.split(',')
        gene_ls5.append(gene5)

      genes_ls = list(chain(*gene_ls5))

      genes_ls_ = []
      rem_ls = ['nogenes']
      [genes_ls_.append(x) for x in genes_ls if x not in rem_ls]
      print('number of genes associated:', len(genes_ls_))

      ge_count_dic = {}
      for ge in genes_ls_:
        ge_count = genes_ls_.count(ge)
        ge_count_dic[ge] = ge_count
      sorted_ge_count_dic = sorted(ge_count_dic.items(), key=lambda x: x[1], reverse=True)
      print('sorted genes count:\n', sorted_ge_count_dic)
      print('')

      # drugs for each side effect using the PathFX results (prediction)

      drug_ls = phg.loc[(phg['CUI'] == cui_term) & (phg['Error'] == 'NoErr'), 'Drug_name'].tolist()
      print('# of drugs associated with %s in the PathFX prediction:' %phen, len(drug_ls))
      print('drugs associated with %s in the PathFX prediction:' %phen, (drug_ls))
      print('')

      # phenotypes analysis results

      header1 = ('Phenotype', 'CUI', '# of Associated Genes', 'Sorted Associated Genes',
                '# of Associated Drugs', 'Associated Drugs')
      met1 = [phen, cui_term, len(genes_ls_), sorted_ge_count_dic, len(drug_ls),
             drug_ls]
      rows1.append(met1)
      dfmet1 = pd.DataFrame.from_records(rows1, columns=header1)
      os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))
      dfmet1.to_csv('./PhenotypesAnalysisResults.csv', index=False)

      # collect genes and drugs information from all phenotypes

      genes_all_ls.extend(genes_ls_)
      pred_drugs_ls.extend(drug_ls)

    print('#####')
    print('')

    print('prediction:')
    print('')

    genes_all_ls_ = []
    rem_ls = ['nogenes']
    [genes_all_ls_.append(x) for x in genes_all_ls if x not in rem_ls]
    genes_all_lss = list(set(genes_all_ls_))
    print('# of genes associated with %s in the PathFX prediction:' %sideeffect, len(genes_all_ls_))
    print('all genes associated with %s in the PathFX prediction:' %sideeffect, genes_all_ls_)

    ge_count_dic_all = {}
    for ge in genes_all_ls_:
      ge_count = genes_all_ls_.count(ge)
      ge_count_dic_all[ge] = ge_count
    sorted_ge_all_count_dic = sorted(ge_count_dic_all.items(), key=lambda x: x[1], reverse=True)
    print('sorted genes count:\n', sorted_ge_all_count_dic)
    print('')

    print('# of drugs associated with %s in the PathFX prediction:' %sideeffect, len(pred_drugs_ls))
    print('all drugs associated with %s in the PathFX prediction:' %sideeffect, pred_drugs_ls)

    # drugs for each side effect using the drug toxicity data (truth)

    print('')
    print('#####')
    print('')
    print('truth:')
    os.chdir(main_dir)

    df = pd.read_csv('Drugs_labeled_for_AEs.txt', sep='\t', low_memory=False)

    # map all phenotypes with no direct matches to the side effects

    orig_se_ls = ['seizures', 'lung cyst', 'anemia', 'tachycardia', 'sleep disorders',
    'sleep apnea syndromes', 'anaphylaxis', 'Prolonged QT interval']

    map_to_orig_se_dict = {

      'lung cyst': 'Interstitial lung disease',
      'anemia': 'Hemolytic anemia',
      #'anemia of chronic disease': 'Hemolytic anemia',
      'tachycardia': 'Ventricular tachycardia',
      #'tachyarrhythmia': 'Ventricular tachycardia',
      'seizures': 'Generalized tonic-clonic seizure',
      #'tonic seizures': 'Generalized tonic-clonic seizure',
      #'tonic - clonic seizures': 'Generalized tonic-clonic seizure',
      #'seizures, clonic': 'Generalized tonic-clonic seizure',
      #'convulsive seizures': 'Generalized tonic-clonic seizure',
      #'generalized seizures': 'Generalized tonic-clonic seizure',
      'anaphylaxis': 'Anaphylactic reaction',
      #'anaphylaxis (non medication)': 'Anaphylactic reaction',
      'Prolonged QT interval': 'QT Prolongation',
      'sleep disorders': 'Sleep disorder',
      #'sleep wake disorders': 'Sleep disorder',
      #'sleep disturbances': 'Sleep disorder',
      #'sleeplessness': 'Sleep disorder',
      'sleep apnea syndromes': 'Sleep apnea syndrome'
      #'sleep apnea central': 'Sleep apnea syndrome',
      #'sleep apnea obstructive': 'Sleep apnea syndrome'

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

    print('')
    print('#####')
    print('')

    # all possible drugs for the analysis (890 drugs)

    dir1 = '/u/home/m/malidoos/Desktop/Drug_Toxicity_Network_Predictions/SE2D_analysis_files/'

    all_files = os.listdir(dir1)
    print("drugs in the analysis results folder: ", len(all_files))
    all_drugs_ls = [drug_name.capitalize() for drug_name in all_files]

    # remove Db##### and add DB#####

    db49_ls1 = ['DB00030', 'DB00046', 'DB00855', 'DB00008', 'DB03793', 'DB03619', 'DB00047',
                   'DB09539', 'DB00974', 'DB00591', 'DB00784', 'DB00513', 'DB06695', 'DB13151',
                   'DB00056', 'DB05773', 'DB01309', 'DB00158', 'DB00060', 'DB04398',
                   'DB03166', 'DB05990', 'DB00233', 'DB06770', 'DB00548', 'DB00313', 'DB01024',
                   'DB00302', 'DB01294', 'DB09269', 'DB00068', 'DB01306', 'DB00936', 'DB03017',
                   'DB09395', 'DB00033', 'DB00069', 'DB00883', 'DB00939', 'DB00078', 'DB09422',
                   'DB00992', 'DB11256', 'DB01307', 'DB00630', 'DB08904', 'DB02325', 'DB01272']
    db49_ls2 = ['Db00030', 'Db00046', 'Db00855', 'Db00008', 'Db03793', 'Db03619', 'Db00047',
                   'Db09539', 'Db00974', 'Db00591', 'Db00784', 'Db00513', 'Db06695', 'Db13151',
                   'Db00056', 'Db05773', 'Db01309', 'Db00158', 'Db00060', 'Db04398',
                   'Db03166', 'Db05990', 'Db00233', 'Db06770', 'Db00548', 'Db00313', 'Db01024',
                   'Db00302', 'Db01294', 'Db09269', 'Db00068', 'Db01306', 'Db00936', 'Db03017',
                   'Db09395', 'Db00033', 'Db00069', 'Db00883', 'Db00939', 'Db00078', 'Db09422',
                   'Db00992', 'Db11256', 'Db01307', 'Db00630', 'Db08904', 'Db02325', 'Db01272']

    os.chdir(os.path.join(main_dir, 'PathFX/rscs/'))

    db = pickle.load(open('pfxDB050620_dbid2name.pkl','rb'))
    drug_name = [name for value, name in db.items() if value in db49_ls1]
    for item in db49_ls2:
      if item in all_drugs_ls:
        all_drugs_ls.remove(item)
    all_drugs_ls.extend(drug_name)

    print('')
    print('#####')
    print('')

    # results

    '''
    results: TP & FP & FN & TN
    '''

    print('results: confusion matrix')
    print('')

    pred_ls_ = [x for x in pred_drugs_ls if x in all_drugs_ls]
    pred_ls = list(set(pred_ls_))
    print('prediction list:', len(pred_ls))
    truth_ls_ = [x for x in dfse_ls if x in all_drugs_ls]
    truth_ls = list(set(truth_ls_))
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

    TN_ls = [x for x in all_drugs_ls if x not in ls_predtruth]
    TN = len(TN_ls)
    print('TNs:', TN)
    print('TN drugs:', TN_ls)
    print('')

    # genes list including the TP predicted drug (pathway) & FP ones

    os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))

    tpdrugs_ls_ = []
    fpdrugs_ls_ = []
    genes_in_tpdrugs = []
    rows4 = []
    genes_in_fpdrugs = []
    rows5 = []
    for i in range(len(keys_cui_ls[0])):
      cui_term = keys_cui_ls[0][i]
      cui_to_genes (cui_term)
      phg = pd.read_csv('phene_gene_table.csv')
      if TP>0:
        for drugname in TP_ls:
          genes_in_tpdrug = phg.loc[(phg['CUI'] == cui_term) &
          (phg['Error'] == 'NoErr') &
          (phg['Drug_name'] == drugname), 'Genes'].tolist()
          tpdrugs_ls_.append(drugname)
          genes_in_tpdrugs.append(genes_in_tpdrug)

          met4 = [drugname, genes_in_tpdrug]
          rows4.append(met4)

        print('TP drugs:', tpdrugs_ls_)
        print('genes list with TP pathways:', genes_in_tpdrugs)
        Err = 'NoErr'
      else:
        genes_in_tpdrugs = ['NoGenes']
        print('no genes because TP=0')
        Err = 'TPIs0'

      header4 = ('Drug', 'Genes')
      dfmet4 = pd.DataFrame.from_records(rows4, columns=header4)
      dfmet4.to_csv('./TPDrGen_' + ''.join(sideeffect) + '.csv', index=False)

      print('')

      if FP>0:
        for drugname in FP_ls:
          genes_in_fpdrug = phg.loc[(phg['CUI'] == cui_term) &
          (phg['Error'] == 'NoErr') &
          (phg['Drug_name'] == drugname), 'Genes'].tolist()
          fpdrugs_ls_.append(drugname)
          genes_in_fpdrugs.append(genes_in_fpdrug)

          met5 = [drugname, genes_in_fpdrug]
          rows5.append(met5)

        print('FP drugs:', fpdrugs_ls_)
        print('genes list with FP pathways:', genes_in_fpdrugs)
        Err = 'NoErr'
      else:
        genes_in_fpdrugs = ['NoGenes']
        print('no genes because FP=0')
        Err = 'FPIs0'

      header5 = ('Drug', 'Genes')
      dfmet5 = pd.DataFrame.from_records(rows5, columns=header5)
      dfmet5.to_csv('./FPDrGen_' + ''.join(sideeffect) + '.csv', index=False)

    tpdrugs_ls = list(set(tpdrugs_ls_))
    fpdrugs_ls = list(set(fpdrugs_ls_))
    print('')

    # find the shared genes between the tpdrugs_ls & fpdrugs_ls

    genes_in_tpdrugs = [x for x in genes_in_tpdrugs if x]
    tpgene = []
    for i in range(len(genes_in_tpdrugs)):
      gene_ls1 = genes_in_tpdrugs[i][0].strip("[]")
      gene_ls2 = gene_ls1.strip("''")
      gene_ls3 = (gene_ls2.split(','))
      tpgene.append(gene_ls3)
      tpgenes = sum(tpgene, [])
    tpgenes_ = list(set(tpgenes))
    tpgenes_len = len(tpgenes_)
    ge_count_dic_all2 = {}
    for ge in tpgenes:
      ge_count = tpgenes.count(ge)
      ge_count_dic_all2[ge] = ge_count
    sorted_tp_genes = sorted(ge_count_dic_all2.items(), key=lambda x: x[1], reverse=True)

    genes_in_fpdrugs = [x for x in genes_in_fpdrugs if x]
    fpgene = []
    for i in range(len(genes_in_fpdrugs)):
      gene_ls4 = genes_in_fpdrugs[i][0].strip("[]")
      gene_ls5 = gene_ls4.strip("''")
      gene_ls6 = (gene_ls5.split(','))
      fpgene.append(gene_ls6)
      fpgenes = sum(fpgene, [])
    fpgenes_ = list(set(fpgenes))
    fpgenes_len = len(fpgenes_)
    ge_count_dic_all3 = {}
    for ge in fpgenes:
      ge_count = fpgenes.count(ge)
      ge_count_dic_all3[ge] = ge_count
    sorted_fp_genes = sorted(ge_count_dic_all3.items(), key=lambda x: x[1], reverse=True)

    shared_genes = [x for x in tpgenes if x in fpgenes]
    shared_genes_ = list(set(shared_genes))
    shared_genes_len = len(shared_genes_)

    distinct_tp_genes = [x for x in tpgenes if x not in fpgenes]
    distinct_tp_genes_ = list(set(distinct_tp_genes))
    distinct_tp_genes_len = len(distinct_tp_genes_)

    distinct_fp_genes = [x for x in fpgenes if x not in tpgenes]
    distinct_fp_genes_ = list(set(distinct_fp_genes))
    distinct_fp_genes_len = len(distinct_fp_genes_)

    print('')

    '''
    results: accuracy & sensitivity & specificity & precision & f_1 score
    '''

    print('results: metrics')
    print('')

    conf_sensitivity = (TP / float(TP + FN + 0.00000001))
    conf_specificity = (TN / float(TN + FP + 0.00000001))
    conf_precision = (TP / float(TP + FP + 0.00000001))

    print(f'Sensitivity: {round(conf_sensitivity,2)}')
    print(f'Specificity: {round(conf_specificity,2)}')
    print(f'Precision: {round(conf_precision,2)}')

    '''
    table of results
    '''

    os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))

    rows2 = []
    header2 = ('Side-Effect', 'TP', 'TN', 'FP', 'FN', 'Sensitivity', 'Specificity',
    'Precision', 'Truth_ToxData_Drugs', 'Truth_Model_Drugs',
     'Pred_PathFX_Drugs', 'Pred_Model_Drugs',
    'PathFX_Genes_all', 'PathFX_Genes', ' # of Genes_with_True_Pathways',
    '# of Genes_with_False_Pathways',
    '# of Shared_Genes_TP_FP_Pathways', '#_Distinct_TP_Genes', '#_Distinct_FP_Genes')
    met2 = [sideeffect, TP, TN, FP, FN, conf_sensitivity, conf_specificity,
    conf_precision, len(dfse_ls), len(truth_ls), len(pred_drugs_ls),
    len(pred_ls), len(genes_all_ls_), len(genes_all_lss), tpgenes_len,
    fpgenes_len, shared_genes_len, distinct_tp_genes_len, distinct_fp_genes_len]
    rows2.append(met2)
    dfmet2 = pd.DataFrame.from_records(rows2, columns=header2)
    dfmet2.to_csv('./SideeffectMetrics.csv', index=False)

    rows3 = []
    header3 = ('Side-Effect', '#_TP_Drugs', 'TP_Drugs', '#_TP_Genes', 'Sorted_TP_Genes',
    '#_FP_Drugs', 'FP_Drugs', '#_FP_Genes', 'Sorted_FP_Genes', '#_Shared_Genes', 'Shared_Genes',
    '#_Distinct_TP_Genes', 'Distinct_TP_Genes', '#_Distinct_FP_Genes', 'Distinct_FP_Genes')
    met3 = [sideeffect, len(tpdrugs_ls), tpdrugs_ls, tpgenes_len, sorted_tp_genes,
    len(fpdrugs_ls), fpdrugs_ls, fpgenes_len, sorted_fp_genes, shared_genes_len, shared_genes_,
    distinct_tp_genes_len, distinct_tp_genes_, distinct_fp_genes_len, distinct_fp_genes_]
    rows3.append(met3)
    dfmet3 = pd.DataFrame.from_records(rows3, columns=header3)
    dfmet3.to_csv('./GenesResults.csv', index=False)

    return

##################################################

# run the evaluation function for different side effects
# combine all results into one csv file

import os
from csv import reader
import pandas as pd
import sys
import csv

csv.field_size_limit(sys.maxsize)

main_dir = '/u/home/m/malidoos/Desktop/Drug_Toxicity_Network_Predictions/'
os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))

se_name = ['edema', 'gastric ulcer', 'neuroleptic malignant syndrome', 'delirium',
'hyperlipidemia', 'completed suicide', 'hepatic necrosis', 'tardive dyskinesia',
'proteinuria', 'hemorrhage', 'myocardial infarction', 'hypertension',
'deep vein thrombosis', 'sepsis', 'cardiac arrest', 'thrombocytopenia',
'agranulocytosis', 'stevens-johnson syndrome', 'cerebral infarction', 'pancreatitis',
'peripheral neuropathy', 'pulmonary edema', 'myopathy', 'pneumonia',
'seizures', 'lung cyst', 'tachycardia', 'sleep disorders',
'anemia', 'sleep apnea syndromes', 'anaphylaxis', 'Prolonged QT interval']

labels1 = ('Phenotype', 'CUI', '# of Associated Genes', 'Sorted Associated Genes',
              '# of Associated Drugs', 'Associated Drugs')
rows1 = []

labels2 = ('Side-Effect', 'TP', 'TN', 'FP', 'FN', 'Sensitivity', 'Specificity',
'Precision', 'Truth_ToxData_Drugs', 'Truth_Model_Drugs',
 'Pred_PathFX_Drugs', 'Pred_Model_Drugs',
'PathFX_Genes_all', 'PathFX_Genes', ' # of Genes_with_True_Pathways',
'# of Genes_with_False_Pathways',
'# of Shared_Genes_TP_FP_Pathways', '#_Distinct_TP_Genes', '#_Distinct_FP_Genes')
rows2 = []

labels3 = ('Side-Effect', '#_TP_Drugs', 'TP_Drugs', '#_TP_Genes', 'Sorted_TP_Genes',
'#_FP_Drugs', 'FP_Drugs', '#_FP_Genes', 'Sorted_FP_Genes', '#_Shared_Genes', 'Shared_Genes',
'#_Distinct_TP_Genes', 'Distinct_TP_Genes', '#_Distinct_FP_Genes', 'Distinct_FP_Genes')
rows3 = []

for item in se_name:
  eval_se2d(item)
  with open('PhenotypesAnalysisResults.csv', 'r') as read_obj1:
    csv_reader1 = reader(read_obj1)
    header1 = next(csv_reader1)
    if header1 != None:
      for row1 in csv_reader1:
        rows1.append(row1)
  with open('SideeffectMetrics.csv', 'r') as read_obj2:
    csv_reader2 = reader(read_obj2)
    header2 = next(csv_reader2)
    if header2 != None:
      for row2 in csv_reader2:
        rows2.append(row2)
  with open('GenesResults.csv', 'r') as read_obj3:
    csv_reader3 = reader(read_obj3)
    header3 = next(csv_reader3)
    if header3 != None:
      for row3 in csv_reader3:
        rows3.append(row3)

  print('')
  print('##########')
  print('')

dfmetall1 = pd.DataFrame.from_records(rows1, columns = labels1)
dfmetall2 = pd.DataFrame.from_records(rows2, columns = labels2)
dfmetall3 = pd.DataFrame.from_records(rows3, columns = labels3)
os.chdir(os.path.join(main_dir, 'Eval_SE2D_Res/'))
dfmetall1.to_csv('./all_PhenotypesAnalysisResults.csv', index=False)
dfmetall2.to_csv('./all_SideeffectMetrics.csv', index=False)
dfmetall3.to_csv('./all_GenesResults.csv', index=False)

print('')
print('##################################################')

##################################################
