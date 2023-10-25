# ToxPathEngine

## Overview

This repository, **ToxPathEngine**, contains the code and associated resources for the research presented in the manuscript titled "**Preclinical Side Effect Prediction through Pathway Engineering of Protein Interaction Network Models**". The research focuses on enhancing drug side effect predictions using protein interaction networks and pathway engineering, by including key pathway genes and omics data measurements. **ToxPathEngine** uses another tool called [**PathFX**](https://github.com/jenwilson521/PathFX), which is a protein-protein interaction tool. **ToxPathEngine** can identify distinct pathways, define novel gene pathways, and generate PathFX results using those new pathways along with the side effect evaluation metrics.

## Table of Contents

- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Citation](#Citation)

## Getting Started

These instructions will help you get a copy of the project up and running on your local machine for development and testing purposes. See [Usage](#usage) for notes on how to use the code in a research context.

## Prerequisites

Prerequisites that are necessary to run the code:

- **ToxPathEngine** was developed using the Python programming language.
- For setup guidance, please visit [**PathFX**](https://github.com/jenwilson521/PathFX) and follow the provided instructions.

## Usage

Detailed instructions on how to use the code:

- Clone this repository to your local machine/cluster.
- You can use the text file in the "data" folder. It contains the dataset we used for our analysis, the drug toxicity dataset, consisting of pairs of drugs and their associated side effects obtained from drug labels.
- To be able to run the last version of PathFX made in our analyses, first, you need to clone [**PathFX**](https://github.com/jenwilson521/PathFX). Afterward, you should add (copy/paste) the files available in the "pathfx" folder here in this GitHub repository to the same folder names (scripts/rscs/results) in your cloned PathFX folder on your local drive. Subsequently, you can use the "runpathfx_scr.py" script in our "scripts" folder to run the last version of PathFX on your operating system and re-generate the results.
-  The "scripts" folder includes all the scripts needed to re-generate our analyses:
   - "map_scr.py": Map drugs and map/match side effects.
   - "runpathfx_scr.py": Run the last version of PathFX to generate the results of the paper.
   - "sepred_scr.py": Evaluate per side effect in the baseline analysis (calculate the sensitivity and specificity values) and identify the key pathway genes.
   - "defpath_scr.py": Define novel pathways using the distinct identified genes and omics data in addition to the old associated pathways.
   - "evalnewpath_scr.py": Evaluate (the novel defined) pathways per side effect and their corresponding phenotypes and produce the evaluation plots.
- Important Note: Make sure to update the directory paths in all scripts to match your local environment before running them.
- We parsed gene signature data from the PharmOmics dataset which was acquired from the authors [**Chen et al., Iscience, 2022**](http://mergeomics.research.idre.ucla.edu/runpharmomics.php). We have included our Python code for extracting signature lists, but this requires the user to have obtained the original data.
- This analysis also requires multiple files related to DrugBank which can be obtained from [**DrugBank**](https://go.drugbank.com/releases). We downloaded release version 5.1.6, and parsed the XML database using the [**DrugBank/Parse**](https://github.com/dhimmel/drugbank/blob/gh-pages/parse.ipynb) notebook as a guide. Specifically, the user will need a pickled dictionary "drugbankid_to_name.pkl".
- You can find the PathFX outcome for the drug Alteplase, as an example, in the "pathfx/results/" directory, which represents our final analysis results.

## Citation

If you use this code or the associated research in your work, please consider citing our manuscript:

Alidoost, Mohammadali and Wilson, L. Jennifer, "Preclinical Side Effect Prediction through Pathway Engineering of Protein Interaction Network Models", Submitted (2023).
