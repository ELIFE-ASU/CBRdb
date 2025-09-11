# CBRdb: A Curated Biochemical Reaction database for precise Biochemical Reaction Analysis
[![DOI](https://zenodo.org/badge/804095458.svg)](https://doi.org/10.5281/zenodo.14948472) 

We present CBR-db, a curated biochemical database that integrates and refines data from KEGG and ATLAS databases to support precise analyses of biochemical reaction data. This curation effort addresses key limitations, such as the presence of malformed chemical representations, inaccurate stoichiometry, and ambiguous or incomplete reaction entries. These key improvements in CBR-db could facilitate the accurate analysis of physicochemical and thermodynamic properties, which are essential for modeling the evolution of chemical reactions in research areas ranging from prebiotic chemistry to metabolic network expansion. Altogether, CBR-db features a high number of high-quality reactions and compounds. The CBR-db is designed to be continuously updated, incorporating the latest releases from the KEGG and ATLAS databases. Furthermore, it provides a framework so the network can be extended and further issues can be addressed. 


To start everything locally, you can look at run_all.py, which will take you through the whole data pipeline. Otherwise, you can run each part by itself.

For further information, see the preprint on [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/67c28c046dde43c908f7aa37).


## Installation
### Fresh environment
It is recommended that you start from a fresh environment to prevent issues.
```
conda create -n cbrdb_env python=3.13
```
Activate the new env.
```
conda activate cbrdb_env
```
Add the conda-forge channel.
```
conda config --env --add channels conda-forge
```
Best to make them strict
```
conda config --set channel_priority true
```
Make sure to upgrade the conda env to force the channel priority.
```
conda update conda --all -y
```

### Install the requirements
```
conda install numpy sympy matplotlib networkx pandas rdkit chempy requests urllib3 chemparse ase pymatgen -y
```

### Optional extras
The main optional extra to install is ORCA and/or MACE. You will need these if you are interested in doing ab initio chemistry calculations.
For ORCA, head to their downloads [page](https://orcaforum.kofo.mpg.de/app.php/dlext/?view=detail&df_id=251).
For MACE, you will need to make sure you have PyTorch. Head to the [official PyTorch installation](https://pytorch.org/get-started/locally/) instructions page. An example might look like:
```
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu129
```
Then, proceed to install MACE. When you first run MACE, downloading a model might take a while.
```
pip3 install mace-torch
```
MACE offers a massive speed-up on GPUs but can run on a CPU. 
For more thermodynamic calculations, we use equilibrator-api.
```
conda install equilibrator-api
```

### CBR-db install
Then install CBRdb:
```
pip3 install git+https://github.com/ELIFE-ASU/CBRdb.git
```

## References
If you use this code in your work, you must reference the following:

- Hafner, J., MohammadiPeyhani, H., Sveshnikova, A., Scheidegger, A., & Hatzimanikatis, V. (2020). Updated ATLAS of biochemistry with new metabolites and improved enzyme prediction power. ACS synthetic biology, 9(6), 1479-1482.

- Hadadi, N., Hafner, J., Shajkofci, A., Zisaki, A., & Hatzimanikatis, V. (2016). ATLAS of biochemistry: a repository of all possible biochemical reactions for synthetic biology and metabolic engineering studies. ACS synthetic biology, 5(10), 1155-1166.

- Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.

Please take a look at the .bib file.
