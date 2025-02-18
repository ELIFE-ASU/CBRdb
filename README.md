# CBRdb: A Curated Biochemical Reaction database for precise Biochemical Reaction Analysis

We present CBR-db, a curated biochemical database that integrates and refines data from KEGG and ATLAS databases to support precise analyses of biochemical reaction data. This curation effort addresses key limitations, such as the presence of malformed chemical representations, inaccurate stoichiometry, and ambiguous or incomplete reaction entries. These key improvements in CBR-db could facilitate accurate analysis of physicochemical and thermodynamic properties, which are essential for modeling the evolution of chemical reactions in research areas ranging from prebiotic chemistry to metabolic network expansion. Altogether, CBR-db features a high number of high-quality reactions and compounds. CBR-db is designed to be continuously updated, incorporating the latest releases from KEGG and ATLAS databases. Furthermore, it provides a framework so that the network can be extended and further issues can be improved.


To start everything locally, you can look at run_all.py, which will take you over the whole data pipeline. Otherwise, you can run each part by themselves.


## Installation
### Fresh environment
It is recommended that you start from a fresh environment to prevent issues.
```
conda create -n cbrdb_env python=3.12
```
Activate the new env.
```
conda activate cbrdb_env
```
Add channels in this order.
```
conda config --env --add channels conda-forge
```
Best to make them strict
```
conda config --set channel_priority true
```
To check your updated channel list, run:
```
conda config --show channels
```
Make sure to upgrade the conda env to force the channel priority.
```
conda update conda --all -y
```
### Install the requirements
```
conda install conda-forge::numpy conda-forge::sympy conda-forge::matplotlib conda-forge::networkx conda-forge::pandas conda-forge::rdkit conda-forge::chempy conda-forge::requests conda-forge::urllib3 conda-forge::chemparse conda-forge::swifter -y
```
Then install CBRdb:
```
pip install git+https://github.com/ELIFE-ASU/CBRdb.git@v1.0.0
```

## References
If you use this code in your work, you must reference the following:

- Hafner, J., MohammadiPeyhani, H., Sveshnikova, A., Scheidegger, A., & Hatzimanikatis, V. (2020). Updated ATLAS of biochemistry with new metabolites and improved enzyme prediction power. ACS synthetic biology, 9(6), 1479-1482.

- Hadadi, N., Hafner, J., Shajkofci, A., Zisaki, A., & Hatzimanikatis, V. (2016). ATLAS of biochemistry: a repository of all possible biochemical reactions for synthetic biology and metabolic engineering studies. ACS synthetic biology, 5(10), 1155-1166.

- Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.

See .bib file.

## Examples of issues with KEGG and ATLAS
Here are a couple of example problems that this code attempts to solve.

### Compounds
- Cases where the R group notation is unstandardized, for example, R#, R, *, is used interchangeably.
- Case with no mol file: C00028, C00030, C00505, C00923, ... etc
- Case with broken mol files, for example, in some cases, OH is an element: C20442, C21011, C21012, C21013, C22197, ... etc.
- Halogen generalized mol file: C00462, C01322, C01365, C01706, ... etc

### Reactions
- Unbalanced; nothing is missing, but the elemental stoichiometry is incorrect
- One or more elements are missing from a side of the equation: , ... etc
- One or more compounds are missing from a side of the equation: , ... etc
- Zero population (for ATLAS KEGG FILE): R04242, R04254, ... etc
- Reactions that are part of shortcuts and are represented by chains of other reactions: , ... etc
- Unclear reactions: R10743, ... etc
- reactions with unclear compounds: , ... etc
