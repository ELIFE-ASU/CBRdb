# CBRdb: A Curated Biochemical Reaction database for precise Biochemical Reaction Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17459478.svg)](https://doi.org/10.5281/zenodo.17459478)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ELIFE-ASU/CBRdb/HEAD)

## 🗺️ Overview

We present CBR-db, a curated biochemical database that integrates and refines data from KEGG and ATLAS databases to
support precise analyses of biochemical reaction data. This curation effort addresses key limitations, such as the
presence of malformed chemical representations, inaccurate stoichiometry, and ambiguous or incomplete reaction entries.
These key improvements in CBR-db could facilitate the accurate analysis of physicochemical and thermodynamic properties,
which are essential for modeling the evolution of chemical reactions in research areas ranging from prebiotic chemistry
to metabolic network expansion. Altogether, CBR-db features a high number of high-quality reactions and compounds. The
CBR-db is designed to be continuously updated, incorporating the latest releases from the KEGG and ATLAS databases.
Furthermore, it provides a framework so the network can be extended and further issues can be addressed.

To start everything locally, you can look at run_all.py, which will take you through the whole data pipeline. Otherwise,
you can run each part by itself. To sidestep installation
checkout [Binder](https://mybinder.org/v2/gh/ELIFE-ASU/CBRdb/HEAD).

For further information, see the preprint
on [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/67c28c046dde43c908f7aa37).

## 📋 Features

- Integration of KEGG and ATLAS biochemical reaction databases.
- Curation of compounds and reactions.
- Thermodynamic properties.
- Reaction pruning based on reaction feasibility and similarity.

Example fixes and receipts of changes can be found in the data directory. See for example `C_IDs_bad.dat`
for a list of KEGG compounds that were removed or adjusted.

## 💡 Getting Started

To get started with the data, you can download the compounds `CBRdb_C.csv` and reactions `CBRdb_R.csv` from the
repository. For convenience, the files are also hosted on Zenodo and are stored as a compressed format to reduce
download time.

```
import pandas as pd
import os

# Download CBRdb compounds
os.system('wget https://github.com/ELIFE-ASU/CBRdb/blob/main/CBRdb_C.csv.zip')

# Load CBRdb compounds
compounds = pd.read_csv('CBRdb_C.csv.zip')
compounds.head()

# Download CBRdb reactions
os.system('wget https://github.com/ELIFE-ASU/CBRdb/blob/main/CBRdb_R.csv.zip')

# Load CBRdb reactions
reactions = pd.read_csv('CBRdb_R.csv.zip')
reactions.head()
```

Below is a breakdown description of the columns in each of the data files. Each description table is collapsible for
easier viewing, but the subcategories are shown for clarity. In the tables below, the `Entry Title` column corresponds
to the column headers
in the CSV files and the `Description` column provides a brief explanation of each entry.
<details>
<summary>CBRdb_C</summary>
<br>

General information:

| Entry Title      | Description                                                |
|------------------|------------------------------------------------------------|
| compound_id      | Unique identifier for each compound in the CBRdb database. |
| smiles           | SMILES representation of the compound.                     |
| inchi            | InChI representation of the compound.                      |
| formula          | Molecular formula of the compound.                         |
| molecular_weight | Molecular weight of the compound.                          | 
| n_heavy_atoms    | number of heavy atoms in the compound                      |
| n_chiral_centers | number of chiral centers in the compound.                  |
| formal_charge    | Formal charge of the compound.                             |
| smiles_capped    | SMILES representation with hydrogen capped R groups.       |
| inchi_capped     | InChI representation with hydrogen capped R groups.        |
| nickname         | Common alias for the compound; first name in KEGG.         |
| comment          | KEGG's comments with context about the compound.           |
| CBRdb_R_ids      | List of reaction IDs in CBRdb that involve this compound.  |

Complexity scores and indices:

| Entry Title     | Description                             |
|-----------------|-----------------------------------------|
| assembly_index  | Assembly index of the compound.         |
| wiener_index    | Wiener index of the compound.           | 
| unique_bonds    | Number of unique bonds in the compound. |
| spacial_score   | Spacial score of the compound.          |
| randic_index    | Randic index of the compound.           |
| proudfoot       | Proudfoot score of the compound.        |                                                                                                               
| mc2             | Molecular complexity score 2.           |
| mc1             | Molecular complexity score 1.           |
| kirchhoff_index | Kirchhoff index of the compound.        |
| fcfp4           | FCFP4 fingerprint of the compound.      |
| bertz           | Bertz complexity of the compound.       |
| balaban_index   | Balaban index of the compound.          |

Properties:

| Entry Title       | Description                                                              |
|-------------------|--------------------------------------------------------------------------|
| ionization_states | Ionization states of the compound. For a pH range of 4 to 10.            |
| std_dgf           | Free energy of formation using experimental values.                      |
| std_dgf_error     | Error in the free energy of formation using experimental values.         |
| std_dgf_p         | Free energy of formation at pressure using experimental values.          |
| std_dgf_p_error   | Error in free energy of formation at pressure using experimental values. |
| d_free            | Free energy of formation using quantum methods.                          |
| d_enthalpy        | Enthalpy of formation using quantum methods.                             |
| d_entropy         | Entropy of formation using quantum methods.                              |
| free              | Free energy using quantum methods.                                       |
| enthalpy          | Enthalpy using quantum methods.                                          |
| entropy           | Entropy energy using quantum methods.                                    |
| vib_energies      | List of vibration modes using quantum methods.                           |

KEGG specific information:

| Entry Title     | Description                                                |
|-----------------|------------------------------------------------------------|
| kegg_type       | Type for peptide or polyketide entry.                      |
| kegg_sequence   | Sequence for peptide or polyketide entry.                  |
| kegg_reaction   | KEGG REACTION entries associated with the compound.        |
| kegg_pathway    | KEGG PATHWAY entries associated with the compound.         |
| kegg_organism   | For peptides or polyketides, KEGG organism with sequence.  |
| kegg_network    | KEGG NETWORK entries associated with the compound.         |
| kegg_mol_weight | Molecular weight provided by KEGG.                         |
| kegg_module     | KEGG MODULE entries associated with the compound.          |
| kegg_glycan     | KEGG GLYCAN entries associated with the compound.          |
| kegg_gene       | KEGG GENES entries associated with the compound.           |
| kegg_formula    | Molecular formula provided by KEGG.                        |
| kegg_exact_mass | Exact (single-isotope) mass provided by KEGG.              |
| kegg_enzyme     | KEGG ENZYME entries associated with the compound.          |
| kegg_drug       | KEGG DRUG entries associated with the compound.            |
| kegg_brite_full | Full KEGG BRITE hierarchies with context, all levels.      |
| kegg_brite      | KEGG BRITE hierarchy entries associated with the compound. |

Database cross-references:

| Entry Title | Description                                    |
|-------------|------------------------------------------------|
| PubChem     | Linked PubChem identifiers for the compound.   |
| PDB_CCD     | Linked PDB CCD identifiers for the compound.   |
| NIKKAJI     | Linked NIKKAJI identifiers for the compound.   |
| LIPIDMAPS   | Linked LIPIDMAPS identifiers for the compound. |
| KNApSAcK    | Linked KNApSAcK identifiers for the compound.  |
| Drug_group  | Linked KEGG DGROUP entries for the compound.   |
| ChEBI       | Linked ChEBI identifiers for the compound.     |
| CAS         | Linked CAS identifiers for the compound.       |
| ATC_code    | Linked ATC codes of the compound.              |

</details>


<details>
<summary>CBRdb_R</summary>
<br>

General information:

| Entry Title | Description                        |
|-------------|------------------------------------|
| id          | Reaction identifier.               |
| reaction    | Reaction equation.                 |
| ec          | Enzyme Commission numbers.         |
| name        | Reaction name.                     |
| smarts      | SMARTS representation.             |
| reac_sim    | Closest reaction similarity.       |
| rhea        | Rhea identifiers if provided.      |
| CBRdb_C_ids | Corresponding CBRdb C identifiers. |
| id_orig     | Original reaction identifier(s).   |

KEGG specific information:

| Entry Title | Description                       |
|-------------|-----------------------------------|
| module      | KEGG MODULE identifiers.          |
| orthology   | KEGG ORTHOLOGY identifiers.       |
| pathway     | KEGG PATHWAY identifiers.         |
| rclass      | KEGG RCLASS IDs & COMPOUND pairs. |
| remark      | KEGG remark.                      |
| comment     | Unstructured KEGG comment field.  |
| overall     | Flag: Overall reaction (br08210). |

ATLAS specific information:

| Entry Title   | Description                          |
|---------------|--------------------------------------|
| most_sim_kegg | Most similar KEGG reaction (MSK). \* |
| bridgit_score | BridgIT score for most_sim_kegg. \*  |

Data flags:

| Entry Title             | Description                          |
|-------------------------|--------------------------------------|
| flags                   | From CBRdb.df_of_suspect_reactions   |
| balancer_failed         | Flag: Balancer failed.               |
| bool_missing_data       | Flag: Missing structure(s).          |
| bool_var_list           | Flag: Variable in coefficients.      |
| cpd_starred             | Flag: Has starred compounds.         |
| els_and_stars_balance   | Flag: Elements & stars are balanced  |

\* Provided by ATLAS;
see [ATLAS User Guide](https://lcsb-databases.epfl.ch/pathways/atlas/files/ATLAS_UserGuide.pdf#page=6.00) for
most_similar_kegg header details.
</details>

To get started with deeper components of CBRdb, follow the installation instructions below. After installation, you can
start using the CBRdb package in your Python scripts or Jupyter notebooks.

## 🔧 Installation

Here you can find the full installation instructions. See
also [Binder](https://mybinder.org/v2/gh/ELIFE-ASU/CBRdb/HEAD).

<details>
<summary>Using conda env</summary>
<br>

Using conda is the recommended way to install the required packages to help with dependencies. First, lets clone the
repo so that we have a local copy.

```
git clone https://github.com/ELIFE-ASU/CBRdb.git
```

Change into the CBRdb directory.

```
cd CBRdb
```

Change into the build tools directory.

```
cd build_tools
```

Create the conda environment.

```
conda env create -f environment.yml
```

Activate the conda environment.

```
conda activate cbrdb
```

Lets leave the build tools directory.

```
cd ..
```

Install CBRdb.

```
pip install -e .
```

If you want to install manually, follow the instructions below.
</details>

<details>
<summary>Manually</summary>
<br>

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

The main optional extra to install is ORCA and/or MACE. You will need these if you are interested in doing ab initio
chemistry calculations.
For ORCA, head to their downloads [page](https://orcaforum.kofo.mpg.de/app.php/dlext/?view=detail&df_id=251).
For MACE, you will need to make sure you have PyTorch. Head to
the [official PyTorch installation](https://pytorch.org/get-started/locally/) instructions page. An example might look
like:

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

</details>

## ⚠️ Issue Tracker

Found a bug? Have an enhancement request? Head over to the [GitHub issue
tracker](https://github.com/ELIFE-ASU/CBRdb/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as possible about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

## 👥 Contributors

Louie Slocombe, Camerian Millsaps, M.Reza Shahjahan, Kamesh Narasimhan, and Sara Walker.

## ⚖️ License

MIT License. We ask that you cite the relevant papers, please!

## 📚 References

If you use this code in your work, you must reference the following:

- Hafner, J., MohammadiPeyhani, H., Sveshnikova, A., Scheidegger, A., & Hatzimanikatis, V. (2020). Updated ATLAS of
  biochemistry with new metabolites and improved enzyme prediction power. ACS synthetic biology, 9(6), 1479-1482.

- Hadadi, N., Hafner, J., Shajkofci, A., Zisaki, A., & Hatzimanikatis, V. (2016). ATLAS of biochemistry: a repository of
  all possible biochemical reactions for synthetic biology and metabolic engineering studies. ACS synthetic biology, 5(
  10), 1155-1166.

- Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.

Please take a look at the .bib file.
