# CBRdb: A Curated Biochemical Reaction database for precise Biochemical Reaction Analysis

[![DOI](https://zenodo.org/badge/804095458.svg)](https://doi.org/10.5281/zenodo.14948472) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ELIFE-ASU/CBRdb/HEAD)

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
you can run each part by itself.

For further information, see the preprint
on [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/67c28c046dde43c908f7aa37).

## 📋 Features

- Integration of KEGG and ATLAS biochemical reaction databases.
- Curation of compounds and reactions.
- Thermodynamic properties.
- Reaction pruning based on reaction feasibility and similarity.

## 💡 Getting Started

To get started with the data, you can download the compounds `CBRdb_C.csv` and reactions `CBRdb_R.csv` from the
repository.

To get started with deeper components of CBRdb, follow the installation instructions below. After installation, you can
start using the
CBRdb package in your Python scripts or Jupyter notebooks.

<details>
<summary>CBRdb_C</summary>
<br>

| Syntax            | Description                                                |
|-------------------|------------------------------------------------------------|
| compound_id       | Unique identifier for each compound in the CBRdb database. |
| smiles            | SMILES representation of the compound.                     |
| formula           | Molecular formula of the compound.                         |
| molecular_weight  | Molecular weight of the compound.                          |
| n_heavy_atoms     | number of heavy atoms in the compound                      |
| n_chiral_centers  | number of chiral centers in the compound.                  |
| smiles_capped     | SMILES representation with hydrogen capped R groups.       |
| inchi_capped      | InChI representation with hydrogen capped R groups.        |
| nickname          | Common nickname or alias for the compound.                 |
| comment           | Additional comments or notes about the compound.           |
| wiener_index      | Wiener index of the compound.                              |
| unique_bonds      | Number of unique bonds in the compound.                    |
| spacial_score     | Spacial score of the compound.                             |
| randic_index      | Randic index of the compound.                              |
| proudfoot         | Proudfoot score of the compound.                           |
| name              | Common name of the compound.                               |
| mc2               | Molecular complexity score 2.                              |
| mc1               | Molecular complexity score 1.                              |
| kirchhoff_index   | Kirchhoff index of the compound.                           |
| kegg_type         | KEGG type of the compound.                                 |
| kegg_sequence     | KEGG sequences of the compound.                            |
| kegg_reaction     | KEGG reactions associated with the compound.               |
| kegg_pathway      | KEGG pathways associated with the compound.                |
| kegg_organism     | KEGG organisms associated with the compound.               |
| kegg_network      | KEGG networks associated with the compound.                |
| kegg_mol_weight   | KEGG molecular weight of the compound.                     |
| kegg_module       | KEGG modules associated with the compound.                 |
| kegg_glycan       | KEGG glycans associated with the compound.                 |
| kegg_gene         | KEGG genes associated with the compound.                   |
| kegg_formula      | KEGG formula of the compound.                              |
| kegg_exact_mass   | KEGG exact mass of the compound.                           |
| kegg_enzyme       | KEGG enzymes associated with the compound.                 |
| kegg_drug         | KEGG drugs associated with the compound.                   |
| kegg_brite_full   | KEGG (full) brite information of the compound.             |
| kegg_brite        | KEGG brite information of the compound.                    |
| ionization_states | Ionization states of the compound.                         |
| inchi             | InChI representation of the compound.                      |
| formal_charge     | Formal charge of the compound.                             |
| fcfp4             | FCFP4 fingerprint of the compound.                         |
| bertz             | Bertz complexity of the compound.                          |
| balaban_index     | Balaban index of the compound.                             |
| PubChem           | PubChem identifier of the compound.                        |
| PDB_CCD           | PDB CCD identifier of the compound.                        |
| NIKKAJI           | NIKKAJI identifier of the compound.                        |
| LIPIDMAPS         | LIPIDMAPS identifier of the compound.                      |
| KNApSAcK          | KNApSAcK identifier of the compound.                       |
| Drug_group        | Drug group of the compound.                                |
| ChEBI             | ChEBI identifier of the compound.                          |
| CAS               | CAS identifier of the compound.                            |
| ATC_code          | ATC code of the compound.                                  |
| CBRdb_R_ids       | List of reaction IDs in CBRdb that involve this compound.  |

</details>


<details>
<summary>CBRdb_R</summary>
<br>

| Syntax                  | Description                        |
|-------------------------|------------------------------------|
| id                      | Reaction identifier.               |
| reaction                | Reaction equation.                 |
| ec                      | Enzyme Commission number.          |
| module                  | Metabolic module.                  |
| orthology               | Orthology identifiers.             |
| pathway                 | Pathway information.               |
| rclass                  | Reaction classification.           |
| rhea                    | Rhea identifier.                   |
| balancer_failed         | Balancer failed flag.              |
| bool_missing_data       | Missing data flag.                 |
| bool_var_list           | Variable list flag.                |
| bridgit_score           | Bridgit score.                     |
| comment                 | Comments.                          |
| cpd_starred             | Starred compounds.                 |
| flags                   | Flags.                             |
| id_orig                 | Original reaction identifier.      |
| is_balanced_except_star | Balanced except starred flag.      |
| kegg_id                 | KEGG reaction identifier.          |
| most_sim_kegg           | Most similar KEGG reaction.        |
| msk_ecs                 | -                                  |
| msk_metacyc             | -                                  |
| msk_mnxr                | -                                  |
| msk_rhea                | -                                  |
| msk_rns                 | -                                  |
| name                    | Reaction name.                     |
| overall                 | Overall information.               |
| remark                  | Remarks.                           |
| smarts                  | SMARTS representation.             |
| CBRdb_C_ids             | Corresponding CBRdb C identifiers. |

</details>

## 🔧 Installation

<details>
<summary>Using conda env</summary>
<br>

Using conda is the recommended way to install the required packages.

```
git clone https://github.com/ELIFE-ASU/CBRdb.git
```

Change into the CBRdb directory.

```
cd CBRdb
```

Create the conda environment.

```
conda env create -f environment.yml
```

Activate the conda environment.

```
conda activate cbrdb
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

Louie Slocombe, Camerian Millsaps, Reza Shahjahan, and Kamesh Narasimhan.

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
