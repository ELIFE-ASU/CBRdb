import pandas as pd
import os

from .tools_files import reaction_csv, compound_csv, space_sep_str_cols_cps

out_fmt = {'encoding': 'utf-8', 'index': False}


def id_indexed(df: pd.DataFrame) -> pd.DataFrame:
    """ Convenience fxn: Sets 'id' or 'compound_id' column as a pd.DataFrame's index. Returns a copy (not modified inplace) """
    if df.index.name in ['id', 'compound_id']:
        return df
    elif 'id' in df.columns and df.index.name != 'id':
        return df.set_index('id')
    elif 'compound_id' in df.columns and df.index.name != 'compound_id':
        return df.set_index('compound_id')
    else:
        return df


def identify_duplicate_compounds(C_main):
    """
    Identifies duplicate compounds in a DataFrame and maps them to new unique IDs.

    Parameters:
    C_main (pd.DataFrame): The main DataFrame containing compound data with 'compound_id' and 'smiles' columns.

    Returns:
    pd.DataFrame: A DataFrame mapping old compound IDs to new unique compound IDs.
    """
    # ID the number to start counting back from when instantiating new IDs
    id_num = C_main.reset_index()['compound_id'].str.lstrip('C').astype(int)
    count_back_from = id_num.loc[id_num.diff().idxmax()] - 1
    # Flag duplicates among well-defined structures
    possible_dupes = C_main.dropna(subset='smiles').query(
        'smiles.duplicated(keep=False) & ~smiles.str.contains("*", regex=False)')
    # Group duplicate structures together
    C_dupemap = possible_dupes.groupby(by='smiles')['compound_id'].apply(sorted)
    # Assign a new ID for each dupe-group
    C_dupemap = C_dupemap.reset_index(drop=True).explode().rename('old_id')
    C_dupemap.index = ('C' + (count_back_from - C_dupemap.index).astype(str)).set_names('new_id')
    # Sort the duplicate compound groups (within and between)
    C_dupemap = C_dupemap.reset_index().sort_values(by=['new_id', 'old_id']).set_index('old_id')
    # Generate dupe-map output file for easy conversion from old IDs to new IDs
    C_dupemap.to_csv('../data/kegg_data_C_dupemap.csv', encoding='utf-8')
    return C_dupemap


def merge_duplicate_compounds(C_main: pd.DataFrame, C_dupemap: pd.DataFrame) -> pd.DataFrame:
    """
    Merges duplicate compounds in a compound DataFrame using the duplication-map DataFrame.
    Combines the lists of alternate IDs and relationships for all duplicated compounds.

    Parameters:
    C_main (pd.DataFrame): A DataFrame with compound data and compound_id column or index.name
    C_dupemap (pd.DataFrame): A DataFrame mapping duplicate compounds (old_id index) to their merged new_id.

    Returns:
    pd.DataFrame: a copy of the main compound DataFrame, but with the duplicate compounds merged.
    """
    # define functions for combining entries
    line_set_str = lambda x: '; '.join(sorted(set(x)))  # unique strings, sorted
    lines_set_str = lambda x: ';~'.join(dict.fromkeys(x.str.split(';~').sum()).keys())
    name_funcs = {'kegg_name': lines_set_str, 'nickname': line_set_str}
    tripleslash_join = lambda x: '///'.join(sorted(set(x.dropna()))) if not x.dropna().empty else None
    list_unique = lambda x: ' '.join(sorted(set(x.dropna().str.split(' ').explode()))) if not x.dropna().empty else None

    # standardize data format
    C_main_copy = id_indexed(C_main.copy(deep=True))
    # Assign new IDs to dupes
    duped_compounds = C_main_copy.loc[C_dupemap.index].assign(compound_id=C_dupemap).reset_index()
    # Flag presence/absence of repeated subunit indicator "n" in kegg_formula
    duped_compounds['n'] = duped_compounds['kegg_formula'].str.contains('n', na=False, regex=False)
    # Group by new IDs
    dupes = duped_compounds.groupby('compound_id')
    # Identify which fields do/don't have conflicting values
    n_vals = dupes.nunique()
    # Take the first non-NaN entry in fields/IDs without conflicting values
    group_attrs = dupes.first()[n_vals < 2].drop(columns=['old_id', 'n'])
    # For fields where values reflect space-separated lists...
    sscol = group_attrs.columns.intersection(set(space_sep_str_cols_cps.keys()))
    # Take the union of everything listed
    dupes_ss = dupes[sscol].agg(list_unique)[n_vals > 1]
    group_attrs.fillna(dupes_ss, inplace=True)

    # Procedurally combine nickname + name fields
    group_attrs.fillna(dupes.agg(name_funcs), inplace=True)

    # For the kegg_brite_full and comment columns...
    label_and_merge = ['comment', 'kegg_brite_full']
    for col in label_and_merge:
        # Store provenance to ease interpretation
        duped_compounds.update((duped_compounds['old_id'] + ': ' + duped_compounds[col]).rename(col))
        # Then merge those fields
        group_attrs.fillna(duped_compounds.groupby('compound_id')[col].apply(tripleslash_join), inplace=True)

    # The only columns left w/conflicts should be 'kegg_mol_weight', 'kegg_exact_mass', and 'kegg_formula'
    one_val = n_vals.columns[n_vals.lt(2).all()]
    cols_done = sscol.union(label_and_merge).union(one_val).union(name_funcs.keys())
    # Only one dupe group (C98917) has multiple values for mass/weight.
    # For self-consistency WRT kegg_formula (also multiple values), get its first entry whole-cloth.
    first_weight = id_indexed(dupes[['compound_id', 'kegg_mol_weight', 'kegg_exact_mass', 'kegg_formula']].head(1))[
        n_vals > 1].dropna(how='any')
    group_attrs.fillna(first_weight, inplace=True)

    # Prioritize retaining kegg_formula values with "n" in them
    sorted_formulas = duped_compounds.sort_values(by=['n', 'old_id'], ascending=[False, True]).groupby('compound_id')
    sorted_formulas = sorted_formulas[['kegg_formula']].first().drop(first_weight.index)[n_vals > 1].dropna()
    group_attrs.fillna(sorted_formulas, inplace=True)

    # Ensure no columns remain un-merged
    unmerged = C_main_copy.columns.difference(
        cols_done.union(['n', 'old_id', 'kegg_formula', 'kegg_mol_weight', 'kegg_exact_mass']))
    if len(unmerged) > 0:
        print("Warning: The following columns were not merged:", unmerged.tolist())
        print("Assuming first non-NaN entry is acceptable.")
        group_attrs.fillna(dupes[unmerged].first(), inplace=True)

    # Now, remove the duplicate entries from the main DataFrame and add the merged entries
    C_main_copy = pd.concat([C_main_copy.drop(C_dupemap.index), group_attrs])
    # sort by compound ID
    C_main_copy.sort_index(inplace=True)
    # return index to initial state
    C_main_copy.reset_index(inplace=True)

    return C_main_copy


def add_R_col_to_C_file(reaction_df : pd.DataFrame, c_file='../CBRdb_C_metadata.csv.zip'):
    """
    Adds a column to the compound metadata indicating which reactions each compound is involved in.
    Parameters:
    reaction_df (pd.DataFrame): Reaction DataFrame.
    c_file (str): File path for the compound metadata DataFrame.
    Returns:
    None: The function modifies the compound and reaction DataFrames in place and saves them to the specified files. 
    """
    print("Adding column in compound metadata to indicate associated reactions", flush=True)
    c2r = list_reactions_per_compound(reaction_df).map(' '.join)
    join_col_to_csv(col=c2r, file=c_file, col_name='CBRdb_R_ids')
    print("Finished appending reaction lists to compound metadata", flush=True)
    return None


def merge_hpc_calculations(final_output_Cs_fp: str|pd.DataFrame = '../CBRdb_C.csv.zip',
                           formation_energies_fp='../hpc/CBRdb_C_formation_energies.csv.gz',
                           mace_spectrum_fp='../hpc/CBRdb_C_mace_spectrum.csv.gz',
                           assembly_index_fp='../hpc/CBRdb_C_assembly_index.csv.zip',
                           c_dupemap_fp='../data/kegg_data_C_dupemap.csv'):
    """ 
    Merges HPC calculations into the main compounds data, overwriting those columns if present.

    """
    print("Merging HPC calculations into compound file or DataFrame", flush=True)

    # Specify file reading and writing parameters
    f_params_in = dict(index_col=0, low_memory=False, dtype=space_sep_str_cols_cps)
    f_params_out = dict(encoding='utf-8', index=True, compression='infer')

    # Import the datasets
    if not isinstance(final_output_Cs_fp, (str, pd.DataFrame)):
        raise ValueError("final_output_Cs_fp must be one of (str, pd.DataFrame)")
    elif isinstance(final_output_Cs_fp, str):
        data_c = pd.read_csv(os.path.abspath(final_output_Cs_fp), **f_params_in)
    else:
        not_id_indexed = final_output_Cs_fp.index.astype(str).str.isnumeric().all()
        data_c = id_indexed(final_output_Cs_fp)

    data_c_dupes = pd.read_csv(os.path.abspath(c_dupemap_fp), **f_params_in)
    data_c_formation = pd.read_csv(os.path.abspath(formation_energies_fp), **f_params_in)
    data_c_mace = pd.read_csv(os.path.abspath(mace_spectrum_fp), **f_params_in)
    data_c_aix = pd.read_csv(os.path.abspath(assembly_index_fp), **f_params_in)

    # For MACE and AI, prepare to merge on SMILES string.
    smi2cp = {v: k for k,v in data_c['smiles'].items()}
    data_c_mace.set_index('smiles', inplace=True)
    data_c_aix.set_index('smiles', inplace=True)
    mace_eligible = data_c_mace.filter(items=data_c['smiles'], axis=0).rename(smi2cp)
    aix_eligible = data_c_aix[['assembly_index']].filter(items=data_c['smiles'], axis=0).rename(smi2cp)
    data_c_to_add = mace_eligible.join(aix_eligible, how='outer').rename_axis('compound_id')

    # For formation energies, prepare to merge on ID by syncing IDs.
    data_c_formation.rename(index=data_c_dupes.iloc[:,0], inplace=True)

    # Remove previous compound calculations, if present.
    data_c.drop(columns=data_c_to_add.columns, errors='ignore', inplace=True)
    data_c.drop(columns=data_c_formation.columns, errors='ignore', inplace=True)

    # Merge the datasets
    data_c = data_c.join(data_c_to_add, how='left').join(data_c_formation, how='left')

    # Save the updated compounds file (note that floats are not rounded unlike in CBRdb.compound_csv)
    data_c.to_csv('../CBRdb_C.csv.zip', **f_params_out)

    print("HPC calculation merger complete", flush=True)

    if isinstance(final_output_Cs_fp, pd.DataFrame):
        if not_id_indexed:
            return data_c.reset_index()
        else:
            return data_c
    else:
        return None


def merge_hpc_reaction_calculations(final_output_Rs_fp : str|pd.DataFrame = '../CBRdb_R.csv.zip',
                                    similarity_fp='../hpc/ATLAS_R_sim.csv.gz',
                                    thermo_params_fp='../hpc/CBRdb_R_reaction_energies.csv.gz'):
    """ 
    Merges reaction thermodynamic parameters into the reactions data file or DataFrame, overwriting it.

    """
    print("Merging HPC reaction calculations into reaction file", flush=True)
    f_params_in = dict(index_col=0, low_memory=False)
    f_params_out = dict(encoding='utf-8', index=True, compression='infer')
    
    # Import the datasets
    if not isinstance(final_output_Rs_fp, (str, pd.DataFrame)):
        raise ValueError("final_output_Rs_fp must be one of (str, pd.DataFrame)")
    elif isinstance(final_output_Rs_fp, str):
        reactions = id_indexed(pd.read_csv(os.path.abspath(final_output_Rs_fp), **f_params_in))
    else:
        not_id_indexed = final_output_Rs_fp.index.astype(str).str.isnumeric().all()
        reactions = id_indexed(final_output_Rs_fp)

    # Assume thermodynamic and similarity calculations were performed on the same CBRdb version
    thermo_params = id_indexed(pd.read_csv(thermo_params_fp, **f_params_in))
    similarity_params = id_indexed(pd.read_csv(similarity_fp, **f_params_in))
    hpc_calcs = similarity_params.join(thermo_params, how='outer')

    # Remove the field we're about to add, if present
    merging_on = ['reaction', 'smarts']
    reactions.drop(columns=hpc_calcs.columns.difference(merging_on), inplace=True, errors='ignore')

    # Merge the datasets
    reactions = id_indexed(reactions.reset_index().merge(hpc_calcs, on=merging_on, how='left'))
    print("HPC reaction calculation merger complete", flush=True)
                     
    if isinstance(final_output_Rs_fp, pd.DataFrame):
        if not_id_indexed:
            return reactions.reset_index()
        else:
            return reactions
    else:
        reactions.to_csv(final_output_Rs_fp **f_params_out)
        return None


def merge_atom_mapping_calculations(final_output_Rs_fp : str|pd.DataFrame = '../CBRdb_R.csv.zip',
                                    atom_tracking_fp='../atom_tracking/combined_output_all.csv.gz'):
    """ 
    Merges atom mapping into the reactions data file or DataFrame, overwriting existing columns if present.

    """
    print("Merging atom-mapping calculations into reaction file", flush=True)
    f_params_in = dict(index_col=0, low_memory=False)
    f_params_out = dict(encoding='utf-8', index=True, compression='infer')
    
    # Import the datasets
    if not isinstance(final_output_Rs_fp, (str, pd.DataFrame)):
        raise ValueError("final_output_Rs_fp must be one of (str, pd.DataFrame)")
    elif isinstance(final_output_Rs_fp, str):
        reactions = id_indexed(pd.read_csv(os.path.abspath(final_output_Rs_fp), **f_params_in))
    else:
        not_id_indexed = final_output_Rs_fp.index.astype(str).str.isnumeric().all()
        reactions = id_indexed(final_output_Rs_fp)

    # Import the atom-mapped reactions
    usecols = ['reaction_id', 'mapped_rxns', 'reaction_no_stoich']
    calc_df = pd.read_csv(atom_tracking_fp, usecols=usecols, **f_params_in)
    calc_df.rename(columns={'mapped_rxns': 'atom_mapping'}, inplace=True)

    # Make the field we're merging on
    calc_df['sorted_no_stoich'] = calc_df['reaction_no_stoich'].str.split(' <=> ').map(sorted).map(' <=> '.join)
    reactions['sorted_no_stoich'] = (reactions['reaction'].str.split(' <=> ', expand=True)
                          .apply(lambda x: x.str.findall(r'(C\d+)'))
                          .map(' + '.join)
                          .map(lambda x: [x])
                          .sum(axis=1)
                          .map(sorted)
                          .map(' <=> '.join))
    
    # Remove the field we're about to add, if present
    if 'atom_mapping' in reactions.columns:
        reactions.drop(columns='atom_mapping', inplace=True)

    # Merge datasets
    reactions = id_indexed(reactions.reset_index().merge(calc_df, on='sorted_no_stoich', how='left'))
    reactions.drop(columns=['sorted_no_stoich', 'reaction_no_stoich'], inplace=True)

    print("Atom-mapping calculation merger complete", flush=True)
                     
    if isinstance(final_output_Rs_fp, pd.DataFrame):
        if not_id_indexed:
            return reactions.reset_index()
        else:
            return reactions
    else:
        reactions.to_csv(final_output_Rs_fp **f_params_out)
        return None
    


def export_reaction_metadata(main_fp : str | pd.DataFrame = '../CBRdb_R.csv.zip',
                              meta_fp = '../CBRdb_R_metadata.csv.zip'):
    
    print("Organizing reaction data and metadata", flush=True)

    f_params_in = dict(index_col=0, low_memory=False)
    f_params_out = dict(encoding='utf-8', index=True, compression='infer')
    
    # Import the datasets
    if not isinstance(main_fp, (str, pd.DataFrame)):
        raise ValueError("main_fp must be one of (str, pd.DataFrame)")
    elif isinstance(main_fp, str):
        reactions = id_indexed(pd.read_csv(os.path.abspath(main_fp), **f_params_in))
    else:
        not_id_indexed = main_fp.index.astype(str).str.isnumeric().all()
        reactions = id_indexed(main_fp)

    data_cols = ['reaction', 'ec', 'smarts', 'sim_max', 'sim_max_id', 'atom_mapping', 'balanced_els_stars']
    preserve = [i for i in data_cols if i in reactions.columns]

    R_meta = reactions.drop(columns=preserve)
    R_meta.to_csv(meta_fp, **f_params_out)
    print(f'Metadata file written to: {meta_fp}')

    reactions.drop(columns=R_meta.columns, inplace=True)
    reactions = reactions[preserve]

    if isinstance(main_fp, pd.DataFrame):
        if not_id_indexed:
            return reactions.reset_index()
        else:
            return reactions
    else:
        reactions.to_csv(main_fp, **f_params_out)
        return None


def list_reactions_per_compound(data_r : pd.DataFrame):
    cps = list_compounds_per_reaction(data_r)
    r2c = cps.explode().reset_index(name=0)
    c2r = r2c.groupby(0)[r2c.columns[0]].apply(set).map(sorted)
    return c2r


def list_compounds_per_reaction(data_r : pd.DataFrame):
    cps = id_indexed(data_r)['reaction'].str.findall(r'(C\d+)')
    r2c = cps.map(set).map(sorted)
    return r2c


def join_col_to_csv(col : pd.Series,
                    file : str,
                    col_name: str|None = None,
                    write_params=None):
    
    df = pd.read_csv(file, index_col=0, low_memory=False)

    icol = col.copy(deep=True)
    if isinstance(col_name, str):
        icol.rename(col_name, inplace=True)

    df.drop(columns=icol.name, inplace=True, errors='ignore')

    if df.index.intersection(icol.index).any():
        df = df.join(icol, how='left')
    elif len(df.index) == len(icol.index):
        df.loc[:,icol.name] = icol.tolist()
    else:
        raise ValueError('column and destination file must share indices or index length')
        return None

    if write_params is None:
        write_params = {'encoding': 'utf-8'}

    df.to_csv(file, index=True, **write_params)
    
    return None


