import pandas as pd

from .merge_data_sets import identify_duplicate_compounds, merge_duplicate_compounds, merge_hpc_calculations, id_indexed
from .tools_eq import standardise_eq, get_eq_all_cids, ordered_reaction_series
from .tools_files import reaction_csv, compound_csv


def all_entries(dbs):
    """
    Retrieves all reaction entries and compound IDs from the given databases.

    Parameters:
    dbs (dict): A dictionary containing the following keys:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
                - 'atlas_data_R': A pandas DataFrame with Atlas reaction data.
                - 'CBRdb_C': A pandas DataFrame with compound data.

    Returns:
    tuple: A tuple containing:
           - all_rns (pd.Series): A Series where each element is a set of compound IDs found in the reactions.
           - all_cps (set): A set of all compound IDs for which structures are available.
    """
    kegg_rns = dbs['kegg_data_R'].set_index('id')['reaction']
    atlas_rns = dbs['atlas_data_R'].set_index('id')['reaction']
    all_rns = pd.concat([kegg_rns, atlas_rns]).str.findall(r"(C\d{5})").map(set)
    all_cps = set(dbs['CBRdb_C'].dropna(subset='smiles')['compound_id'])
    return all_rns, all_cps


def all_kegg_comments(dbs):
    """
    Retrieves and processes KEGG comments from the given databases.

    Parameters:
    dbs (dict): A dictionary containing the following key:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.

    Returns:
    pd.DataFrame: A DataFrame with processed KEGG comments and reaction references.
    """
    kegg_cmts = dbs['kegg_data_R'].set_index('id')['comment'].str.replace("reaction;see", "reaction (see")
    kegg_cmts = kegg_cmts.str.split(';').explode().dropna().to_frame()
    kegg_cmts['rn_refs'] = kegg_cmts['comment'].str.upper().str.findall(r"(R\d{5})").map(set).sub(
        kegg_cmts.index.map(lambda x: {x}))
    kegg_cmts.at['R10693', 'comment'] = 'part of ' + kegg_cmts.at['R10693', 'comment']
    return kegg_cmts


def list_reactions_with_halogen_dupes(dbs):
    """
    Identifies reactions that contain duplicate halogen compounds.

    Parameters:
    dbs (dict): A dictionary containing the following keys:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
                - 'atlas_data_R': A pandas DataFrame with Atlas reaction data.
                - 'CBRdb_C': A pandas DataFrame with compound data.

    Returns:
    pd.Index: An index of reaction IDs that contain duplicate halogen compounds.
    """
    all_rns, all_cps = all_entries(dbs)
    cps_duped = dbs['CBRdb_C'].query('compound_id.str.contains("C99") & smiles.duplicated(keep=False)')[
        'compound_id'].values
    rns = all_rns[~all_rns.map(lambda x: x.isdisjoint(cps_duped))].index.drop_duplicates()
    return rns


def list_reactions_missing_structures(dbs):
    """
    Identifies reactions that are missing structures.

    Parameters:
    dbs (dict): A dictionary containing the following keys:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
                - 'atlas_data_R': A pandas DataFrame with Atlas reaction data.
                - 'CBRdb_C': A pandas DataFrame with compound data.

    Returns:
    pd.Index: An index of reaction IDs that are missing structures.
    """
    all_rns, all_cps = all_entries(dbs)
    rns = all_rns[~all_rns.map(lambda x: x.issubset(all_cps))].index.drop_duplicates()
    return rns


def list_reactions_with(dbs, phrase):
    """
    Identifies reactions that contain a specific phrase in their comments.

    Parameters:
    dbs (dict): A dictionary containing the following key:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
    phrase (str): The phrase to search for in the comments.

    Returns:
    pd.Index: An index of reaction IDs that contain the specified phrase in their comments.
    """
    kegg_cmts = all_kegg_comments(dbs)
    query = f'comment.str.lower().str.contains("{phrase}")'
    rns = kegg_cmts.query(query).index.drop_duplicates()
    return rns


def list_multistep_parts(dbs):
    """
    Identifies reactions that are part of multistep processes.

    Parameters:
    dbs (dict): A dictionary containing the following key:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.

    Returns:
    pd.Index: An index of reaction IDs that are part of multistep processes.
    """
    cmts = all_kegg_comments(dbs)['comment']
    rns = cmts[cmts.str.contains(" of ") & cmts.str.contains("step")].index.drop_duplicates()
    return rns


def list_multistep_enumerated(dbs, verbose=False):
    """
    Identifies reactions that are part of multistep processes and enumerates them.

    Parameters:
    dbs (dict): A dictionary containing the following key:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
    verbose (bool): If True, saves a file enumerating each reaction's constituent steps. Default is False.
    Returns:
    pd.Index: An index of reaction IDs that are part of multistep processes.
    """
    if dbs['kegg_data_R'].index.name == 'id':
        dbs['kegg_data_R'].reset_index(inplace=True)
    reactions_multistep_parts = list_multistep_parts(dbs)  # KEEP THESE 
    rns = all_kegg_comments(dbs).drop(reactions_multistep_parts).query(
        'rn_refs.str.len()>1 & comment.str.contains("step") & ~comment.str.contains("possibl|probabl|similar|Overall Reaction")')
    if verbose:
        multistep_enum = pd.concat([dbs['kegg_data_R'].set_index('id').loc[['R10671'], ['comment']],
                                    rns])
        multistep_enum['rn_refs'] = multistep_enum['rn_refs'].map(lambda x: ' '.join(list(x)), na_action='ignore')
        multistep_enum.to_csv('../data/multistep_enumerated.csv', encoding='utf-8')
    rns = rns.index.drop_duplicates().union({'R10671'})  # false negative: "similar" in line
    return rns


def df_of_suspect_reactions(dbs):
    """
    Identifies and categorizes suspect reactions from the given databases.

    Parameters:
    dbs (dict): A dictionary containing the following keys:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
                - 'atlas_data_R': A pandas DataFrame with Atlas reaction data.
                - 'CBRdb_C': A pandas DataFrame with compound data.

    Returns:
    pd.DataFrame: A DataFrame with reaction IDs as the index and a column 'reason'
                  indicating the reasons why each reaction is considered suspect.
    """
    reactions_general = list_reactions_with(dbs, 'general reaction')
    reactions_incomplete = list_reactions_with(dbs, 'incomplete|incomplete reaction')
    reactions_missing_structures = list_reactions_missing_structures(dbs)
    reactions_multistep = list_multistep_enumerated(dbs)
    reactions_unclear = list_reactions_with(dbs, 'unclear reaction')

    sus = {'general': reactions_general,
           'incomplete': reactions_incomplete,
           'unclear': reactions_unclear,
           'structure_missing': reactions_missing_structures,
           'shortcut': reactions_multistep}

    sus = pd.concat([pd.Series(index=v, data=len(v) * [k]) for k, v in sus.items()])
    sus = sus.groupby(level=0).apply(list).map(lambda x: '+'.join(sorted(x))).to_frame('reason')

    return sus


def suspect_reaction_subset(sus, matching):
    """
    Filters and returns a subset of suspect reactions that match a given pattern.

    Parameters:
    sus (pd.DataFrame): A DataFrame with reaction IDs as the index and a column 'reason'
                        indicating the reasons why each reaction is considered suspect.
    matching (str): A regex pattern to match against the 'reason' column.

    Returns:
    pd.Index: An index of reaction IDs that match the given pattern in their 'reason' column.
    """
    return sus.map(lambda x: x.split('+')).explode('reason').query(
        "reason.str.contains(@matching)").index.drop_duplicates()


def add_suspect_reactions_to_existing_bad_file(sus_new, R_IDs_bad_file="../data/R_IDs_bad.dat"):
    """
    Adds new suspect reactions to an existing file of bad reactions.

    Parameters:
    sus_new (pd.DataFrame): A DataFrame with reaction IDs as the index and a column 'reason'
                            indicating the reasons why each reaction is considered suspect.
    R_IDs_bad_file (str): The file path to the existing bad reactions file.

    Returns:
    pd.DataFrame: A DataFrame with combined suspect reactions from the new and existing files.
    """
    sus_old = pd.read_csv(R_IDs_bad_file, header=0, index_col=0)['reason'].str.split('+').explode()
    sus_new = sus_new['reason'].str.split('+').explode()
    sus_combo = (pd.concat([sus_old, sus_new]).groupby(level=0)
                 .apply(lambda x: sorted(list(set(x)))).map('+'.join).to_frame('reason')).sort_index()
    sus_combo.to_csv(R_IDs_bad_file)
    return sus_combo


def quarantine_suspect_reactions_matching(dbs, sus, matching="shortcut|structure_missing"):
    """
    Quarantines suspect reactions that match a given pattern from the databases.

    Parameters:
    dbs (dict): A dictionary containing the following keys:
                - 'kegg_data_R': A pandas DataFrame with KEGG reaction data.
                - 'atlas_data_R': A pandas DataFrame with Atlas reaction data.
    sus (pd.DataFrame): A DataFrame with reaction IDs as the index and a column 'reason'
                        indicating the reasons why each reaction is considered suspect.
    matching (str): A regex pattern to match against the 'reason' column. Default is "shortcut|structure_missing".

    Returns:
    dict: The updated databases dictionary with quarantined reactions.
    """
    to_quarantine = suspect_reaction_subset(sus, matching)
    dbs['sus'] = sus
    dbs['quarantine_categories'] = matching
    dbs['quarantined'] = pd.concat([dbs['kegg_data_R'], dbs['atlas_data_R']]).query('id.isin(@to_quarantine)').copy(
        deep=True)
    dbs['kegg_data_R_orig'] = dbs['kegg_data_R'].copy(deep=True)
    dbs['atlas_data_R_orig'] = dbs['atlas_data_R'].copy(deep=True)
    dbs['kegg_data_R'] = dbs['kegg_data_R'].query('~id.isin(@to_quarantine)')
    dbs['atlas_data_R'] = dbs['atlas_data_R'].query('~id.isin(@to_quarantine)')
    return dbs


def iteratively_prune_entries(kegg_data_R, atlas_data_R, C_main):
    """
    Prunes reactions before attempting to dedupe compounds.

    Parameters:
    kegg_data_R (pd.DataFrame): A pandas DataFrame with KEGG reaction data.
    atlas_data_R (pd.DataFrame): A pandas DataFrame with Atlas reaction data.
    C_main (pd.DataFrame): A pandas DataFrame with main compound data.

    Returns:
    dict: The updated databases dictionary with pruned and deduplicated compound entries.
    """
    # turn datasets into a dictionary
    dbs = {'kegg_data_R': kegg_data_R, 'atlas_data_R': atlas_data_R, 'kegg_data_C': C_main}

    # identify all duplicated compounds
    dbs['CBRdb_C'] = dbs['kegg_data_C'].copy(deep=True)
    dbs['C_dupemap'] = identify_duplicate_compounds(dbs['CBRdb_C'])
    dbs['CBRdb_C'] = merge_duplicate_compounds(dbs['CBRdb_C'], dbs['C_dupemap'])

    for db in ['kegg_data_R', 'atlas_data_R']:
        # replace duped compound IDs in reaction DataFrames
        dbs[db].loc[:, 'reaction'] = dbs[db].loc[:, 'reaction'].str.split(expand=True).replace(
            dbs['C_dupemap']['new_id']).fillna('').apply(lambda x: ' '.join(x), axis=1).str.strip()
        # standardize reaction format again
        dbs[db]['reaction'] = dbs[db]['reaction'].map(standardise_eq)
        dbs[db]['CBRdb_C_ids'] = dbs[db]['reaction'].map(get_eq_all_cids)

    # identify suspect reactions and those matching them after compound de-duping
    sus = df_of_suspect_reactions(dbs)
    sus = add_sus_reaction_dupes(sus, dbs)
    sus = add_suspect_reactions_to_existing_bad_file(sus)

    # write CSV output files for the reaction balancer to read in.
    for k in ['kegg_data_R', 'atlas_data_R']:
        reaction_csv(dbs[k], f'../data/{k}_dedupedCs.csv')

    # Separate compound data from metadata
    kegg_meta, xref_meta = dbs['CBRdb_C'].filter(like='kegg_'), dbs['CBRdb_C'].filter(like='xref_')
    dbs['CBRdb_C_metadata'] = dbs['CBRdb_C'][['comment']].join(kegg_meta).join(xref_meta).copy(deep=True)
    dbs['CBRdb_C'].drop(columns=dbs['CBRdb_C_metadata'].columns, inplace=True)
    dbs['CBRdb_C_metadata'].loc[:,'compound_id'] = dbs['CBRdb_C'].loc[:,'compound_id']

    # Merge in compound property calculations from HPC
    dbs['CBRdb_C'] = merge_hpc_calculations(id_indexed(dbs['CBRdb_C'])).reset_index()

    # Write preliminary compound data and metadata files
    compound_csv(dbs['CBRdb_C'], '../CBRdb_C.csv.zip')
    compound_csv(dbs['CBRdb_C_metadata'], '../CBRdb_C_metadata.csv.zip')

    return dbs


def add_sus_reaction_dupes(sus, dbs):
    """
    Identifies and adds duplicate reactions of suspect reactions to the suspect list.

    This function finds reactions in the KEGG and Atlas databases that have the
    exact same chemical equation as reactions already marked as suspect. The
    comparison ignores reaction directionality by sorting the reactant and product
    sides of the equation string. These newly found duplicates are then added to
    the suspect DataFrame.

    Parameters
    ----------
    sus : pd.DataFrame
        A DataFrame of initial suspect reactions, indexed by reaction ID, with a
        'reason' column.
    dbs : dict
        A dictionary containing the primary DataFrames, including 'kegg_data_R'
        and 'atlas_data_R'.

    Returns
    -------
    pd.DataFrame
        The updated suspect DataFrame, now including reactions that are
        duplicates of the original suspect reactions.
    """
    sort_sides = lambda df: df['reaction'].str.split(' <=> ').map(lambda x: ' <=> '.join(sorted(x)))

    # Add equations for suspect Atlas reactions (above func only checks for missing structures))
    sus = sus.join(dbs['atlas_data_R'].set_index('id')[['reaction']])
    # Add equations for suspect KEGG reactions (also checks for multi-step, incomplete, etc.)
    sus.update(dbs['kegg_data_R'].set_index('id')[['reaction']])

    # Sort sides to remove directionality-only differences
    sus['reaction'] = sort_sides(sus)
    # Do the same for Atlas and KEGG reactions
    atlas_eqs = sort_sides(dbs['atlas_data_R'].set_index('id'))
    kegg_eqs = sort_sides(dbs['kegg_data_R'].set_index('id'))
    # Consider the not-yet-suspect reactions
    pre_sus_reactions = pd.concat([atlas_eqs, kegg_eqs]).drop(sus.index)
    # For each equation, find matching not-yet-suspect reactions
    pre_sus_reactions = pre_sus_reactions.reset_index().groupby('reaction')['id'].apply(list)
    # If a suspect reaction's equation matches a not-yet-suspect reaction's equation, add it to the suspect list
    s = 'matches'
    sus[s] = sus.reaction.map(pre_sus_reactions)
    sus = pd.concat([sus, sus.dropna(subset=s).explode(s).set_index(s).rename_axis('id')])
    sus = sus.drop(columns=['reaction', s])
    return sus



def generate_reaction_dupemap(df, prefix="T"):
    """
    Transforms the 'reaction' column of a DataFrame to identify duplicate reactions and returns a map of oldID to newID with a user-specified prefix.

    Parameters:
    df (pd.DataFrame): The DataFrame containing reaction data with an 'id' column and a 'reaction' column.
    prefix (str): The prefix to use for the new IDs. Defaults to "T".

    Returns:
    pd.Series: A Series mapping old IDs to new IDs with the specified prefix.
    """
    # Group the reactions by sorting and standardizing them
    equation_groups = ordered_reaction_series(df.set_index('id')['reaction']).to_frame(name='eq_grp')

    # Identify duplicated groups and assign a group number
    duped_groups = equation_groups.query('eq_grp.duplicated(keep=False)').assign(
        grp_num=lambda x: x.eq_grp.factorize()[0])

    # Create a map of old IDs to new IDs with the specified prefix
    dupemap = prefix + duped_groups['grp_num'].astype(str).str.zfill(5)

    # Save the map to a CSV file
    dupemap.sort_index().to_csv(f'../data/{prefix}_reaction_dupemap.csv', encoding='utf-8')

    return dupemap



def sync_reaction_dupemap(df, prefix="T"):

    dupemap_filename = f'../data/{prefix}_reaction_dupemap.csv'
    try:
        existing_dupemap = pd.read_csv(dupemap_filename, index_col=0)
    except:
        existing_dupemap = generate_reaction_dupemap(df, prefix=prefix)
        return existing_dupemap
    
    # Standardize the formatting of current reaction equations
    grps = ordered_reaction_series(id_indexed(df)['reaction'])

    # Identify the reactions that are currently duplicated
    grps = grps[grps.duplicated(keep=False)].to_frame(name='reaction')

    # At least some reactions have previously been assigned a dupe-group ID
    grps = grps.join(existing_dupemap)

    # Flag newly-identified dupes (i.e. those lacking a pre-existing dupe-group ID)
    dupe_is_new = grps['grp_num'].isna()

    # Flag members of internally inconsistent dupe-groups (where group and eqn are not a 1:1 mapping)
    eqs_per_grp = grps.groupby('grp_num')['reaction'].nunique()
    grps_per_eq = grps.groupby('reaction')['grp_num'].nunique()
    gt1_eq_per_grp = grps['grp_num'].map(eqs_per_grp).gt(1)
    gt1_grp_per_eq = grps['reaction'].map(grps_per_eq).gt(1)

    # Dupe-groups with a 1:1 mapping should be assigned the same dupe-group ID as before
    dupes_known = ~(gt1_eq_per_grp | gt1_grp_per_eq | dupe_is_new)

    # Use this information to start making the new mapping scheme
    dupemap = grps.loc[dupes_known].copy(deep=True)
    eqn_to_grp = dict(zip(dupemap['reaction'], dupemap['grp_num']))
    
    # Check the unknown subset to see if any of their reactions correspond with a known dupe-group
    to_assign = grps.loc[~dupes_known].copy(deep=True)
    to_assign.update({'grp_num': to_assign['reaction'].map(eqn_to_grp)}, overwrite=True)

    # Define each mixed dupe-group ID by its most common eqn, ensuring each eqn only points to one group
    largest_mixed_groups = (to_assign.groupby(by='grp_num', as_index=False)['reaction']
                            .value_counts().sort_values(by='count', ascending=False)
                            .drop_duplicates(subset='reaction')
                            .set_index('reaction')['grp_num'])
    to_assign.update({'grp_num': to_assign['reaction'].map(largest_mixed_groups)}, overwrite=True)
    # Everything that can be assigned to an existing group number now has been. So update known entries 
    dupemap = pd.concat([to_assign.dropna(subset='grp_num'), dupemap], axis=0)
    # Now focus on remaining unknowns (i.e. equations without an existing group)
    to_assign = to_assign.query('grp_num.isna()')
    # Find the largest assigned dupe-group ID number
    largest_id = dupemap['grp_num'].str.strip(prefix).astype(int).max()
    # Start one up from there
    to_assign.loc[:,'grp_num'] = to_assign['reaction'].factorize()[0] + (largest_id + 1)
    # Format by the same scheme as before
    zfill_len = list(set(dupemap['grp_num'].map(len)))[0] - len(prefix)
    to_assign.loc[:,'grp_num'] = prefix + to_assign['grp_num'].astype(str).str.zfill(zfill_len)
    # Update dupemap
    dupemap = pd.concat([dupemap, to_assign], axis=0)
    # Sort again by ID 
    dupemap = dupemap.sort_index()['grp_num']
    # Save as a CSV
    dupemap.to_csv(dupemap_filename, encoding='utf-8')
    return dupemap

