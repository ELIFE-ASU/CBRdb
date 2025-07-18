import pandas as pd

from .merge_data_sets import identify_duplicate_compounds
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
           - all_cps (set): A set of all compound IDs.
    """
    kegg_rns = dbs['kegg_data_R'].set_index('id')['reaction']
    atlas_rns = dbs['atlas_data_R'].set_index('id')['reaction']
    all_rns = pd.concat([kegg_rns, atlas_rns]).str.findall(r"(C\d{5})").map(set)
    all_cps = set(dbs['CBRdb_C']['compound_id'])
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
    reactions_overall = dbs['kegg_data_R'].set_index('id')['overall'].dropna().index
    reactions_multistep_parts = list_multistep_parts(dbs)  # KEEP THESE 
    rns = all_kegg_comments(dbs).drop(reactions_multistep_parts).query(
        'rn_refs.str.len()>1 & comment.str.contains("step") & ~comment.str.contains("possibl|probabl|similar")')
    if verbose:
        multistep_enum = pd.concat([dbs['kegg_data_R'].set_index('id').loc[['R10671'], ['comment']],
                                     dbs['kegg_data_R'].set_index('id').loc[reactions_overall, ['comment']],
                                     rns])
        multistep_enum['rn_refs'] = multistep_enum['rn_refs'].map(lambda x: ' '.join(list(x)), na_action='ignore')
        multistep_enum.to_csv('../data/multistep_enumerated.csv', encoding='utf-8')
    rns = rns.index.union(reactions_overall).drop_duplicates().union({'R10671'})  # false negative: "similar" in line
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


def iteratively_prune_entries(kegg_data_R, atlas_data_R, C_main, to_quarantine="shortcut|structure_missing"):
    """
    Prunes reactions before attempting to dedupe compounds.

    Parameters:
    kegg_data_R (pd.DataFrame): A pandas DataFrame with KEGG reaction data.
    atlas_data_R (pd.DataFrame): A pandas DataFrame with Atlas reaction data.
    C_main (pd.DataFrame): A pandas DataFrame with main compound data.
    to_quarantine (str): A regex pattern to match against the 'reason' column for quarantining reactions. Default is "shortcut|structure_missing".

    Returns:
    dict: The updated databases dictionary with pruned and deduplicated entries.
    """
    # turn datasets into a dictionary
    dbs = {'kegg_data_R': kegg_data_R, 'atlas_data_R': atlas_data_R, 'kegg_data_C': C_main}
    dbs['CBRdb_C'] = dbs['kegg_data_C'].copy(deep=True)

    sus = df_of_suspect_reactions(dbs)  # identify suspect reactions
    sus = add_sus_reaction_dupes(sus, dbs)  # add reactions whose equations match that of a suspect reaction
    sus = add_suspect_reactions_to_existing_bad_file(sus)  # add to log and import log

    # remove reactions matching quarantine specs
    dbs = quarantine_suspect_reactions_matching(dbs, sus, matching=to_quarantine)

    # note what compounds are actually used in remaining reactions
    all_rns, all_cps = all_entries(dbs)

    # a few compounds were only used in those reactions; remove those. kegg_data_C retains quarantined entries
    all_rns, all_cps = all_entries(dbs)
    dbs['CBRdb_C'] = dbs['CBRdb_C'].query('compound_id.isin(@all_cps)')

    # now identify and log duplicate compounds among what remains
    dbs['C_dupemap'] = identify_duplicate_compounds(dbs['CBRdb_C'])

    # replace compound entries with dupe names
    dbs['CBRdb_C']['compound_id'] = dbs['CBRdb_C']['compound_id'].replace(dbs['C_dupemap']['new_id'])
    # de-duplicate compounds
    dbs['CBRdb_C'] = dbs['CBRdb_C'].sort_values(by='compound_id').drop_duplicates(subset='compound_id', keep='first')

    # replace compound IDs in reaction dfs. kegg_data_R_orig and atlas_data_R_orig retain original entries.
    dbs['kegg_data_R'].loc[:, 'reaction'] = dbs['kegg_data_R'].loc[:, 'reaction'].str.split(expand=True).replace(
        dbs['C_dupemap']['new_id']).fillna('').apply(lambda x: ' '.join(x), axis=1).str.strip()
    dbs['atlas_data_R'].loc[:, 'reaction'] = dbs['atlas_data_R'].loc[:, 'reaction'].str.split(expand=True).replace(
        dbs['C_dupemap']['new_id']).fillna('').apply(lambda x: ' '.join(x), axis=1).str.strip()

    # write CSV output files for the reaction balancer to read in.
    for k in ['kegg_data_R', 'atlas_data_R']:
        reaction_csv(dbs[k], f'../data/{k}_dedupedCs.csv')
    compound_csv(dbs['CBRdb_C'], '../CBRdb_C.csv')

    return dbs


def add_sus_reaction_dupes(sus, dbs):
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
