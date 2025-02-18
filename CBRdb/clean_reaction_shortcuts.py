import os

import pandas as pd
from rdkit import RDLogger

from .tools_files import add_suffix_to_file

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def load_reactions_data(target_dir):
    """
    Loads reaction data from a specified directory.

    Parameters:
    target_dir (str): The directory containing the reaction data files.

    Returns:
    DataFrame: A pandas DataFrame containing the loaded reaction data.
    """
    a = os.listdir(target_dir)
    reactions = {}
    for i in a:
        file = os.path.abspath(f"{target_dir}/{i}/{i}.data")
        try:
            with open(file, 'r') as f:
                entries = {}
                for line in f.readlines()[:-1]:
                    if not line.startswith(' '):
                        field = line.split(' ')[0]
                        entries[field] = line.lstrip(field).lstrip().rstrip('\n')
                    else:
                        entries[field] += '; ' + line.lstrip().rstrip('\n')
                reactions[i] = entries
        except:
            print(f"Error reading {i}.data", flush=True)
            raise FileExistsError
    return pd.DataFrame(reactions).fillna('').T


def get_reaction_ids_substr(reactions, substr="incomplete reaction"):
    """
    Retrieves reaction IDs from a DataFrame where the 'COMMENT' column contains a specific substring.

    Parameters:
    reactions (DataFrame): The pandas DataFrame containing reaction data.
    substr (str): The substring to search for within the 'COMMENT' column. Default is "incomplete reaction".

    Returns:
    list: A list of reaction IDs where the 'COMMENT' column contains the specified substring.
    """
    incomplete_reaction_ids = reactions[
        reactions['comment'].str.contains(substr, case=False, na=False)].index.tolist()
    return incomplete_reaction_ids


def find_suspect_reactions(r_file='../data/kegg_data_R.csv', data_dir='../data/', verbose=False):
    """
    Identifies and flags KEGG reactions that are suspect, including shortcuts
    (summaries of multi-step processes) and reactions with incomplete or general data.

    Parameters:
    r_file (str): The path to the directory containing the preprocessed KEGG reaction csv.
    data_dir (str): The path to the directory where processed data files will be saved.
    verbose (bool, default False): Whether to create files preserving info on suspect reactions found

    Returns:
    pd.Series: A Series with index=ID and value=reason(s) reaction was flagged as suspect
    """
    reaction_pattern = r'r\d{5}'
    keep_acceptable_nonalnum_chars = lambda x: ''.join([i for i in x if i.isalnum() or i in '-+ ;'])
    str_replacement_order = {'step': '-step', '--': '-', ' -step': ' step', 'two': '2', 'three': '3', 'four': '4',
                             'multi': '0', ' + ': '+',
                             'similar': ';similar', 'incomplete reaction': '', 'unclear reaction': '',
                             'probably': 'possibly',
                             'possibly': ';possibly', ';possibly ;similar': ';possibly similar'}

    def _replace_strings(s: str, replacements: dict):
        for k, v in replacements.items():
            s = s.replace(k, v)
        return s

    def _strip_col1_from_col2(df: pd.DataFrame, col1: str, col2: str):
        for i in df.index:
            df.at[i, col2] = df.at[i, col2].replace(df.at[i, col1], '').strip()
        return df

    global reactions, OK
    OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']
    # Prepare the full path of the files
    r_file = os.path.abspath(r_file)
    reactions = pd.read_csv(r_file, header=0)

    # Trickiest "shortcuts" to identify: multi-step reactions, compressing other reaction IDs into a net reaction.
    rmulti = ((reactions.dropna(subset='comment')[['id', 'comment']].copy(
        deep=True)  # info is contained in unstructured comment field
               .apply(lambda x: x.str.lower().str.replace('eaction;',
                                                          'eaction, '))  # compress two-line reaction attributes into one line
               .map(keep_acceptable_nonalnum_chars, na_action='ignore')  # remove uninformative formatting
               .set_index('id')['comment'].str.split(
        ';').explode().reset_index()  # one reaction attribute per line --> split by lines
               .query(
        'comment.str.contains(@reaction_pattern) & comment.str.lower().str.contains("step")')  # keep only multi-step parameters
               .sort_values(by='id').set_index('id'))['comment'].apply(  # sort by reaction ID
        lambda x: _replace_strings(x, str_replacement_order))  # standardize formatting
    .explode().str.split('reaction',
                         expand=True).reset_index()  # general format per line: N-step reaction, see: RXXXXX+RXXXXX
    .rename(columns={'id': 'id', 0: 'n_steps', 1: 'parts'})  # number of steps, parts themselves
    .query(
        '~n_steps.str.contains("of|possibly|one")'))  # remove entries that are "part of" a multi-step reaction or "possibly" multi-step
    rmulti = _strip_col1_from_col2(rmulti, 'id', 'parts').apply(
        lambda x: x.str.strip())  # ensure that parts column contains only constituent steps
    rmulti['parts'] = rmulti['parts'].str.replace(' or ', ';').str.split(';')  # split step batches into separate lines
    rmulti = (rmulti.explode('parts').query(
        "parts.str.count(@reaction_pattern)>1 & ~parts.str.contains('possibly|similar')")  # remove lines with comparisons not step lists
              .sort_values(by='n_steps', ascending=False)).reset_index(drop=True).map(lambda x: x.upper()).query(
        'id!="R10693"')  # this is a comparison instance
    rmulti['parts'] = rmulti['parts'].str.replace('R08637', 'R11101+R11098')  # this is a multistep reaction itself

    if verbose:  # extract stepwise details of multi-step reactions
        rmulti['parts'] = rmulti['parts'].apply(
            lambda x: ' '.join([i.lstrip('+') for i in x.split() if 'R1' in i or 'R0' in i]))
        rmulti.update(rmulti.query('~parts.str.contains("+", regex=False)').assign(
            parts=lambda x: x['parts'].str.split().apply('+'.join)))
        rmulti['parts'] = rmulti['parts'].str.split()
        rmulti = rmulti.explode('parts').reset_index(drop=True)
        rmulti['n_step_sets'] = rmulti['id'].map(rmulti['id'].value_counts())
        rmulti = (rmulti.sort_values(by=['n_step_sets', 'id'], ascending=[False, True])
                  .reset_index(drop=True).drop('n_step_sets', axis=1)
                  .apply(lambda x: x.str.rstrip('-STEP')).replace('0', 'N'))
        rmulti.to_csv(add_suffix_to_file(r_file, 'suspect_steps'), index=False)

    reactions_multistep = list(rmulti['id'].unique())
    print('Reactions that are multi-step:', len(reactions_multistep), flush=True)

    # Easiest "shortcuts" to identify: "overall" (i.e. top-level) reactions in
    # [this table](https://www.kegg.jp/kegg/tables/br08210.html), tagged in "overall" field of kegg_data_R.
    reactions_overall = reactions.query('overall.notna()')['id'].tolist()
    print('Reactions that are overall:', len(reactions_overall), flush=True)

    reactions = reactions.set_index('id').replace('incomplete', 'incomplete reaction')
    reactions_incomplete = get_reaction_ids_substr(reactions, substr="incomplete reaction")
    print('Reactions with incomplete data:', len(reactions_incomplete), flush=True)

    reactions_general = get_reaction_ids_substr(reactions, substr="general reaction")
    print('General reactions:', len(reactions_general), flush=True)

    data = pd.DataFrame({
        'id': reactions_multistep + reactions_overall + reactions_incomplete + reactions_general,
        'reason': ['shortcut'] * len(reactions_multistep)
                  + ['shortcut'] * len(reactions_overall)
                  + ['incomplete'] * len(reactions_incomplete)
                  + ['general'] * len(reactions_general)}).sort_values(by=['id']).reset_index(drop=True)
    # open existing R_IDs_bad.dat file and append new data
    if os.path.exists(os.path.join(data_dir, 'R_IDs_bad.dat')):
        data_old = pd.read_csv(os.path.join(data_dir, 'R_IDs_bad.dat'), header=0)
        data_old['reason'] = data_old['reason'].str.split('+')
        data_old = data_old.explode('reason')
        data = pd.concat([data_old, data], ignore_index=True)
    data = data.groupby(by='id')['reason'].apply(lambda x: '+'.join(sorted(list(set(x)))))
    data = data.reset_index().drop_duplicates().sort_values(by='id').set_index('id')
    data.to_csv(os.path.join(data_dir, 'R_IDs_bad.dat'))

    if verbose:
        reactions.query('compound_id.isin(@data.id)').to_csv(add_suffix_to_file(r_file, 'suspect'), index=False)

    return data


def remove_suspect_reactions(r_file='../data/kegg_data_R.csv', data_dir='../data/', verbose=False):
    """
    Removes suspect reactions from a reaction data file.
    NOTE: Does NOT remove reactions whose definitions (but not ID) matches a suspect reaction.
    """
    sus = find_suspect_reactions(r_file, data_dir, verbose)
    rns = pd.read_csv(r_file, index_col=0)
    if 'kegg_id' in rns.columns:
        rns = rns.query('kegg_id.isin(@sus.index)==False')
    else:
        rns = rns.drop(sus.index, axis=0, errors='ignore').sort_index()
    rns.to_csv(r_file, encoding='utf-8')
    return rns.reset_index()
