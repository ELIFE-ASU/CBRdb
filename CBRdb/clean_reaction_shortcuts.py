import os

import pandas as pd
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def scan_for_kegg_pattern(s, prefix):
    """
    Scans a string for KEGG patterns with a specific prefix.

    Parameters:
    s (str): The string to scan for KEGG patterns.
    prefix (str): The prefix to look for in the KEGG patterns.

    Returns:
    list: A list of KEGG patterns that match the specified prefix.
    """
    if prefix in s:
        return [i for i in s.replace(';', ' ').replace('[', ' ').replace(']', ' ').split()
                if i.startswith(prefix) and len(i) == len(prefix) + 5]
    return []


def make_lists_printable(df):
    """
    Converts list-like columns in a DataFrame to printable strings.

    Parameters:
    df (DataFrame): The pandas DataFrame to be processed.

    Returns:
    DataFrame: A copy of the DataFrame with list-like columns converted to strings.
    """
    list_like_cols = df.loc[:, df.iloc[0].apply(lambda x: isinstance(x, (set, list)))]
    printable_df = df.copy(deep=True)
    printable_df[list_like_cols.columns] = list_like_cols.map(' '.join)
    return printable_df


def scan_for_terms(s):
    """
    Scans a string for specific terms defined in the global variable OK.

    Parameters:
    s (str): The string to scan for terms.

    Returns:
    str or None: The original string if any term from OK is found or if the string starts with 'R1' or 'R0', otherwise None.
    """
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return s
    return None


def remove_terms(s):
    """
    Removes specific terms defined in the global variable OK from the string.

    Parameters:
    s (str): The string to remove terms from.

    Returns:
    str or None: None if any term from OK is found or if the string starts with 'R1' or 'R0', otherwise the original string.
    """
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return None
    return s


def flag_missing_reactions(x):
    """
    Identifies missing reactions from a given list.

    Parameters:
    x (list): A list of reaction identifiers to check.

    Returns:
    float or list: NaN if no reactions are missing, otherwise a list of missing reactions.
    """
    missing_reactions = [i for i in x if i not in reactions.index]
    return float('NaN') if not missing_reactions else missing_reactions


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


def find_suspect_reactions(r_file='../data/kegg_data_R.csv.zip', data_dir='../data/'):
    """
    Identifies and flags KEGG reactions that are suspect, including shortcuts
    (summaries of multi-step processes) and reactions with incomplete or general data.

    Parameters:
    r_file (str): The path to the directory containing the preprocessed KEGG reaction csv.
    data_dir (str): The path to the directory where processed data files will be saved.

    Returns:
    None
    """
    reaction_pattern = r'r\d{5}'
    keep_acceptable_nonalnum_chars = lambda x: ''.join([i for i in x if i.isalnum() or i in '-+ ;'])
    str_replacement_order = {'step':'-step', '--':'-', ' -step':' step', 'two':'2', 'three':'3', 'four':'4', 'multi':'0', ' + ':'+',
                            'similar': ';similar', 'incomplete reaction':'', 'unclear reaction':'', 'probably':'possibly', 
                            'possibly':';possibly', ';possibly ;similar': ';possibly similar'}

    def _replace_strings(s:str, replacements:dict):
        for k,v in replacements.items():
            s = s.replace(k,v)
        return s

    def _strip_col1_from_col2(df:pd.DataFrame, col1:str, col2:str):
        for i in df.index:
            df.at[i,col2] = df.at[i,col2].replace(df.at[i,col1], '').strip()
        return df
    
    global reactions, OK
    OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']
    # Prepare the full path of the files
    r_file = os.path.abspath(r_file)
    reactions = pd.read_csv(r_file, header=0)
    get_multistep_details = False

    # Trickiest "shortcuts" to identify: multi-step reactions, compressing other reaction IDs into a net reaction.
    rmulti = ((reactions.dropna(subset='comment')[['id','comment']].copy(deep=True) # info is contained in unstructured comment field
        .apply(lambda x: x.str.lower().str.replace('eaction;', 'eaction, ')) # compress two-line reaction attributes into one line
        .map(keep_acceptable_nonalnum_chars, na_action='ignore') # remove uninformative formatting
        .set_index('id')['comment'].str.split(';').explode().reset_index() # one reaction attribute per line --> split by lines
        .query('comment.str.contains(@reaction_pattern) & comment.str.lower().str.contains("step")') # keep only multi-step parameters
        .sort_values(by='id').set_index('id'))['comment'].apply( # sort by reaction ID
            lambda x: _replace_strings(x, str_replacement_order)) # standardize formatting
            .explode().str.split('reaction', expand=True).reset_index() # general format per line: N-step reaction, see: RXXXXX+RXXXXX
            .rename(columns={'id': 'id', 0: 'n_steps', 1: 'parts'}) # number of steps, parts themselves
            .query('~n_steps.str.contains("of|possibly|one")')) # remove entries that are "part of" a multi-step reaction or "possibly" multi-step
    rmulti = _strip_col1_from_col2(rmulti, 'id', 'parts').apply(lambda x: x.str.strip()) # ensure that parts column contains only constituent steps
    rmulti['parts'] = rmulti['parts'].str.replace(' or ',';').str.split(';') # split step batches into separate lines
    rmulti = (rmulti.explode('parts').query("parts.str.count(@reaction_pattern)>1 & ~parts.str.contains('possibly|similar')") # remove lines with comparisons not step lists
            .sort_values(by='n_steps', ascending=False)).reset_index(drop=True).map(lambda x: x.upper()).query('id!="R10693"') # this is a comparison instance
    rmulti['parts'] = rmulti['parts'].str.replace('R08637','R11101+R11098') # this is a multistep reaction itself

    if get_multistep_details: # optional: extract details of multi-step reactions
        rmulti['parts'] = rmulti['parts'].apply(lambda x: ' '.join([i.lstrip('+') for i in x.split() if 'R1' in i or 'R0' in i]))
        rmulti.update(rmulti.query('~parts.str.contains("+", regex=False)').assign(parts=lambda x: x['parts'].str.split().apply('+'.join)))
        rmulti['parts'] = rmulti['parts'].str.split()
        rmulti = rmulti.explode('parts').reset_index(drop=True)
        rmulti['n_step_sets'] = rmulti['id'].map(rmulti['id'].value_counts())
        rmulti = (rmulti.sort_values(by=['n_step_sets','id'], ascending=[False,True])
                .reset_index(drop=True).drop('n_step_sets', axis=1)
                .apply(lambda x: x.str.rstrip('-STEP')).replace('0','N'))
        rmulti.to_csv(data_dir+'multi_step_reactions.csv.zip', index=False, compression='zip')
    
    reactions_multistep = list(rmulti['id'].unique())
    print('Reactions that are multi-step:', len(reactions_multistep), flush=True)

    # Easiest "shortcuts" to identify: "overall" (i.e. top-level) reactions in
    # [this table](https://www.kegg.jp/kegg/tables/br08210.html), tagged in "overall" field of kegg_data_R.
    reactions_overall = reactions.query('overall.notna()')['id'].tolist()
    print('Reactions that are overall:', len(reactions_overall), flush=True)

    reactions = reactions.set_index('id')
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
    data = data.groupby(by='id')['reason'].apply(lambda x: '+'.join(set(sorted(list(x)))))
    data = data.reset_index().drop_duplicates().sort_values(by='id').set_index('id')
    data.to_csv(os.path.join(data_dir, 'R_IDs_bad.dat'))
    return data

def remove_suspect_reactions(r_file='../data/kegg_data_R.csv.zip', data_dir='../data/'):
    """
    Removes suspect reactions from a reaction data file.
    NOTE: Does NOT remove reactions whose definitions (but not ID) matches a suspect reaction.
    """
    sus = find_suspect_reactions(r_file, data_dir)
    rns = pd.read_csv(r_file, index_col=0)
    if 'kegg_id' in rns.columns:
        rns = rns.query('kegg_id.isin(@sus.index)==False')
    else:
        rns = rns.drop(sus.index, axis=0, errors='ignore')
    rns.to_csv(r_file, compression='zip', encoding='utf-8')
    return rns.reset_index()
