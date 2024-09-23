# Input: An unzipped, locally-downloaded folder of raw KEGG reaction data files.
# Output: CSVs with relevant reaction info, e.g. ties to other databases and flags for possible removal.
import os

import pandas as pd


def scan_for_kegg_pattern(s, prefix):
    if prefix in s:
        as_list = s.replace(';', ' ').replace('[', ' ').replace(']', ' ').split(' ')
        return [i for i in as_list if i.startswith(prefix) and len(i) == len(prefix) + 5]
    else:
        return []


def make_lists_printable(df):
    list_like_cols = df.loc[:, df.iloc[0].apply(type).apply(lambda x: x in [set, list])]
    printable_cols = list_like_cols.applymap(' '.join)
    printable_df = df.copy(deep=True)
    printable_df.loc[:, printable_cols.columns] = printable_cols
    return printable_df


def scan_for_terms(s):
    for j in OK:
        if (j in s) or s.startswith('R1') or s.startswith('R0'):
            return s
    else:
        return None


def remove_terms(s):
    for j in OK:
        if (j in s) or s.startswith('R1') or s.startswith('R0'):
            return None
    else:
        return s


# remove reactions for which - for all step-sequences available - not all steps shown are actually present in our dataset.
def flag_missing_reactions(x):
    a = [i for i in x if i not in reactions.index]
    if len(a) == 0:
        return float('Nan')
    else:
        return a


prefix = 'kegg_data_R/'
a = os.listdir(prefix)
reactions = dict()
for i in a:
    try:
        f = open(prefix + i + '/' + i + '.data', 'r')
        entries = dict()
        for line in f.readlines()[:-1]:
            if not line.startswith(' '):
                field = line.split(' ')[0]
                entries[field] = line.lstrip(field).lstrip().rstrip('\n')
            else:
                entries[field] += '; ' + line.lstrip().rstrip('\n')
        reactions[i] = entries
    except:
        pass

reactions = pd.DataFrame(reactions).fillna('').T

# Parse each reaction's ties to other databases.
# Flag reactions with glycan participants and "overall" reactions (see [this table](https://www.kegg.jp/kegg/tables/br08210.html) for info).

reactions = reactions.assign(
    OVERALL_FLAG=reactions['ENTRY'].str.contains('Overall', case=False),
    GLYCAN_FLAG=reactions['EQUATION'].str.contains('G'),
    EQUATION_ENTRIES=reactions['EQUATION'].str.findall(r'C\d{5}|G\d{5}'),
    ENZYME_ENTRIES=reactions['ENZYME'].str.replace(';', ' ').str.split(),
    DBLINKS_RHEA_ENTRIES=reactions['DBLINKS'].str.lstrip('RHEA: ').str.split())
prefixes = {'RCLASS': r'RC\d{5}', 'BRITE': r'br\d{5}', 'PATHWAY': r'rn\d{5}', 'MODULE': r'M\d{5}',
            'ORTHOLOGY': r'K\d{5}'}
for col, prefix in prefixes.items():
    reactions[col + '_ENTRIES'] = reactions[col].str.findall(prefix)

# KEGG no longer maintains the RCLASS database, BUT some edges in KEGG pathways
# are still tagged with RCLASS entries instead of their constituent reactions.
# Let's reconstruct RCLASS from REACTION (since we get that for free).
r_rclass = reactions['RCLASS'].str.split('; ').explode().to_frame()
r_rclass = r_rclass.assign(RC=r_rclass['RCLASS'].str.extract(r'(RC\d{5})'),
                           CPAIR=r_rclass['RCLASS'].str.findall(r'(C\d{5}_C\d{5})')).drop('RCLASS', axis=1).explode(
    'CPAIR')
rclass_db = r_rclass.reset_index(names='REACTION').groupby(by=['RC', 'CPAIR'])['REACTION'].apply(set).apply(
    list).to_frame()
make_lists_printable(rclass_db).to_csv('../RCLASS_DB.csv')

# Use comments to identify multi-step reactions - i.e., artificial shortcuts compressing multiple existing reaction IDs into a single net reaction.

to_replace = {'MULTISTEP': 'MULTI-STEP', '3STEP': '3-STEP', '; SEE R': ' (SEE R', 'REACTION;': 'REACTION '}
dm = reactions['COMMENT'].replace('', float('NaN')).dropna().str.upper().reset_index()
dm['ORIGINAL_COMMENT'] = dm['COMMENT'].copy(deep=True)
for k, v in to_replace.items(): dm['COMMENT'] = dm['COMMENT'].str.replace(k, v)
dm['COMMENT'] = dm['COMMENT'].str.split(';')
dm = dm.explode('COMMENT').reset_index(drop=True).query('COMMENT.str.contains("STEP REACTION") and \
    COMMENT.str.contains("R\d{5}\+|R\d{5} \+|R\d{5} \,") and \
    not COMMENT.str.contains("STEP OF|PART OF|THE LAST|THE FIRST|THE FORMER|THE LATTER|POSSIBLY|PROBABLY IDENTICAL")')
to_replace = {'TWO': '2', 'THREE': '3', 'FOUR': '4', 'MULTI': 'N', ' + ': '+'}
for k, v in to_replace.items():
    dm['COMMENT'] = dm['COMMENT'].str.replace(k, v).str.strip()
dm = dm.sort_values('COMMENT')

format_to_remove = r'SIMILAR TO' + r' R\d{5}, R\d{5}\+R\d{5}'
for nstep in range(2, 15):
    to_remove = dm['COMMENT'].str.findall(format_to_remove + r'\)').explode().dropna()
    for i, str_to_remove in to_remove.items():
        dm.loc[i, 'COMMENT'] = dm.loc[i, 'COMMENT'].replace(str_to_remove, '')
    format_to_remove += r'\+R\d{5}'
format_to_remove = r'SIMILAR TO' + r' R\d{5}\+R\d{5}'
for nstep in range(2, 15):
    to_remove = dm['COMMENT'].str.findall(format_to_remove + r'\)').explode().dropna()
    for i, str_to_remove in to_remove.items():
        dm.loc[i, 'COMMENT'] = dm.loc[i, 'COMMENT'].replace(str_to_remove, '')
    format_to_remove += r'\+R\d{5}'
reaction_string_format = r'R\d{5}'
dm = dm.loc[dm['COMMENT'].str.contains(reaction_string_format)]

OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']

dm['COMMENT_FILTERED'] = dm['COMMENT'].str.replace(' OR ', ', ').str.split().explode().apply(
    scan_for_terms).dropna().groupby(level=0).apply(' '.join)
dm['OTHER_COMMENTS'] = dm['COMMENT'].str.split().explode().apply(remove_terms).dropna().groupby(level=0).apply(' '.join)
dm = dm.query('~COMMENT_FILTERED.str.contains("SIMILAR TO")')

to_remove = ['SEE ', '(', ')', 'REACTION ', '-STEP']
dm['info'] = dm['COMMENT_FILTERED'].copy()
for i in to_remove:
    dm['info'] = dm['info'].str.replace(i, '').str.lstrip().str.rstrip(', ')
dm['n_steps_cited'] = dm['info'].apply(lambda x: x[0]).replace('N', float('nan')).astype(float)
removed = dict()
for i in dm.index:
    noself = dm.loc[i, 'info'][1:].replace(dm.loc[i, 'index'] + '+', '').replace('+' + dm.loc[i, 'index'], '').replace(
        dm.loc[i, 'index'], '')
    if '+' in noself:
        removed[i] = [j.strip() for j in noself.split(',')]
dm['steps_shown'] = removed
dm = dm.explode('steps_shown').dropna(subset=['steps_shown']).sort_index().reset_index(drop=True)
dm['steps_shown'] = dm['steps_shown'].str.split('+')
dm['n_steps_shown'] = dm['steps_shown'].apply(len)

has_rns_missing = dm.loc[dm['steps_shown'].apply(flag_missing_reactions).dropna().index]['index'].to_list()
# the below looks weird - it's written so we don't remove reactions with multiple N-step sequences shown, at least one of which is complete.
incomplete_multistep = dm.loc[dm['index'].isin(has_rns_missing)]['index'].unique()
dm = dm.query('index not in @incomplete_multistep')
dm['step_group'] = dm['steps_shown'].apply('+'.join)

# uncomment to see some example raw-data formats:
# print(dm.loc[dm['OTHER_COMMENTS'].dropna().apply(len).sort_values(ascending=False).head(9).index]['ORIGINAL_COMMENT'])

# uncomment to see "passing" examples from each filtering step.
# dm.tail()


# Check: do all step counts match (where provided)?


# individual_reactions = dm.drop_duplicates(subset='index')
# print(individual_reactions['n_steps_shown'].value_counts().sort_index())
# print(((individual_reactions['n_steps_shown']-individual_reactions['n_steps_cited']).dropna()).value_counts())

# a = dm.dropna(subset='n_steps_cited').groupby(by='index')
# print('All step counts match: ', all(a['n_steps_cited'].unique() == a['n_steps_shown'].unique()))
# print('Step count distribution:\n', (a['n_steps_cited'].unique().explode()).value_counts())


# Flag these multistep "shortcuts" in the dataset, then generate data files.


reactions['MULTISTEP_FLAG'] = reactions.index.isin(dm['index'].unique())
reactions['STEP_ENTRIES'] = dm.groupby(by='index')['step_group'].apply(list)
reactions['STEP_ENTRIES'] = reactions['STEP_ENTRIES'].fillna('')

# make_lists_printable(dm).to_csv('../REACTIONS_MULTISTEP_INTERMEDIATE_PROCESSING.csv.zip', compression='zip', encoding='utf-8')

reactions_printable = make_lists_printable(reactions)
reactions_printable.to_csv('../REACTIONS_PROCESSED_FULL.csv.zip', compression='zip', encoding='utf-8')

reactions_processed = pd.concat([reactions_printable.filter(like='FLAG'), reactions_printable.filter(like='ENTRIES')],
                                axis=1)
reactions_processed.to_csv('../REACTIONS_PROCESSED_BASIC.csv.zip', compression='zip', encoding='utf-8')
