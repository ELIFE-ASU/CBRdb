# Input: An unzipped, locally-downloaded folder of raw KEGG reaction data files.
# Output: CSVs with relevant reaction info, e.g. ties to other databases and flags for possible removal.
import os

import pandas as pd

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def scan_for_kegg_pattern(s, prefix):
    if prefix in s:
        as_list = s.replace(';', ' ').replace('[', ' ').replace(']', ' ').split(' ')
        return [i for i in as_list if i.startswith(prefix) and len(i) == len(prefix) + 5]
    else:
        return []


def make_lists_printable(df):
    list_like_cols = df.loc[:, df.iloc[0].apply(type).apply(lambda x: x in [set, list])]
    printable_cols = list_like_cols.map(' '.join)
    printable_df = df.copy(deep=True)
    printable_df.loc[:, printable_cols.columns] = printable_cols
    return printable_df


def scan_for_terms(s):
    for j in OK:
        if j in s or s.startswith('R1') or s.startswith('R0'):
            return s
    else:
        return None


def remove_terms(s):
    for j in OK:
        if j in s or s.startswith('R1') or s.startswith('R0'):
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


def load_reactions_data(target_dir):
    a = os.listdir(os.path.abspath(target_dir))
    reactions = {}
    for i in a:
        try:
            with open(f"{target_dir}{i}/{i}.data", 'r') as f:
                entries = {}
                for line in f.readlines()[:-1]:
                    if not line.startswith(' '):
                        field = line.split(' ')[0]
                        entries[field] = line.lstrip(field).lstrip().rstrip('\n')
                    else:
                        entries[field] += '; ' + line.lstrip().rstrip('\n')
                reactions[i] = entries
        except:
            pass
    return pd.DataFrame(reactions).fillna('').T


def main():
    global reactions
    global OK

    f_save_files = False

    # This needs to point to the directory where the KEGG reaction data files are stored.
    target_dir = '../data/kegg_data_R/'
    data_dir = "Data/"
    reactions = load_reactions_data(target_dir)

    # Parse each reaction's ties to other databases.
    # Flag reactions with glycan participants and "overall" reactions (see [this table](https://www.kegg.jp/kegg/tables/br08210.html) for info).
    reactions = reactions.assign(
        OVERALL_FLAG=reactions['ENTRY'].str.contains('Overall', case=False),
        GLYCAN_FLAG=reactions['EQUATION'].str.contains('G'),
        EQUATION_ENTRIES=reactions['EQUATION'].str.findall(r'C\d{5}|G\d{5}'),
        ENZYME_ENTRIES=reactions['ENZYME'].str.replace(';', ' ').str.split(),
        DBLINKS_RHEA_ENTRIES=reactions['DBLINKS'].str.lstrip('RHEA: ').str.split())
    prefixes = {'RCLASS': r'RC\d{5}',
                'BRITE': r'br\d{5}',
                'PATHWAY': r'rn\d{5}',
                'MODULE': r'M\d{5}',
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

    if f_save_files:
        rclass_db = r_rclass.reset_index(names='REACTION').groupby(by=['RC', 'CPAIR'])['REACTION'].apply(set).apply(
            list).to_frame()
        make_lists_printable(rclass_db).to_csv(os.path.join(data_dir, 'rclass_db.csv'))

    # Use comments to identify multi-step reactions
    # i.e., artificial shortcuts compressing multiple existing reaction IDs into a single net reaction.
    to_replace = {'MULTISTEP': 'MULTI-STEP',
                  '3STEP': '3-STEP',
                  '; SEE R': ' (SEE R',
                  'REACTION;': 'REACTION '}
    dm = reactions['COMMENT'].replace('', float('NaN')).dropna().str.upper().reset_index()
    dm['ORIGINAL_COMMENT'] = dm['COMMENT'].copy(deep=True)
    for k, v in to_replace.items():
        dm['COMMENT'] = dm['COMMENT'].str.replace(k, v)
    dm['COMMENT'] = dm['COMMENT'].str.split(';')
    s1 = r'COMMENT.str.contains("STEP REACTION")'
    s2 = r'COMMENT.str.contains(r"R\d{5}\+|R\d{5} \+|R\d{5} \,")'
    s3 = r'COMMENT.str.contains("STEP OF|PART OF|THE LAST|THE FIRST|THE FORMER|THE LATTER|POSSIBLY|PROBABLY IDENTICAL")'
    dm = dm.explode('COMMENT').reset_index(drop=True).query(s1 + ' and ' + s2 + ' and not ' + s3)
    to_replace = {'TWO': '2',
                  'THREE': '3',
                  'FOUR': '4',
                  'MULTI': 'N',
                  ' + ': '+'}
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
    dm['OTHER_COMMENTS'] = dm['COMMENT'].str.split().explode().apply(remove_terms).dropna().groupby(level=0).apply(
        ' '.join)
    dm = dm.query('~COMMENT_FILTERED.str.contains("SIMILAR TO")')

    to_remove = ['SEE ', '(', ')', 'REACTION ', '-STEP']
    dm['info'] = dm['COMMENT_FILTERED'].copy()
    for i in to_remove:
        dm['info'] = dm['info'].str.replace(i, '').str.lstrip().str.rstrip(', ')
    dm['n_steps_cited'] = dm['info'].apply(lambda x: x[0]).replace('N', float('nan')).astype(float)
    removed = dict()
    for i in dm.index:
        noself = dm.loc[i, 'info'][1:].replace(dm.loc[i, 'index'] + '+', '').replace('+' + dm.loc[i, 'index'],
                                                                                     '').replace(
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
    print('Reactions with incomplete multistep shortcuts:', len(incomplete_multistep), flush=True)
    print(incomplete_multistep, flush=True)

    dm = dm.query('index not in @incomplete_multistep')
    dm['step_group'] = dm['steps_shown'].apply('+'.join)

    # Flag these multistep "shortcuts" in the dataset, then generate data files.
    reactions['MULTISTEP_FLAG'] = reactions.index.isin(dm['index'].unique())
    reactions['STEP_ENTRIES'] = dm.groupby(by='index')['step_group'].apply(list)
    reactions['STEP_ENTRIES'] = reactions['STEP_ENTRIES'].fillna('')
    reactions_printable = make_lists_printable(reactions)
    if f_save_files:
        make_lists_printable(dm).to_csv('../reactions_multistep_intermediate_processing.csv.zip', compression='zip',
                                        encoding='utf-8')
        reactions_printable.to_csv(os.path.join(data_dir, 'reactions_processed_full.csv.zip'), compression='zip',
                                   encoding='utf-8')

    reactions_processed = pd.concat(
        [reactions_printable.filter(like='FLAG'),
         reactions_printable.filter(like='ENTRIES')],
        axis=1)
    if f_save_files:
        reactions_processed.to_csv(os.path.join(data_dir, 'reactions_processed_basic.csv.zip'), compression='zip',
                                   encoding='utf-8')

    # select only the reactions that are overall_flagged or glycan_flagged
    reactions_overall = reactions_printable.query('OVERALL_FLAG == True')
    reactions_glycan = reactions_printable.query('GLYCAN_FLAG == True')
    reactions_overall_glycan = pd.concat([reactions_overall, reactions_glycan])
    # print out the ids of the reactions that are overall_flagged or glycan_flagged
    bad_list = reactions_overall_glycan['ENTRY'].to_list()
    bad_list = [i.split()[0].strip() for i in bad_list]

    out = ['R13204', 'R13205', 'R13206', 'R13207', 'R13208', 'R13209', 'R13210', 'R13211', 'R13222', 'R01331', 'R05901',
           'R05902', 'R05903', 'R05904', 'R05905', 'R05906', 'R05907', 'R05908', 'R05909', 'R05910', 'R05911', 'R05912',
           'R05913', 'R05914', 'R05915', 'R05916', 'R05917', 'R05918', 'R05919', 'R05922', 'R05923', 'R05924', 'R05925',
           'R05926', 'R05927', 'R05928', 'R05929', 'R05930', 'R05931', 'R05932', 'R05933', 'R05934', 'R05935', 'R05936',
           'R05937', 'R05938', 'R05939', 'R05940', 'R05941', 'R05942', 'R05943', 'R05944', 'R05945', 'R05946', 'R05947',
           'R05948', 'R05949', 'R05950', 'R05951', 'R05952', 'R05953', 'R05954', 'R05955', 'R05956', 'R05957', 'R05958',
           'R05959', 'R05960', 'R05961', 'R05962', 'R05963', 'R05964', 'R05965', 'R05966', 'R05967', 'R05968', 'R05969',
           'R05970', 'R05971', 'R05972', 'R05973', 'R05974', 'R05975', 'R05976', 'R05977', 'R05978', 'R05979', 'R05980',
           'R05981', 'R05983', 'R05984', 'R05985', 'R05986', 'R05987', 'R05988', 'R05989', 'R05990', 'R05991', 'R05992',
           'R05994', 'R05995', 'R05996', 'R05998', 'R05999', 'R06000', 'R06001', 'R06003', 'R06004', 'R06005', 'R06006',
           'R06007', 'R06009', 'R06010', 'R06011', 'R06012', 'R06013', 'R06014', 'R06015', 'R06016', 'R06018', 'R06020',
           'R06021', 'R06022', 'R06023', 'R06024', 'R06025', 'R06026', 'R06027', 'R06028', 'R06029', 'R06030', 'R06031',
           'R06032', 'R06033', 'R06034', 'R06035', 'R06036', 'R06037', 'R06038', 'R06039', 'R06040', 'R06041', 'R06042',
           'R06043', 'R06044', 'R06045', 'R06049', 'R06050', 'R06051', 'R06052', 'R06053', 'R06055', 'R06056', 'R06057',
           'R06058', 'R06059', 'R06061', 'R06062', 'R06066', 'R06067', 'R06068', 'R06069', 'R06070', 'R06071', 'R06072',
           'R06073', 'R06074', 'R06075', 'R06076', 'R06077', 'R06078', 'R06079', 'R06080', 'R06081', 'R06082', 'R06083',
           'R06084', 'R06085', 'R06086', 'R06087', 'R06088', 'R06089', 'R06090', 'R06091', 'R06092', 'R06093', 'R06094',
           'R06095', 'R06096', 'R06097', 'R06098', 'R06099', 'R06100', 'R06101', 'R06102', 'R06103', 'R06104', 'R06105',
           'R06106', 'R06107', 'R06108', 'R06109', 'R06110', 'R06112', 'R06113', 'R06114', 'R06115', 'R06116', 'R06118',
           'R06121', 'R06122', 'R06125', 'R06127', 'R06128', 'R06129', 'R06130', 'R06135', 'R06139', 'R06140', 'R06141',
           'R06142', 'R06143', 'R06144', 'R06145', 'R06147', 'R06149', 'R06150', 'R06151', 'R06152', 'R06153', 'R06155',
           'R06156', 'R06157', 'R06158', 'R06160', 'R06161', 'R06162', 'R06163', 'R06164', 'R06165', 'R06167', 'R06168',
           'R06169', 'R06170', 'R06172', 'R06173', 'R06174', 'R06175', 'R06176', 'R06177', 'R06179', 'R06181', 'R06182',
           'R06183', 'R06184', 'R06185', 'R06186', 'R06187', 'R06188', 'R06189', 'R06190', 'R06191', 'R06192', 'R06193',
           'R06197', 'R06198', 'R06199', 'R06200', 'R06201', 'R06202', 'R06203', 'R06204', 'R06206', 'R06207', 'R06208',
           'R06209', 'R06210', 'R06211', 'R06213', 'R06214', 'R06215', 'R06216', 'R06217', 'R06218', 'R06219', 'R06221',
           'R06222', 'R06224', 'R06225', 'R06226', 'R06227', 'R06228', 'R06229', 'R06230', 'R06232', 'R06233', 'R06234',
           'R06236', 'R06237', 'R06238', 'R06239', 'R06240', 'R06241', 'R06242', 'R06243', 'R06244', 'R06245', 'R06246',
           'R06247', 'R06249', 'R06250', 'R06251', 'R06253', 'R06258', 'R06259', 'R06260', 'R06261', 'R06262', 'R06263',
           'R06264', 'R06273', 'R06274', 'R06275', 'R06276', 'R06277', 'R06278', 'R06279', 'R06280', 'R06283', 'R06285',
           'R06289', 'R06290', 'R06623', 'R06722', 'R06986', 'R07129', 'R07131', 'R07609', 'R07610', 'R07611', 'R07614',
           'R07615', 'R07616', 'R07617', 'R07619', 'R07620', 'R07621', 'R07622', 'R07623', 'R07624', 'R07625', 'R07628',
           'R07805', 'R07806', 'R07807', 'R07808', 'R07809', 'R07810', 'R07811', 'R07812', 'R07813', 'R07814', 'R07815',
           'R07816', 'R07817', 'R07818', 'R07819', 'R07820', 'R07821', 'R07822', 'R07823', 'R07824', 'R07825', 'R08107',
           'R08599', 'R08717', 'R08718', 'R09290', 'R09295', 'R09296', 'R09297', 'R09298', 'R09299', 'R09300', 'R09301',
           'R09302', 'R09303', 'R09304', 'R09315', 'R09316', 'R09318', 'R09319', 'R09320', 'R09323', 'R09324', 'R09764',
           'R09766', 'R10138', 'R10139', 'R10830', 'R10831', 'R10905', 'R11316', 'R11317', 'R11318', 'R11320', 'R11321',
           'R11322', 'R11375', 'R11376', 'R11377', 'R11378', 'R11379', 'R11407', 'R11994', 'R11995', 'R11997', 'R11998',
           'R12038', 'R12039', 'R12040', 'R12041', 'R12042', 'R12075', 'R12077', 'R12078', 'R12079', 'R12080', 'R12244',
           'R12292', 'R12293', 'R12294', 'R12295', 'R12296', 'R12297', 'R12298', 'R12316', 'R12350', 'R12351', 'R12477',
           'R12479', 'R12482', 'R12637', 'R12702', 'R12704', 'R12705', 'R12738', 'R12739', 'R12740', 'R12741', 'R12742',
           'R12743', 'R12744', 'R12745', 'R12808', 'R12810', 'R12811', 'R12812', 'R12813', 'R12814', 'R12815', 'R12816',
           'R12817', 'R12818', 'R12819', 'R12820', 'R12821', 'R12822', 'R12823', 'R12859', 'R12860', 'R12861', 'R12862',
           'R12863', 'R12864', 'R12865', 'R12866', 'R12867', 'R12868', 'R12869', 'R12870', 'R12871', 'R12872', 'R12873',
           'R12874', 'R12875', 'R12898', 'R12899', 'R12900', 'R12901', 'R12902', 'R12903', 'R12904', 'R12905', 'R12906',
           'R12944', 'R12951', 'R13003', 'R13004', 'R13005', 'R13169', 'R13171', 'R13172', 'R13173', 'R13174', 'R13175',
           'R13248', 'R13390', 'R13391', 'R13392', 'R13393']
    assert bad_list == out, f"Expected does not match actual"


if __name__ == "__main__":
    print("Program started", flush=True)
    main()
    print("Program finished", flush=True)
