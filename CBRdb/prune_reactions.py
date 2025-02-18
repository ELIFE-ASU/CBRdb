import pandas as pd 

def all_entries(dbs):
    kegg_rns = dbs['kegg_data_R'].set_index('id')['reaction']
    atlas_rns = dbs['atlas_data_R'].set_index('id')['reaction']
    all_rns = pd.concat([kegg_rns, atlas_rns]).str.findall(r"(C\d{5})").map(set)
    all_cps = set(dbs['CBRdb_C']['compound_id'])
    return(all_rns, all_cps)

def all_kegg_comments(dbs):
    kegg_cmts = dbs['kegg_data_R'].set_index('id')['comment'].str.replace("reaction;see", "reaction (see")
    kegg_cmts = kegg_cmts.str.split(';').explode().dropna().to_frame()
    kegg_cmts['rn_refs'] = kegg_cmts['comment'].str.upper().str.findall(r"(R\d{5})").map(set).sub(
        kegg_cmts.index.map(lambda x: {x}))
    kegg_cmts.at['R10693','comment'] = 'part of '+kegg_cmts.at['R10693','comment']
    return kegg_cmts

def list_reactions_with_halogen_dupes(dbs):
    all_rns, all_cps = all_entries(dbs)
    cps_duped = dbs['CBRdb_C'].query('compound_id.str.contains("C99") & smiles.duplicated(keep=False)')['compound_id'].values
    rns = all_rns[~all_rns.map(lambda x: x.isdisjoint(cps_duped))].index.drop_duplicates()
    return rns

def list_reactions_missing_structures(dbs):
    all_rns, all_cps = all_entries(dbs)
    rns = all_rns[~all_rns.map(lambda x: x.issubset(all_cps))].index.drop_duplicates()
    return rns

def list_reactions_with(dbs, phrase):
    kegg_cmts = all_kegg_comments(dbs)
    query = f'comment.str.lower().str.contains("{phrase}")'
    rns = kegg_cmts.query(query).index.drop_duplicates()
    return rns

def list_multistep_parts(dbs): #KEEP THESE
    cmts = all_kegg_comments(dbs)['comment']
    rns = cmts[cmts.str.contains(" of ") & cmts.str.contains("step")].index.drop_duplicates()
    return rns

def list_multistep_enumerated(dbs):
    reactions_overall = dbs['kegg_data_R'].set_index('id')['overall'].dropna().index
    reactions_multistep_parts = list_multistep_parts(dbs) #KEEP THESE
    rns = all_kegg_comments(dbs).drop(reactions_multistep_parts).query(
        'rn_refs.str.len()>1 & comment.str.contains("step") & ~comment.str.contains("possibl|probabl|similar")')
    rns = rns.index.union(reactions_overall).drop_duplicates().union({'R10671'}) # false negative: "similar" in line
    return rns

def df_of_suspect_reactions(dbs):
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
    sus = pd.concat([pd.Series(index=v, data=len(v)*[k]) for k, v in sus.items()])
    sus = sus.groupby(level=0).apply(list).map(lambda x: '+'.join(sorted(x))).to_frame('reason')
    return sus

def suspect_reaction_subset(sus, matching):
    return sus.map(lambda x: x.split('+')).explode('reason').query("reason.str.contains(@matching)").index.drop_duplicates()

def add_suspect_reactions_to_existing_bad_file(sus_new, R_IDs_bad_file="../data/R_IDs_bad.dat"):
    sus_old = pd.read_csv(R_IDs_bad_file, header=0, index_col=0)['reason'].str.split('+').explode()
    sus_new = sus_new['reason'].str.split('+').explode()
    sus_combo = (pd.concat([sus_old, sus_new]).groupby(level=0)
                 .apply(lambda x: sorted(list(set(x)))).map('+'.join).to_frame('reason')).sort_index()
    sus_combo.to_csv(R_IDs_bad_file)
    return sus_combo

def quarantine_suspect_reactions_matching(dbs, sus, matching="shortcut|structure_missing"):
    to_quarantine = suspect_reaction_subset(sus, matching)
    dbs['sus'] = sus
    dbs['quarantine_categories'] = matching
    dbs['quarantined'] = (dbs['kegg_data_R'].merge(dbs['atlas_data_R'], how='outer', on='id')
                          .query('id.isin(@to_quarantine)').copy(deep=True))
    dbs['kegg_data_R_raw'] = dbs['kegg_data_R'].copy(deep=True)
    dbs['atlas_data_R_raw'] = dbs['atlas_data_R'].copy(deep=True)
    dbs['kegg_data_R'] = dbs['kegg_data_R'].query('~id.isin(@to_quarantine)')
    dbs['atlas_data_R'] = dbs['atlas_data_R'].query('~id.isin(@to_quarantine)')
    return dbs

#Usage:
#dbs = CBRdb.dedupe_compounds()  #after deduping...
#sus = CBRdb.df_of_suspect_reactions(dbs) # identify suspect reactions
#sus = CBRdb.add_suspect_reactions_to_existing_bad_file(sus) #optional: add to log and import log
#dbs = CBRdb.quarantine_suspect_reactions_matching(dbs, sus, matching="shortcut|structure_missing")
#write CSVs for intermediate output files as desired, otherwise continue to use DFs