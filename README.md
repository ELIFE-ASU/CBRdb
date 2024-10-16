# lets-get-kegg
Gets key information from the  mol files from the KEGG database. 

You will want to start with lets_get_kegg.py. You can get KEGG C, KEGG D, and KEGG R. Note that by default, it will obtain the mol files for C or D. However, if you pass C_full or D_full, it will get the webpage info instead.

Next, you might want to convert the mole files into something more portable, such as smiles. mol_to_smiles.py will help you with this.

Next, let us say you want to get a working format of reactions. Use fix_reactions_data.py

You can get the ATLAS equation info using atlas_converter.py

## File Generation explained
Anything that starts with KEGG is from KEGG. Anything that starts with ATLAS comes from ATLAS.
- kegg_data_C.csv.zip, this is your go-to for the processed and vetted compounds

Intermediate files
- C_IDs_bad
- C_IDs_good
- C_IDs_manual
- R_IDs_bad
