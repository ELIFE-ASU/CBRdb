# lets-get-kegg
Gets mol files from the KEGG database. 

You will want to start with lets_get_kegg.py. You can get KEGG C, KEGG D, and KEGG R. Note that by default, it will obtain the mol files for C or D. However, if you pass C_full or D_full, it will get the webpage info instead.

Next, you might want to convert the mole files into something more portable, such as smiles. mol_to_smiles.py will help you with this.

Next, let us say you want to get a working format of reactions. Use fix_reactions_data.py

You can get the ATLAS equation info using atlas_converter.py
