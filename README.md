# lets-get-kegg
Gets key information from the KEGG database. To start everything locally, you can just look at run_all.py, which will take you over the whole data pipeline. Otherwise, you can run each part by themselves. 



## File Generation explained
Anything that begins with KEGG is from KEGG. Anything that starts with ATLAS comes from ATLAS.
- KEGG reactions (2024)
- ATLAS reactions expansion
- ATLAS's KEGG reactions pull (2018)

## Un-Processed files

These are considered raw and unprocessed. Please be sure to use caution.

| Final files | Description |
| :---        |    :----   |
| kegg_data_R.csv.zip | KEGG reactions (2024) list |
| atlas_data_R.csv.zip | ATLAS reactions expansion |
| atlas_data_kegg_R.csv.zip | ATLAS's KEGG reactions pull (2018) |

## Processed files

These are processed files that have been filtered and processed.

| Final files | Description |
| :---        |    :----   |
| kegg_data_C.csv.zip | This is your go-to for processed and vetted compounds. |
| ec_ids.csv.zip | This is a complete list of ECs for all the reactions. |
| kegg_data_R_processed.csv.zip | Processed KEGG reactions (2024) list |
| atlas_data_R_processed.csv.zip | Processed ATLAS reactions expansion |
| atlas_data_kegg_R_processed.csv.zip | ATLAS's KEGG reactions pull (2018) |

## Merged files

These are files that have been merged for convenience.

| Final files | Description |
| :---        |    :----   |
| atlas_kegg_processed_merged.csv.zip | kegg_data_R_processed (2024) + atlas_data_R_processed = kegg_atlas_processed_merged |
| kegg_atlas_processed_merged.csv.zip | atlas_data_kegg_R_processed (2018) + atlas_data_R_processed = atlas_kegg_processed_merged |



## Intermediate files

Along the way, you might generate a list of files corresponding to bad or good compound or reaction IDs.

| Intermediate files | Description |
| :---        |    :----   |
| C_IDs_bad.dat | Bad compound IDs to skip |
| C_IDs_good.dat | Good compound IDs to follow up on |
| C_IDs_manual.dat | Based off of C_IDs_good these correspond to manual fixes |
| R_IDs_bad.dat | A list of reactions that have been marked as bad and are not included in final data sets |

