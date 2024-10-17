# lets-get-kegg
Gets key information from the KEGG database. To start everything locally, you can just look at run_all.py, which will take you over the whole data pipeline.

Otherwise, you can run each part by themselves. 

## File Generation explained
Anything that begins with KEGG is from KEGG. Anything that starts with ATLAS comes from ATLAS.
- KEGG reactions
- ATLAS expansion
- ATLAS's 2018 KEGG pull

### Processed files
  
| Final files | Description |
| :---        |    :----   |
| kegg_data_C.csv.zip | This is your go-to for the processed and vetted compounds |
| kegg_data_R_processed.csv.zip |  |
| atlas_data_R_processed.csv.zip |  |
| .csv.zip |  |
| .csv.zip |  |
| .csv.zip |  |

### Merged files

| Final files | Description |
| :---        |    :----   |
| full_processed_merged.csv.zip | kegg_data_R_processed (2024) + atlas_data_R_processed = kegg_atlas_processed_merged |
| full_processed_merged.csv.zip | atlas_data_kegg_R_processed (2018) + atlas_data_R_processed = atlas_kegg_processed_merged |
| .csv.zip |  |
| .csv.zip |  |
| .csv.zip |  |

Along the way, you might generate a list of files corresponding to bad or good compound or reaction IDs.

### Intermediate files

| Intermediate files | Description |
| :---        |    :----   |
| C_IDs_bad.dat | Bad compound IDs to skip       |
| C_IDs_good.dat   | Good compound IDs to follow up on        |
| C_IDs_manual.dat   | Based off of C_IDs_good these correspond to manual fixes |
| R_IDs_bad.dat   | List is a list of compounds that have been marked as bad and are not included in final data sets |

