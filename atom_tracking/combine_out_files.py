import os
import csv

# Define the columns and their corresponding line prefixes
columns = [
    ("reaction_id", "reaction_id:"),
    ("tracking", "tracking:"),
    ("mapped_rxns", "mapped_rxns:"),
    ("reaction_no_stoich", "reaction_no_stoich:"),
    ("reaction_after_cofactor", "reaction_after_cofactor:"),
    ("selection_method", "selection_method:"),
    ("mapper_used", "mapper_used:"),
    ("cofactor_handling", "cofactor_handling:"),
    ("mapping_error", "mapping_error:")
]

out_dir = "/home/mshahjah/AtomTracking/Logs/out" 
output_dir = "."  

out_files = [f for f in os.listdir(out_dir) if f.endswith('.out')]
rows = []

for filename in out_files:
    data = {col: "" for col, _ in columns}
    with open(os.path.join(out_dir, filename), 'r', encoding='utf-8') as f:
        lines = f.readlines()[1:]  # skip first line
        for line in lines:
            for col, prefix in columns:
                if line.startswith(prefix):
                    data[col] = line[len(prefix):].strip()
    rows.append([data[col] for col, _ in columns])

with open(os.path.join(output_dir, 'combined_output_Remainings.csv'), 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([col for col, _ in columns])
    writer.writerows(rows)