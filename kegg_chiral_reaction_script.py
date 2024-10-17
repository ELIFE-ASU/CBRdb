
import pandas as pd
import re

# Load the compound and reaction files
def load_files(compound_file, reaction_file):
    kegg_data_C = pd.read_csv(compound_file)
    kegg_data_R = pd.read_csv(reaction_file)
    return kegg_data_C, kegg_data_R

# Create a dictionary from kegg_data_C for quick lookup of compound chiral centers
def create_compound_dict(kegg_data_C):
    return kegg_data_C.set_index('compound_id')['n_chiral_centers'].to_dict()

# Function to parse reactions and extract compounds with their stoichiometry
def parse_reaction_for_stoichiometry(reaction_str):
    def parse_compounds(compound_str):
        compounds = compound_str.split('+')
        parsed_compounds = []
        stoichiometries = []
        for compound in compounds:
            compound = compound.strip()
            # Match numeric stoichiometry and compound ID
            match = re.match(r'(\d*)\s*(C\d{5})', compound)
            if match:
                stoichiometry = int(match.group(1)) if match.group(1) else 1
                compound_id = match.group(2)
                parsed_compounds.append((compound_id, stoichiometry))
                stoichiometries.append(str(stoichiometry))
            else:
                # Handle non-numeric stoichiometry like "n", "n+1"
                parsed_compounds.append((compound, "NA"))
                stoichiometries.append("NA")
        return parsed_compounds, stoichiometries

    reactants_str, products_str = reaction_str.split('<=>')
    reactants, reactant_stoichiometries = parse_compounds(reactants_str.strip())
    products, product_stoichiometries = parse_compounds(products_str.strip())

    return reactants, products, reactant_stoichiometries, product_stoichiometries

# Function to label reaction type based on chiral centers
def label_reaction(reactants_chiral, products_chiral):
    if reactants_chiral == 0 and products_chiral == 0:
        return "achiral conserving"
    elif reactants_chiral == products_chiral:
        return "chiral conserving"
    else:
        return "chiral modulating"

# Function to calculate chiral centers with error handling
def calculate_chiral_centers_with_error_handling(compounds, compound_chiral_dict):
    total_chiral_centers = 0
    missing_compounds = []
    for compound_id, stoichiometry in compounds:
        if compound_id.startswith('C') and compound_id in compound_chiral_dict:
            chiral_centers = compound_chiral_dict[compound_id]
            if stoichiometry != "NA":
                total_chiral_centers += chiral_centers * stoichiometry
        else:
            missing_compounds.append(compound_id)
    return total_chiral_centers, missing_compounds

# Function to calculate the fraction of chiral compounds
def calculate_fraction_chiral(reactants, products, compound_chiral_dict):
    total_chiral_compounds = 0
    total_compounds = 0
    
    for compound_id, stoichiometry in reactants + products:
        chiral_centers = compound_chiral_dict.get(compound_id, 0)
        if chiral_centers > 0:
            total_chiral_compounds += 1 * stoichiometry
        total_compounds += 1 * stoichiometry
    
    if total_compounds == 0:
        return 0  # Handle division by zero gracefully
    
    return total_chiral_compounds / total_compounds

# Main function to process reactions
def process_reactions(kegg_data_R, compound_chiral_dict):
    reaction_labels = []
    additional_data = []
    stoichiometry_column = []
    notes_column = []

    for idx, row in kegg_data_R.iterrows():
        reaction_str = row['reaction']
        
        try:
            # Parse the reaction to extract stoichiometric values
            reactants, products, reactant_stoichiometries, product_stoichiometries = parse_reaction_for_stoichiometry(reaction_str)
            
            # Join stoichiometries of reactants and products into a comma-separated string
            stoichiometry_value = ', '.join(reactant_stoichiometries + product_stoichiometries)
            
            # Calculate chiral centers for reactants and products, handle missing compounds
            reactants_chiral, reactants_missing = calculate_chiral_centers_with_error_handling(reactants, compound_chiral_dict)
            products_chiral, products_missing = calculate_chiral_centers_with_error_handling(products, compound_chiral_dict)
            
            # If there is any issue with missing compounds or non-numeric stoichiometry, mark the rxn_label as "NA"
            if reactants_missing or products_missing:
                reaction_label = "NA"
                # chiral_order_param = "NA"  # Commenting out chiral order param
                # homochirality_order_param = "NA"  # Commenting out homochirality param
                fraction_chiral = "NA"
                note = f"Missing compounds: {reactants_missing + products_missing}"
            else:
                # Calculate labels and parameters normally
                reaction_label = label_reaction(reactants_chiral, products_chiral)
                # chiral_order_param = calculate_chiral_order_parameter(reactants_chiral, products_chiral)  # Commented out
                # homochirality_order_param = calculate_homochirality_order_parameter(reactants_chiral, products_chiral)  # Commented out
                fraction_chiral = calculate_fraction_chiral(reactants, products, compound_chiral_dict)
                note = ""

            # Append results
            reaction_labels.append([reaction_label])
            additional_data.append([fraction_chiral])  # Only adding fraction_chiral to additional data
            stoichiometry_column.append(stoichiometry_value)
            notes_column.append(note)

        except Exception as e:
            reaction_labels.append(["NA"])
            additional_data.append(["NA"])
            stoichiometry_column.append("Error")
            notes_column.append(f"Error: {str(e)}")

    # Add the calculated labels, stoichiometry column, and notes to the dataframe
    kegg_data_R['rxn_labels'] = [r[0] for r in reaction_labels]
    kegg_data_R['fraction_chiral'] = [d[0] for d in additional_data]
    kegg_data_R['Stoichiometry'] = stoichiometry_column
    kegg_data_R['note'] = notes_column

    return kegg_data_R

# Save processed data to CSV
def save_to_csv(kegg_data_R, output_file):
    kegg_data_R.to_csv(output_file, index=False)

# Example usage
if __name__ == "__main__":
    compound_file = input("Enter the path to the compound file (kegg_data_C.csv): ")
    reaction_file = input("Enter the path to the reaction file (kegg_data_R.csv): ")
    output_file = input("Enter the path for the output file: ")

    kegg_data_C, kegg_data_R = load_files(compound_file, reaction_file)
    compound_chiral_dict = create_compound_dict(kegg_data_C)
    processed_data = process_reactions(kegg_data_R, compound_chiral_dict)
    save_to_csv(processed_data, output_file)
    print(f"Processed data saved to {output_file}")
