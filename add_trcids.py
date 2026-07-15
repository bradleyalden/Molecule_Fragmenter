import csv

def process_files(txt_input_path, csv_input_path, txt_output_path):
    # 1. Load CSV into a dictionary for fast lookup
    # Key: smiles string, Value: trcid
    smiles_map = {}
    
    try:
        with open(csv_input_path, mode='r', encoding='utf-8') as csv_file:
            # Using DictReader assumes the first row contains headers "trcid" and "smiles"
            reader = csv.DictReader(csv_file)
            for row in reader:
                smiles_map[row['smiles']] = row['trcid']
    except FileNotFoundError:
        print("Error: CSV file not found.")
        return
    except KeyError:
        print("Error: CSV must have columns named 'trcid' and 'smiles'.")
        return

    # 2. Process the text file and write to output
    try:
        with open(txt_input_path, mode='r', encoding='utf-8') as txt_in, \
             open(txt_output_path, mode='w', encoding='utf-8') as txt_out:
            
            for line in txt_in:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                
                # Split the line at the first occurrence of ' {' 
                # to isolate the 'string' from the '{dict}'
                parts = line.split(' {', 1)
                if len(parts) == 2:
                    smiles_string = parts[0].strip()
                    dict_part = ' {' + parts[1]
                    
                    # Look up the string in our CSV map
                    if smiles_string in smiles_map:
                        trcid = smiles_map[smiles_string]
                        # Write: trcid string {dict}
                        txt_out.write(f"{trcid} {smiles_string}{dict_part}\n")
                    else:
                        # If not found in CSV, you can choose to keep the line as is 
                        # or skip it. Here we keep it as is.
                        txt_out.write(f"{line}\n")
                else:
                    # Line doesn't follow the expected "string {dict}" format
                    txt_out.write(f"{line}\n")
                    
        print(f"Success! Processed file saved to: {txt_output_path}")

    except FileNotFoundError:
        print("Error: TXT file not found.")

# --- Execution ---
# Change these filenames to match your actual files
process_files('fragmented_groups_real_solutions.txt', 'trcid_smiles_6_22_2026.csv', 'fragmented_groups_trcid_and_names.txt')