import ast

def convert_py_list_to_txt(py_file_path, output_txt_path):
    try:
        # 1. Read the content of the .py file
        with open(py_file_path, 'r', encoding='utf-8') as f:
            tree = ast.parse(f.read())

        # 2. Find the assignment to 'UNIFAC'
        unifac_data = None
        for node in tree.body:
            if isinstance(node, ast.Assign):
                # Check if the variable being assigned to is named 'UNIFAC'
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == 'MARGAN':
                        # Convert the AST list node back into a real Python list
                        unifac_data = ast.literal_eval(node.value)
                        break
        
        if unifac_data is None:
            print("Error: Could not find a variable named 'UNIFAC' in the file.")
            return

        # 3. Write the indexed items to the .txt file
        with open(output_txt_path, 'w', encoding='utf-8') as f:
            for index, item in enumerate(unifac_data, start=1):
                # item is a tuple, e.g., ("CH3", "[CH3;X4]")
                f.write(f"{index} {item}\n")

        print(f"Successfully processed {len(unifac_data)} items into {output_txt_path}")

    except FileNotFoundError:
        print("Error: The .py file was not found.")
    except SyntaxError:
        print("Error: The .py file contains invalid Python syntax.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# --- Execution ---
# Replace 'input_file.py' with your actual filename
convert_py_list_to_txt('SMARTS_MARGAN.py', 'SMARTS_MARGAN_groups_list.txt')