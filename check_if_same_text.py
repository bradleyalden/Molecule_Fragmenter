def compare_files(file1_path, file2_path):
    try:
        with open(file1_path, 'r', encoding='utf-8') as f1, \
             open(file2_path, 'r', encoding='utf-8') as f2:
            
            if f1.read() == f2.read():
                print("The files are identical.")
            else:
                print("The files are different.")
                
    except FileNotFoundError:
        print("Error: One or both files were not found.")

# Usage
compare_files("fragmented_groups_sorted_tde_groups.txt", "fragmented_groups_sorted_nolimit.txt")