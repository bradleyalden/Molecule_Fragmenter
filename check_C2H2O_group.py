import sys

def check_line_conditions(line):
    words = line.split()
    
    if len(words) == 1:
        return True

    if "179" not in words:
        return False
          
    if "169" not in words:
        return False
        
    for i, word in enumerate(words):
        if word == "169":
            if i > 0:
                    # Convert the preceding word to a float to handle both ints and decimals
                    prev_number = int(words[i-1])
                    if prev_number > 1:
                        return True
                    
    return False

def process_file(check_from,print_to):
    with open(check_from, 'r', encoding='utf-8') as file:
        lines = []
        for line in file:
            if not line.strip():
                continue
            if check_line_conditions(line):
                print(f"PASS: {line.strip()}")
                lines.append(line)
    with open(print_to, "w") as file:
        for index, line in enumerate(lines):
            file.write(line)


# --- Test Block for your specific examples ---
if __name__ == "__main__":
    process_file("group_lists\\fragmented_groups_sorted_full_margan.txt","group_lists\\C2H2O_possible_fails.txt")