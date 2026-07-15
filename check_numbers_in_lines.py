import sys

# checks if there are two of the first number in one line (id matches a group in the solution)

def check_duplicate_first_word(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            for line in file:
                # Split the line into a list of words based on whitespace
                words = line.split()
                
                # Skip completely blank lines
                if not words:
                    continue
                
                # Identify the first word
                first_word = words[0]
                
                # Count how many times the first word appears in the list of words
                # If it appears less than 2 times, print the original line
                if words.count(first_word) < 2:
                    print(line.strip())
                    
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Ensure the user provided a file to check
    check_duplicate_first_word("group_lists\\fragmented_groups.txt")