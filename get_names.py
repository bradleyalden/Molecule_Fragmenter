import json

def sort_list (sort_from, sort_to):

    resultsList = []

    with open(sort_from, "r") as file:
        for line in file:
            
            name = line.strip(", )\n(").split(", ")[0]
            resultsList.append(name)

    with open(sort_to, "w") as file:
        for result in resultsList:
            file.write(result + "\n")

if __name__ == "__main__":
    sort_list("SMARTS_MARGAN.py", "group_lists\\margan_names_only.txt")