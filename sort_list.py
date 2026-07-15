import json

def sort_list (sort_from, sort_to):

    resultsDict = {}

    with open(sort_from, "r") as file:
        for line in file:
            # if line.find(" ") == -1:
            #     continue
            id = int(line.split()[0])
            try:
                theRest = line.split(None,1)[1]
            except:
                theRest = ""
            resultsDict[id] = theRest

    sortedDict = dict(sorted(resultsDict.items()))

    with open(sort_to, "w") as file:
        for index, (id, theRest) in enumerate(sortedDict.items()):
            file.write(str(id) + " " + theRest.strip())
            if index < len(sortedDict) - 1:
                file.write("\n")

if __name__ == "__main__":
    sort_list("group_lists\\fragmented_groups.txt", "group_lists\\fragmented_groups_sorted_full_margan.txt")