import linecache

with open("Mar Gan Groups/groups_only.txt", "w", encoding="utf-8") as new:
    with open("Mar Gan Groups/groups_raw.txt", "r", encoding="utf-8") as raw:
        for index, raw_line in enumerate(raw):
            corrected_line = linecache.getline("groups_only_corrected.txt", index + 1)
            to_write = raw_line.strip().split(",")
            to_write[1] = corrected_line.strip()
            new.write(to_write[0] + "," + to_write[1])
            new.write("\n")