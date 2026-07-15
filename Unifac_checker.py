import csv

simon_array = []
rchat_array = []

with open('Unifac_simon.csv', mode='r', encoding='utf-8') as simon:
    simon_reader = csv.DictReader(simon)

    for row in simon_reader:
        simon_array.append(row['Name'])


with open('Unifac_rchat.csv', mode='r', encoding='utf-8') as rchat:
    rchat_reader = csv.DictReader(rchat)

    for row in rchat_reader:
        rchat_array.append(row['Name'])

for index, entry in enumerate(rchat_array):
    if entry != simon_array[index]:
        print("Difference at ", index + 1, ":\n", entry, "\n", simon_array[index], "\n")


# simon's descriptors:
# 3 is the number of heavy atoms
# 5 is non C/H atoms
# 6 is true for aromatic/ring structures (see group 27)
# 7 is number of triple bonds
# 8 is number of double bonds

# line 20 of the real solution should be CH=O not CHO

#lines 110-112 of simon's are extra: NCO, (CH2)2SU, CH2CHSU