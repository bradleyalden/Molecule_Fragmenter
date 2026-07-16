import operator
from rdkit import Chem
from fragmenter import fragmenter
# from fragmenter_utils import draw_mol_with_highlights_and_legend
# from pprint import pprint
# from rdkit.Chem import rdmolops
import csv
import SMARTS
import SMARTS_MARGAN
import SMARTS_MARGAN_GEM
import SMARTS_TDE
import sort_list

# choose which method to use
# SMARTS_LIST = SMARTS.UNIFAC.copy()
# SMARTS_LIST = SMARTS_MARGAN.MARGAN.copy()
SMARTS_LIST = SMARTS_MARGAN_GEM.MARGAN.copy()
# SMARTS_LIST = SMARTS_TDE.UNIFAC.copy()

# get the fragmentation scheme in the format necessary
fragmentation_scheme = {i+1: j[1] for i, j in enumerate(SMARTS_LIST)}
# fragmentation_scheme = {j[0]: j[1] for j in SMARTS_LIST}

def function_to_choose_fragmentation(fragmentations):
    best_fragmentation = None
    best_score = None
    
    # Cache to store SMARTS weights so we don't recalculate them across loops
    weight_cache = {}
    group_names = {
        s: name 
        for name, smarts in SMARTS_LIST 
        for s in (smarts if isinstance(smarts, list) else [smarts])
    }
    for frag in fragmentations:
        group_weights = []
        total_groups = 0
        urea_count = 0
        ac_r_count = 0
        
        for smarts, matches in frag.items():
            num_matches = len(matches)
            total_groups += num_matches
            
            name = group_names[smarts]
            if "urea" in name.lower() or any(sub in name for sub in ["NCON", "NHCON", "NH2CON"]):
                urea_count += num_matches
            if name.startswith("aC-"):
                ac_r_count += num_matches

            if smarts not in weight_cache:
                mol = Chem.MolFromSmarts(smarts)
                weight = 0.0
                if mol:
                    for atom in mol.GetAtoms():
                        weight += atom.GetMass()
                        query = atom.DescribeQuery()
                        if "AtomHCount" in query:
                            num_h = int(query.split('AtomHCount')[1][1])
                            weight += num_h * 1.008
                            
                weight_cache[smarts] = weight
            
            group_weights.extend([weight_cache[smarts]] * num_matches)
            
        group_weights.sort(reverse=True)

        current_score = (urea_count, ac_r_count, group_weights, -total_groups)
                
        if best_score is None or current_score > best_score:
            best_score = current_score
            best_fragmentation = frag
            
    return best_fragmentation

    # fragmentations_descriptors = {}
    # i = 0
    # for fragmentation in fragmentations:
    #     fragmentations_descriptors[i] = [len(fragmentation)]
    #     i += 1
    
    # sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

    # return fragmentations[sorted_fragmentations_dict[0][0]]

calculate_from = "csvs\\margan_groups_check.csv"
write_to = "group_lists\\fragmented_groups.txt"
# sort_to = "group_lists\\fragmented_groups_sorted_corrected_margan.txt"
smiles = []
molecule_ids = []
start_offset = 0

with open(calculate_from, mode='r', encoding='utf-8') as file:
    csv_reader = csv.DictReader(file)
    
    for row in csv_reader:
        if row['smiles'] != "NULL":
            smiles.append(row['smiles'])
            molecule_ids.append(row['id'])

n_noSolution = 0
with open(write_to, "w") as file:
       file.write("")

print('complete algorithm 1')
frg = fragmenter(
    fragmentation_scheme,
    algorithm='complete',
    n_heavy_atoms_cuttoff=30,
    function_to_choose_fragmentation=function_to_choose_fragmentation,
    match_hydrogens=False,
    n_max_fragmentations_to_find=400,
    reject_fragmented_molecules=True,
    )
smiles_length = len(smiles)

# ids removed because of errors/not solving: 19929, 13894, 30784, 22893

for index, smi in enumerate(smiles[start_offset:]):
    # print(smi)
    fragmentation, success, fragmentation_matches = frg.fragment(smi)
    fragmentation_sorted = dict(sorted(fragmentation.items()))
    groups = molecule_ids[index + start_offset] + " "
    for group, number in fragmentation_sorted.items():
        groups +=  str(number) + " " + str(group) + " "
    with open(write_to, "a") as file:
        file.write(groups)
        if index < smiles_length - 1:
            file.write("\n")
    if not success:
        n_noSolution += 1


print ("No Solution: ", n_noSolution)

# sort_list.sort_list(write_to, sort_to)

# fragmentation_scheme = {
#     'CH2' : '[CH2]',
#     'OH' : '[OH]',
#     'CH3' : '[CH3]',
#     'CH2-CH2' : '[CH2][CH2]'
# }
# fragmentation_scheme_order1 = ['CH2-CH2', 'CH3', 'CH2', 'OH']

# print('simple algorithm 1')
# frg = fragmenter(fragmentation_scheme, fragmentation_scheme_order=fragmentation_scheme_order1, algorithm='simple')
# for smi in smiles:
#     fragmentation, success, fragmentation_matches = frg.fragment(smi)
#     print(smi, fragmentation)

# # examples of fragmentation drawing
# for i, SMILES in enumerate(smiles):
#     mol = Chem.MolFromSmiles(SMILES)
#     fragmentation, success, fragmentation_matches = frg.fragment(mol)
#     img = draw_mol_with_highlights_and_legend(mol, fragmentation_matches)
#     img.save(f'simple_example{i+1}.png')

# print()
# print('simple algorithm 2')
# fragmentation_scheme_order2 = ['CH3', 'CH2', 'CH2-CH2', 'OH']
# frg = fragmenter(fragmentation_scheme, fragmentation_scheme_order=fragmentation_scheme_order2, algorithm='simple')
# for smi in smiles:
#     fragmentation, success, fragmentation_matches = frg.fragment(smi)
#     print(smi, fragmentation)

# print()
# print('complete algorithm 2')
# frg = fragmenter(fragmentation_scheme,
#                  algorithm='complete',
#                  n_heavy_atoms_cuttoff=30,
#                  function_to_choose_fragmentation=lambda x: x)
# for smi in smiles:
#     fragmentations, success, fragmentations_matches = frg.fragment_complete(smi)
#     print(smi)
#     pprint(fragmentations)
#     pprint(fragmentations_matches) # some of the fragmentations are the same, but the found fragmentation_matches are different.