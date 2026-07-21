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

def create_fragmentation_scheme_order():
    """
    Creates an order to prioritize fragmentations
    """
    scheme_descriptors = []
    smarts_to_name = {
        s: name 
        for name, smarts in SMARTS_LIST 
        for s in (smarts if isinstance(smarts, list) else [smarts])
    }

    for group_id, smarts in fragmentation_scheme.items():
        if isinstance(smarts, list):
            smarts = smarts[0]
        mol_SMARTS = fragmenter.Chem.MolFromSmarts(smarts)
        group_name = smarts_to_name[smarts]
        weight = 0.0
        if mol_SMARTS:
            for atom in mol_SMARTS.GetAtoms():
                weight += atom.GetMass()
                query = atom.DescribeQuery()
                if "AtomHCount" in query:
                    h_match = int(query.split("AtomHCount")[1][1])
                else:
                    h_match = 0
                if h_match:
                    weight += h_match * 1.008

        is_urea = 1 if any(sub in group_name for sub in ["NCON", "NHCON", "NH2CON"]) else 0
        
        is_ac_r = 1 if group_name.startswith("aC-") else 0
        
        scheme_descriptors.append((group_id, is_urea, is_ac_r, weight, len(smarts)))
        
    # Priority: Ureas > aC-R > Heaviest Weight > Longest SMARTS string
    fragmentation_scheme_order = [
        i for i, is_urea, is_ac_r, weight, length in sorted(
            scheme_descriptors, 
            key=lambda x: (x[1], x[2], x[3], x[4]), 
            reverse=True
        )
    ]
    return fragmentation_scheme_order
fragmentation_scheme_order_margan = create_fragmentation_scheme_order()

def function_to_choose_fragmentation(fragmentations):
    """
        Selects a fragmentation from a list of fragmentation for Marrero-Gani method
        checks to ensure fragmentations have used non-split groups where applicable
        selects the fragmentation with the heaviest single group
    """
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
            if any(sub in name for sub in ["NCON", "NHCON", "NH2CON"]):
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

CUSTOM_PROPERTY_MATCH_FUNCTIONS = {
    "IsInRing": (
        lambda query_atom: True if "AtomInNRings" in query_atom.DescribeQuery() or query_atom.GetIsAromatic() else False,
        lambda atom: atom.IsInRing()
    ),
}

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
    fragmentation_scheme_order=fragmentation_scheme_order_margan,
    properties_to_match=["GetTotalNumHs", "GetFormalCharge", "IsInRing"],
    custom_property_match_functions=CUSTOM_PROPERTY_MATCH_FUNCTIONS
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