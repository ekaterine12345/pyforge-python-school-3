from rdkit import Chem


def substructure_search(mols, mol):
    matched_mols = []  # Initialize an empty list to hold the matched molecules
    substructure_mol = Chem.MolFromSmiles(mol)  # Convert the substructure SMILES to an RDKit molecule

    for each_mol in mols:  # To Iterate over the list of molecules
        mol_obj = Chem.MolFromSmiles(each_mol)     # Let's convert the molecule SMILES to an RDKit molecule

        if mol_obj.HasSubstructMatch(substructure_mol):   # Check if the molecule contains the substructure
            matched_mols.append(each_mol)           # If it matches, add the original SMILES string to the result list

    return matched_mols                             # Return the list of matched molecules


print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))
# Output:  ['c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O']
