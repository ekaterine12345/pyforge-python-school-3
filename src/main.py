from typing import List
from fastapi import FastAPI, HTTPException, status, UploadFile, File
from rdkit import Chem
from io import StringIO
import csv
from pydantic import BaseModel
from typing import List, Optional


# Pydantic models for request and response bodies
class Molecule(BaseModel):
    identifier: int
    name: str
    smiles: str


class MoleculeUpdate(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None

app = FastAPI()

molecules = [  # For better testing experience there are same values in list
    {"identifier": 1, "name": "Water", "smiles": "O"},
    {"identifier": 2, "name": "Ethanol", "smiles": "CCO"},
    {"identifier": 3, "name": "Methane", "smiles": "C"},
    {"identifier": 4, "name": "Benzene", "smiles": "c1ccccc1"},
    {"identifier": 5, "name": "Acetic Acid", "smiles": "CC(=O)O"},
    {"identifier": 6, "name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
]  # In-memory storage for molecules


@app.post("/molecules", status_code=status.HTTP_201_CREATED, tags=["molecules"], summary="Create a molecule",
          response_description="molecule is created")
def add_molecule(molecule: Molecule):
    """ Add new molecule """
    if any(m['identifier'] == molecule.identifier for m in molecules):
        raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")
    molecules.append(molecule.dict())
    return molecule


@app.get("/molecules",  status_code=status.HTTP_200_OK, tags=["molecules"], summary="Retrieve all molecules",
         response_description="All molecules got returned")
def get_molecules():
    """ Retrieve all molecules """
    return molecules


@app.get("/molecules/{identifier}", response_model=Molecule, status_code=status.HTTP_200_OK, tags=["molecules"],
         summary="Retrieve molecule by identifier", response_description="specified molecule got returned")
def get_molecule(identifier: int):
    """ Retrieve molecule by identifier """
    for mol in molecules:
        if mol.get("identifier") == identifier:
            return mol
    raise HTTPException(status_code=404, detail="molecule is not found")


@app.put("/molecules/{identifier}", response_model=Molecule, tags=["molecules"],
         summary="Update molecule by identifier", response_description="specified molecule got updated")
def update_molecule(identifier: int, molecule_update: MoleculeUpdate):
    """ Update molecule by identifier """
    for index, mol in enumerate(molecules):
        if mol["identifier"] == identifier:
            if molecule_update.name is not None:
                mol["name"] = molecule_update.name
            if molecule_update.smiles is not None:
                mol["smiles"] = molecule_update.smiles
            molecules[index] = mol  # Update the molecule in the list
            return Molecule(**mol)
    raise HTTPException(status_code=404, detail="molecule is not found")


@app.delete("/molecules/{identifier}", response_model=Molecule, tags=["molecules"],
            summary="Delete molecule by identifier", response_description="specified molecule got deleted")
def delete_molecule(identifier: int):
    """ Delete molecule by identifier """
    for index, mol in enumerate(molecules):
        if mol["identifier"] == identifier:
            deleted_mol = molecules.pop(index)
            return deleted_mol
    raise HTTPException(status_code=404, detail="molecule is not found")


@app.get("/search", response_model=List[Molecule], tags=["molecules"],
         summary="substructure search by substructure", response_description="substructure search got performed")
def substructure_search(substructure: str):
    """ Perform  substructure search by specified substructure """
    matched_mols = []            # Initialize an empty list to hold the matched molecules
    substructure_mol = Chem.MolFromSmiles(substructure)  # Convert the substructure SMILES to an RDKit molecule

    if substructure_mol is None:  # Chck if the conversion was successful
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES:")

    for data in molecules:       # To Iterate over the list of molecules
        mol_obj = Chem.MolFromSmiles(data["smiles"])    # Convert the molecule SMILES to an RDKit molecule
        if mol_obj.HasSubstructMatch(substructure_mol):  # Check if the molecule contains the substructure
            matched_mols.append(Molecule(**data))  # If it matches, add the original SMILES to the result list
    return matched_mols                            # Return the list of matched molecules


# Optional
@app.post("/upload",  tags=["molecules"], summary="Reading molecules from CSV file",
          response_description="file read successfully")
async def upload_molecules(file: UploadFile = File(...)):
    """
    The function reads a csv file, reads the molecules and adds them to our storage.
    Upload a CSV file containing molecules with columns: identifier, name, smiles
    """

    if file.content_type != 'text/csv':
        raise HTTPException(status_code=400, detail="Invalid file format. Only CSV is accepted!")

    content = await file.read()
    file_str = StringIO(content.decode())
    reader = csv.DictReader(file_str)

    new_molecules = []
    for row in reader:
        try:
            identifier = int(row["identifier"])
            name = row["name"]
            smiles = row["smiles"]

            # Check for duplicates against the existing storage
            if any(m["identifier"] == identifier for m in molecules):
                raise ValueError(f"Molecule with identifier {identifier} already exists.")

            # Check for duplicate identifiers within the same file
            if any(m["identifier"] == identifier for m in new_molecules):
                raise ValueError(f"File should not contain molecules with same identifier - {identifier}")

            new_molecules.append({"identifier": identifier, "name": name, "smiles": smiles})
        except (ValueError, KeyError) as e:
            raise HTTPException(status_code=400, detail=f"Invalid data: {e}")

    molecules.extend(new_molecules)
    return {"message": "Molecules uploaded successfully.", "molecules": new_molecules}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
