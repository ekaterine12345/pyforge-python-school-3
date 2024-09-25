from os import getenv
from fastapi import FastAPI, HTTPException, status, UploadFile, File
from src.molecule.router import router as molecule_router
from rdkit import Chem
from io import StringIO
import csv
from pydantic import BaseModel
from typing import List, Optional
from fastapi.logger import logger


# Pydantic models for request and response bodies
class Molecule(BaseModel):
    identifier: int
    name: str
    smiles: str


class MoleculeUpdate(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None


app = FastAPI()

molecules = {
    1: {"name": "Water", "smiles": "O"},
    2: {"name": "Ethanol", "smiles": "CCO"},
    3: {"name": "Methane", "smiles": "C"},
    4: {"name": "Benzene", "smiles": "c1ccccc1"},
    5: {"name": "Acetic Acid", "smiles": "CC(=O)O"},
    6: {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
}
# In-memory storage for molecules. For better testing experience there are some values in the list


@app.post("/molecules", status_code=status.HTTP_201_CREATED, tags=["molecules"], summary="Create a molecule",
          response_description="molecule is created")
def add_molecule(molecule: Molecule):
    """ Add new molecule """
    if molecule.identifier in molecules:
        raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")
    molecules[molecule.identifier] = {"name": molecule.name, "smiles": molecule.smiles}
    return molecule


@app.get("/molecules",  status_code=status.HTTP_200_OK, tags=["molecules"], summary="Retrieve all molecules",
         response_description="All molecules got returned")
def get_molecules():
    """ Retrieve all molecules """
    return [{"identifier": k, **v} for k, v in molecules.items()]


@app.get("/molecules/{identifier}", response_model=Molecule, status_code=status.HTTP_200_OK, tags=["molecules"],
         summary="Retrieve molecule by identifier", response_description="specified molecule got returned")
def get_molecule(identifier: int):
    """ Retrieve molecule by identifier """
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": identifier, **molecules[identifier]}


@app.put("/molecules/{identifier}", response_model=Molecule, tags=["molecules"],
         summary="Update molecule by identifier", response_description="specified molecule got updated")
def update_molecule(identifier: int, molecule_update: MoleculeUpdate):
    """ Update molecule by identifier """

    # Validate that at least one field is provided
    if molecule_update.name is None and molecule_update.smiles is None:
        raise HTTPException(status_code=422, detail="At least one of 'name' or 'smiles' must be provided")

    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")

    mol = molecules[identifier]
    if molecule_update.name is not None:
        mol["name"] = molecule_update.name
    if molecule_update.smiles is not None:
        mol["smiles"] = molecule_update.smiles
    molecules[identifier] = mol  # Update the molecule in the dictionary
    return {"identifier": identifier, **mol}


@app.delete("/molecules/{identifier}", response_model=Molecule, tags=["molecules"],
            summary="Delete molecule by identifier", response_description="specified molecule got deleted")
def delete_molecule(identifier: int):
    """ Delete molecule by identifier """
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    deleted_mol = molecules.pop(identifier)
    return {"identifier": identifier, **deleted_mol}


@app.get("/search", response_model=List[Molecule], tags=["molecules"],
         summary="substructure search by substructure", response_description="substructure search got performed")
def substructure_search(substructure: str):
    """ Perform  substructure search by specified substructure """
    matched_mols = []            # Initialize an empty list to hold the matched molecules
    substructure_mol = Chem.MolFromSmiles(substructure)  # Convert the substructure SMILES to an RDKit molecule

    # print(substructure)
    if substructure_mol is None:  # Chck if the conversion was successful
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES:")

    for identifier, data in molecules.items():  # Iterate over the dictionary of molecules
        try:
            mol_obj = Chem.MolFromSmiles(data["smiles"])  # Convert the molecule SMILES to an RDKit molecule
            if mol_obj is None:
                logger.warning(f"Invalid molecule SMILES: {data['smiles']} for identifier {identifier}")
                continue
            if mol_obj.HasSubstructMatch(substructure_mol):  # Check if the molecule contains the substructure
                matched_mols.append(
                    {"identifier": identifier, **data})  # Add matched original SMILES to the result list
        except Exception as e:
            logger.error(f"Error processing molecule with identifier {identifier}: {str(e)}")

    return matched_mols  # Return the list of matched molecules


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

    new_molecules = {}
    for row in reader:
        try:
            identifier = int(row["identifier"])
            name = row["name"]
            smiles = row["smiles"]

            # Check for duplicates against the existing storage
            if identifier in molecules:
                raise ValueError(f"Molecule with identifier {identifier} already exists.")

            # Check for duplicate identifiers within the same file
            if identifier in new_molecules:
                raise ValueError(f"File should not contain molecules with same identifier - {identifier}")

            new_molecules[identifier] = {"name": name, "smiles": smiles}
        except (ValueError, KeyError) as e:
            raise HTTPException(status_code=400, detail=f"Invalid data: {e}")

    molecules.update(new_molecules)
    return {"message": "Molecules uploaded successfully.", "molecules":
            [{"identifier": k, **v} for k, v in new_molecules.items()]}


@app.get("/", tags=["Server"])
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


app.include_router(molecule_router)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
