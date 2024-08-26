from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from src.molecule.request_body import Molecule, MoleculeUpdate
# from molecule.request_body import Molecule, MoleculeUpdate
from typing import List
from src.molecule.dao import MoleculeDao
# from molecule.dao import MoleculeDao
from io import StringIO
import csv
# from src.molecule.schema import DrugResponse, DrugAdd


router = APIRouter(tags=["molecules"])


@router.post("/molecules", status_code=status.HTTP_201_CREATED, summary="Create a molecule",
             response_description="Molecule created")
async def add_molecule(molecule: Molecule):
    """ Add new molecule """
    # existing_molecule = await MoleculeDao.get_molecule_by_id(molecule.identifier)
    # if existing_molecule:
    #     raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")

    new_molecule_id = await MoleculeDao.add_molecule(**molecule.dict())
    return {"identifier": new_molecule_id, **molecule.dict()}


@router.get("/molecules", status_code=status.HTTP_200_OK, summary="Retrieve all molecules",
            response_description="All molecules returned")
async def get_molecules():
    """ Retrieve all molecules """
    molecules = await MoleculeDao.get_all_molecules()
    return molecules


@router.get("/molecules/{identifier}", response_model=Molecule, status_code=status.HTTP_200_OK,
            summary="Retrieve molecule by identifier", response_description="Specified molecule returned")
async def get_molecule(identifier: int):
    """ Retrieve molecule by identifier """
    molecule = await MoleculeDao.get_molecule_by_id(identifier)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@router.put("/molecules/{identifier}", response_model=Molecule, summary="Update molecule by identifier",
            response_description="Specified molecule updated")
async def update_molecule(identifier: int, molecule_update: MoleculeUpdate):
    """ Update molecule by identifier """
    molecule = await MoleculeDao.get_molecule_by_id(identifier)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    update_data = molecule_update.dict(exclude_unset=True)
    await MoleculeDao.update_molecule(identifier, **update_data)
    updated_molecule = await MoleculeDao.get_molecule_by_id(identifier)
    return updated_molecule


@router.delete("/molecules/{identifier}", response_model=Molecule, summary="Delete molecule by identifier",
               response_description="Specified molecule deleted")
async def delete_molecule(identifier: int):
    """ Delete molecule by identifier """
    molecule = await MoleculeDao.get_molecule_by_id(identifier)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    await MoleculeDao.delete_molecule(identifier)
    return molecule


@router.get("/search", response_model=List[Molecule], summary="Substructure search by substructure",
            response_description="Substructure search got performed")
async def substructure_search(substructure: str):
    matched_molecules = await MoleculeDao.substructure_search(substructure)
    if matched_molecules is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES")
    return matched_molecules


@router.post("/upload", summary="Reading molecules from CSV file",
             response_description="File read successfully")
async def upload_molecules(file: UploadFile = File(...)):
    if file.content_type != 'text/csv':
        raise HTTPException(status_code=400, detail="Invalid file format. Only CSV is accepted!")

    content = await file.read()
    file_str = StringIO(content.decode())
    reader = csv.DictReader(file_str)

    new_molecules = []
    for row in reader:
        try:
            # identifier = int(row["identifier"])
            name = row["name"]
            smiles = row["smiles"]

            # existing_molecule = await MoleculeDao.get_molecule_by_id(name=name)
            # if existing_molecule:
            #     raise ValueError(f"Molecule with name {name} already exists.")

            new_molecules.append({"name": name, "smiles": smiles})
        except (ValueError, KeyError) as e:
            raise HTTPException(status_code=400, detail=f"Invalid data: {e}")

    uploaded_molecules = await MoleculeDao.upload_molecules(new_molecules)
    return {"message": "Molecules uploaded successfully.", "molecules": uploaded_molecules}
