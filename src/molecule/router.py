from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from src.molecule.request_body import Molecule, MoleculeUpdate
from typing import List, Dict
from src.molecule.dao import MoleculeDao
from io import StringIO
import csv
import logging
import os
from typing import Iterator, Optional
import redis
import json

# Connect to Redis
redis_client = redis.Redis(host='redis', port=6379, db=0)

router = APIRouter(tags=["Molecules"])
# logging.basicConfig(level=logging.INFO, filename="logs.log", filemode="w")
logging.basicConfig(
    level=logging.DEBUG,  # Set to DEBUG to capture all log levels
    filename="logs.log",
    filemode="w",  # Overwrites the file every time the program runs
    format="%(asctime)s - %(levelname)s - %(message)s"
)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))


@router.post("/mol/molecules", status_code=status.HTTP_201_CREATED, summary="Create a molecule",
             response_description="Molecule created")
async def add_molecule(molecule: Molecule):
    """ Add new molecule """
    logging.info("\nInside the add_molecule() function in router.py file")
    new_molecule_id = await MoleculeDao.add_molecule(**molecule.dict())
    logging.info(f"Adding a new molecule {molecule} with id = {new_molecule_id}")
    return {"identifier": new_molecule_id, **molecule.dict()}


@router.get("/mol/molecules", status_code=status.HTTP_200_OK, summary="Retrieve all molecules",
            response_model=Dict, response_description="All molecules returned")
async def get_molecules(limit: Optional[int] = None):
    """ Retrieve all molecules """
    logging.info("\n")
    logging.info("Inside the get_molecules() function in router.py file")

    cache_key = f"get_molecules:{limit}"
    cached_result = get_cached_result(cache_key)
    logging.debug(f"cache_key={cache_key}; cached_result={cached_result}")

    if cached_result:
        return {"source": "cache", "data": cached_result}

    # print(os.getcwd())
    molecules = []
    async for molecule in MoleculeDao.get_all_molecules(limit=limit):
        molecules.append(molecule)

    logging.info(f"Molecules are {molecules}")
    # Convert Molecule objects to a dictionary with identifier as the key
    molecules_dicts = {molecule.identifier: molecule.to_dict() for molecule in molecules}

    logging.info(f"Molecules for limit = {limit} are molecules_dicts = {molecules_dicts}")

    set_cache(cache_key, molecules_dicts, expiration=300)
    return {"source": "database", "data": molecules_dicts}


@router.get("/mol/molecules/{identifier}", response_model=Dict, status_code=status.HTTP_200_OK,
            summary="Retrieve molecule by identifier", response_description="Specified molecule returned")
async def get_molecule(identifier: int):
    """ Retrieve molecule by identifier """
    logging.info("\nInside the get_molecule() function in router.py file")

    cache_key = f"get_molecule:{identifier}"
    cached_result = get_cached_result(cache_key)
    logging.debug(f"cache_key={cache_key}; cached_result={cached_result}")

    if cached_result:
        return {"source": "cache", "data": cached_result}

    molecule = await MoleculeDao.get_molecule_by_id(identifier)

    if not molecule:
        logging.debug(f"Molecule with identifier = {identifier} was not found")
        raise HTTPException(status_code=404, detail="Molecule not found")

    logging.info(f"Molecule is {molecule}")
    set_cache(cache_key, molecule.to_dict(), expiration=120)
    return {"source": "database", "data": molecule.to_dict()}


@router.put("/mol/molecules/{identifier}", response_model=Molecule, summary="Update molecule by identifier",
            response_description="Specified molecule updated")
async def update_molecule(identifier: int, molecule_update: MoleculeUpdate):
    """ Update molecule by identifier """
    logging.info("\nInside the update_molecule() function in router.py file")

    molecule = await MoleculeDao.get_molecule_by_id(identifier)
    if not molecule:
        logging.debug(f"Molecule with identifier = {identifier} was not found")
        raise HTTPException(status_code=404, detail="Molecule not found")

    logging.info(f"Provided molecule is molecule_update = {molecule_update}")

    update_data = molecule_update.dict(exclude_unset=True)
    await MoleculeDao.update_molecule(identifier, **update_data)
    updated_molecule = await MoleculeDao.get_molecule_by_id(identifier)
    logging.info(f"Updated Molecule is {updated_molecule}")
    return updated_molecule


@router.delete("/mol/molecules/{identifier}", response_model=Molecule, summary="Delete molecule by identifier",
               response_description="Specified molecule deleted")
async def delete_molecule(identifier: int):
    """ Delete molecule by identifier """
    logging.info("\nInside the delete_molecule() function in router.py file")
    molecule = await MoleculeDao.get_molecule_by_id(identifier)
    if not molecule:
        logging.debug(f"Molecule with identifier = {identifier} was not found")
        raise HTTPException(status_code=404, detail="Molecule not found")

    await MoleculeDao.delete_molecule(identifier)
    logging.info(f"Molecule {molecule} with identifier = {identifier} got deleted")
    return molecule


@router.get("/mol/search", response_model=Dict, summary="Substructure search by substructure",
            response_description="Substructure search got performed")
async def substructure_search(substructure: str):
    logging.info("\nInside the substructure_search() function in router.py file")
    cache_key = f"search:{substructure}"
    cached_result = get_cached_result(cache_key)

    logging.debug(f"cache_key={cache_key}; cached_result={cached_result}")

    if cached_result is not None:
        return {"source": "cache", "data": cached_result}

    matched_molecules = await MoleculeDao.substructure_search(substructure)
    if matched_molecules is None:
        logging.debug(f"Molecules with substructure = {substructure} is not correct")
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES")

    if not matched_molecules:
        logging.debug(f"Molecules with substructure = {substructure} was not found")
        set_cache(cache_key, {}, expiration=120)
        return {"source": "database", "data": []}

    logging.info(f"Matched molecules for substructure = {substructure} are {matched_molecules}")

    # Convert Molecule objects to a dictionary with identifier as the key
    matched_molecules_dicts = {molecule.identifier: molecule.to_dict() for molecule in matched_molecules}

    logging.info(f"Matched molecules for substructure = {substructure} are "
                 f"matched_molecules_dicts = {matched_molecules_dicts}")

    set_cache(cache_key, matched_molecules_dicts, expiration=300)
    return {"source": "database", "data": list(matched_molecules_dicts.values())}


@router.post("/mol/upload", summary="Reading molecules from CSV file",
             response_description="File read successfully")
async def upload_molecules(file: UploadFile = File(...)):
    logging.info("\nInside the upload_molecules() function in router.py file")
    if file.content_type != 'text/csv':
        logging.debug(f"Invalid content_type={file.content_type}; Only text/csv file formats are accepted! ")
        raise HTTPException(status_code=400, detail="Invalid file format. Only CSV is accepted!")

    content = await file.read()
    file_str = StringIO(content.decode())
    reader = csv.DictReader(file_str)

    new_molecules = []
    for row in reader:
        try:
            name = row["name"]
            smiles = row["smiles"]

            # Check for duplicate names within the same file
            if any(m["name"] == name for m in new_molecules):
                logging.debug(f"File should not contain molecules with same name - {name}")
                raise ValueError(f"File should not contain molecules with same name - {name} ")

            # Check for duplicate smiles within the same file
            if any(m["smiles"] == name for m in new_molecules):
                logging.debug(f"File should not contain molecules with same smiles - {smiles}")
                raise ValueError(f"File should not contain molecules with same smiles - {smiles} ")

            existing_molecule_by_name = await MoleculeDao.get_molecule_by_name(name=name)
            if existing_molecule_by_name:
                logging.debug(f"Molecule with name {name} already exists.")
                raise ValueError(f"Molecule with name {name} already exists.")

            existing_molecule_by_smiles = await MoleculeDao.get_molecule_by_smiles(smiles=smiles)
            if existing_molecule_by_smiles:
                logging.debug(f"Molecule with smiles {smiles} already exists.")
                raise ValueError(f"Molecule with smiles {smiles} already exists.")

            new_molecules.append({"name": name, "smiles": smiles})
        except (ValueError, KeyError) as e:
            logging.debug(f"Error occurred during reading a file. error = {e} ")
            raise HTTPException(status_code=400, detail=f"Invalid data: {e}")

    uploaded_molecules = await MoleculeDao.upload_molecules(new_molecules)
    logging.info(f"Uploaded molecules successfully {uploaded_molecules}")
    return {"message": "Molecules uploaded successfully.", "molecules": uploaded_molecules}
