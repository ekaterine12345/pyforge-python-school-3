from celery import Celery
import logging
from src.molecule.dao import MoleculeDao
from rdkit import Chem
from src.molecule.celery_worker import celery
from src.molecule.celery_worker import logging
from src.molecule.redis_utils import redis_client, get_cached_result, set_cache
import asyncio
from fastapi import APIRouter, Depends, HTTPException

# logging.basicConfig(
#     level=logging.DEBUG,  # Set to DEBUG to capture all log levels
#     filename="logs.log",
#     filemode="w",  # Overwrites the file every time the program runs
#     format="%(asctime)s - %(levelname)s - %(message)s"
# )


@celery.task(name='src.molecule.celery_tasks.substructure_search_task')  # Explicitly naming the task
def substructure_search_task(substructure: str):
    logging.info(f"in celery_tasks.py file - Started substructure search task for substructure: {substructure}")
    cache_key = f"search:{substructure}"
    cached_result = get_cached_result(cache_key)

    logging.debug(f"cache_key={cache_key}; cached_result={cached_result}")

    if cached_result is not None:
        return {"source": "cache", "data": cached_result}

    substructure_mol = Chem.MolFromSmiles(substructure)
    if substructure_mol is None:
        logging.error(f"Invalid substructure SMILES: {substructure}")
        raise ValueError(f"Invalid substructure SMILES {substructure}")

    # to run the asynchronous DAO method within the current event loop
    loop = asyncio.get_event_loop()

    # Check if loop is running
    if loop.is_running():
        matched_molecules = loop.run_until_complete(MoleculeDao.substructure_search(substructure))
    else:
        matched_molecules = loop.run_until_complete(MoleculeDao.substructure_search(substructure))

    if matched_molecules is None:
        logging.debug(f"Molecules with substructure = {substructure} is not correct")
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES")

    if not matched_molecules:
        logging.debug(f"Molecules with substructure = {substructure} was not found")
        set_cache(cache_key, {}, expiration=120)
        return {"source": "database", "data": []}

    logging.info(f"Found {len(matched_molecules)} molecules matching the substructure: {substructure}")
    logging.info(f"Matched molecules for substructure = {substructure} are {matched_molecules}")

    matched_molecules_dicts = {molecule.identifier: molecule.to_dict() for molecule in matched_molecules}

    logging.info(f"Matched molecules for substructure = {substructure} are "
                 f"matched_molecules_dicts = {matched_molecules_dicts}")

    set_cache(cache_key, matched_molecules_dicts, expiration=300)
    return {"source": "database", "data": list(matched_molecules_dicts.values())}
