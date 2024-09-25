from sqlalchemy import delete, update
from sqlalchemy.future import select
from src.molecule.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from rdkit import Chem
# from molecule.models import Molecule
# from database import async_session_maker
# from dao.base import BaseDAO
# import psycopg2
from sqlalchemy.exc import IntegrityError
from fastapi import HTTPException
import logging
from typing import Iterator, Optional


class MoleculeDao(BaseDAO):
    model = Molecule


    @classmethod
    async def get_all_molecules(cls, limit: Optional[int] = None) -> Iterator[Molecule]:
        async with async_session_maker() as session:
            logging.info(f"In the dao.py file - get_all_molecules function; limit = {limit} [2]")
            logging.info(f"in the dao.py file - cls.model = {cls.model}")
            query = select(cls.model)
            logging.info(f"In the dao.py file - query - {query}")
            if limit:
                query = query.limit(limit)
            molecules = await session.execute(query)
            logging.info(f"in dao.py file - molecules - {molecules}")
            print(cls)
            # return molecules.scalars().all()
            for molecule in molecules.scalars():
                logging.info(f"Molecule {molecule}")
                yield molecule

    @classmethod
    async def get_molecule_by_id(cls, molecule_id: int):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(identifier=molecule_id)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()

            if not molecule_info:   # If molecule is not found, return None
                return None

            return molecule_info

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                new_molecule = cls.model(**molecule_data)
                session.add(new_molecule)
                try:
                    await session.flush()
                    new_molecule_id = new_molecule.identifier
                    await session.commit()
                    return new_molecule_id
                except IntegrityError as e:
                    await session.rollback()
                    if "unique constraint" in str(e.orig):
                        logging.debug(f"Molecule {molecule_data} with this  name or smiles already exists [From dao.py]")
                        raise HTTPException(status_code=400, detail="Molecule with this name or smiles already exists.")
                    logging.debug(f"An error occurred while adding the molecule. error = {e} [FROM dao.py]")
                    raise HTTPException(status_code=500, detail="An error occurred while adding the molecule.")

    @classmethod
    async def delete_molecule(cls, molecule_id: int):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(identifier=molecule_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                if not molecule_to_delete:
                    return None

                await session.execute(delete(cls.model).filter_by(identifier=molecule_id))
                await session.commit()
                return molecule_id

    @classmethod
    async def update_molecule(cls, molecule_id: int, **molecule_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                query = (
                    update(cls.model)
                    .where(cls.model.identifier == molecule_id)
                    .values(**molecule_data)
                    .execution_options(synchronize_session="fetch")
                )
                try:
                    await session.execute(query)
                    await session.commit()
                    return molecule_id
                except IntegrityError as e:
                    await session.rollback()
                    if "unique constraint" in str(e.orig):
                        logging.debug(
                            f"Molecule {molecule_data} with this  name or smiles already exists [From dao.py]")
                        raise HTTPException(status_code=400, detail="Molecule with this name or smiles already exists.")
                    logging.debug(f"An error occurred while adding the molecule. error = {e} [FROM dao.py]")
                    raise HTTPException(status_code=500, detail="An error occurred while updating the molecule.")

    @classmethod
    async def substructure_search(cls, substructure: str):
        async with async_session_maker() as session:
            query = select(cls.model)
            result = await session.execute(query)
            molecules = result.scalars().all()

            # print(query)
            substructure_mol = Chem.MolFromSmiles(substructure)
            if substructure_mol is None:
                return None  # Handle invalid SMILES case

            matched_molecules = []
            for molecule in molecules:
                mol_obj = Chem.MolFromSmiles(molecule.smiles)
                if mol_obj and mol_obj.HasSubstructMatch(substructure_mol):
                    matched_molecules.append(molecule)

            return matched_molecules

    @classmethod
    async def upload_molecules(cls, molecules_list: list[dict]):
        async with async_session_maker() as session:
            async with session.begin():
                new_molecules = [cls.model(**data) for data in molecules_list]
                session.add_all(new_molecules)
                await session.commit()
                return new_molecules

    @classmethod
    async def get_molecule_by_name(cls, name: str):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(name=name)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()
            print("name = ", name, molecule_info)

            if not molecule_info:   # If molecule is not found, return None
                return None

            return molecule_info

    @classmethod
    async def get_molecule_by_smiles(cls, smiles: str):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(smiles=smiles)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()

            if not molecule_info:   # If molecule is not found, return None
                return None

            return molecule_info

