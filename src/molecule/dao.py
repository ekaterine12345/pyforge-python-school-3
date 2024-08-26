from sqlalchemy import delete, update
from sqlalchemy.future import select
from src.molecule.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from rdkit import Chem
# from molecule.models import Molecule
# from database import async_session_maker
# from dao.base import BaseDAO


class MoleculeDao(BaseDAO):
    model = Molecule

    @classmethod
    async def get_all_molecules(cls):
        async with async_session_maker() as session:
            query = select(cls.model)
            molecules = await session.execute(query)
            return molecules.scalars().all()

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
                await session.flush()
                new_molecule_id = new_molecule.identifier
                await session.commit()
                return new_molecule_id

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
                await session.execute(query)
                await session.commit()
                return molecule_id

    @classmethod
    async def substructure_search(cls, substructure: str):
        async with async_session_maker() as session:
            query = select(cls.model)
            result = await session.execute(query)
            molecules = result.scalars().all()

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
