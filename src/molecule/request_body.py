from pydantic import BaseModel
from typing import List, Optional


class Molecule(BaseModel):
    name: str
    smiles: str


class MoleculeUpdate(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None
