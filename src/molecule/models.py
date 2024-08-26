from sqlalchemy.orm import Mapped
from src.database import Base, str_uniq, int_pk

class Molecule(Base):
    identifier: Mapped[int_pk]
    name: Mapped[str_uniq]
    smiles: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.identifier}, "
            f"name={self.name!r},"
            f"smiles={self.smiles!r})"
        )

    def __repr__(self):
        return str(self)
