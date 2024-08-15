import pytest
from fastapi.testclient import TestClient
from httpx import AsyncClient, ASGITransport
from src.main import app, molecules

client = TestClient(app)


@pytest.fixture(autouse=True)
def setup_molecules():
    """Reset the molecules dictionary before each test."""
    molecules.clear()
    molecules.update({
        1: {"name": "Water", "smiles": "O"},
        2: {"name": "Ethanol", "smiles": "CCO"},
        3: {"name": "Methane", "smiles": "C"},
        4: {"name": "Benzene", "smiles": "c1ccccc1"},
        5: {"name": "Acetic Acid", "smiles": "CC(=O)O"},
        6: {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
    })


@pytest.mark.asyncio
@pytest.mark.parametrize("substructure, expected_count, expected_identifiers", [
    ("C", 5, [2, 3, 4, 5, 6]),
    ("O", 4, [1, 2, 5, 6]),
    ("c1ccccc1", 2, [4, 6]),
    ("N", 0, [])
])
async def test_substructure_search(substructure, expected_count, expected_identifiers):
    async with AsyncClient(transport=ASGITransport(app), base_url="http://test") as ac:
        response = await ac.get(f"/search?substructure={substructure}")
    assert response.status_code == 200
    result = response.json()
    assert len(result) == expected_count
    assert sorted([molecule["identifier"] for molecule in result]) == sorted(expected_identifiers)


@pytest.mark.asyncio
async def test_invalid_substructure_search():
    substructure = "invalid_smiles"
    async with AsyncClient(transport=ASGITransport(app), base_url="http://test") as ac:
        response = await ac.get(f"/search?substructure={substructure}")
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid substructure SMILES:"


def test_add_molecule():
    new_molecule = {
        "identifier": 7,
        "name": "Glucose",
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"
    }
    response = client.post("/molecules", json=new_molecule)
    assert response.status_code == 201
    assert response.json() == new_molecule
    assert molecules[7] == {"name": "Glucose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"}


@pytest.mark.parametrize("molecule_data, expected_status", [
    ({"identifier": 7, "name": "Glucose"}, 422),                # Missing 'smiles'
    ({"identifier": 7, "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"}, 422),  # Missing 'name'
    ({"name": "Glucose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"}, 422),  # Missing 'identifier'
    ({}, 422)                                                   # Missing all fields
])
def test_add_molecule_missing_params(molecule_data, expected_status):
    response = client.post("/molecules", json=molecule_data)
    assert response.status_code == expected_status


def test_add_molecule_duplicate_identifier():
    """ Test for the molecule that has the same identifier as an existing molecule"""
    duplicate_molecule = {
        "identifier": 1,  # Identifier already exists in the initial setup
        "name": "Duplicate Molecule",
        "smiles": "C1CCCCC1"
    }
    response = client.post("/molecules", json=duplicate_molecule)
    assert response.status_code == 400
    assert response.json()["detail"] == "Molecule with this ID already exists."


@pytest.mark.parametrize("update_data, expected_status, expected_response", [
    ({"name": "Updated Methane"}, 200,
     {"identifier": 3, "name": "Updated Methane", "smiles": "C"}),   # Only 'name' provided
    ({"smiles": "C=C"}, 200,
     {"identifier": 3, "name": "Methane", "smiles": "C=C"}),          # Only 'smiles' provided
    ({"name": "Updated Methane", "smiles": "C=C"}, 200,
     {"identifier": 3, "name": "Updated Methane", "smiles": "C=C"}),  # Both provided
    ({}, 422, {"detail": "At least one of 'name' or 'smiles' must be provided"})  # Neither provided
])
def test_update_molecule(update_data, expected_status, expected_response):
    """ Test for update_molecule with various combinations of 'name' and 'smiles' """
    response = client.put("/molecules/3", json=update_data)
    assert response.status_code == expected_status
    assert response.json() == expected_response

    if expected_status == 200:
        assert molecules[3]["name"] == expected_response["name"]
        assert molecules[3]["smiles"] == expected_response["smiles"]


def test_update_molecule_invalid_identifier():
    """ Test for update_molecule with invalid identifier """

    update_data = {
        "name": "Updated Molecule",
        "smiles": "C"
    }
    response = client.put("/molecules/999", json=update_data)  # ID 999 does not exist
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found"


def test_delete_molecule():
    """ Test for delete_molecule with valid identifier """
    response = client.delete("/molecules/4")
    assert response.status_code == 200
    assert "identifier" in response.json()
    assert 4 not in molecules


def test_delete_molecule_invalid_identifier():
    """ Test for delete_molecule with invalid identifier """
    response = client.delete("/molecules/999")  # ID 999 does not exist
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found"


def test_get_molecule():
    """ Test for get_molecule with valid identifier """
    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json() == {"identifier": 1, "name": "Water", "smiles": "O"}


def test_get_molecule_invalid_identifier():
    """ Test for get_molecule with invalid identifier """
    response = client.get("/molecules/999")  # ID 999 does not exist
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found"


def test_upload_molecules():
    """Positive test: Upload valid molecules CSV file."""
    csv_content = "identifier,name,smiles\n7,Glucose,C(C1C(C(C(C(O1)O)O)O)O)O\n8,Fructose,C(C(C1C(C(C(O1)O)O)O)O)O\n"
    response = client.post("/upload", files={"file": ("molecules.csv", csv_content, "text/csv")})
    assert response.status_code == 200
    assert len(molecules) == 8  # 6 initial + 2 uploaded
    assert molecules[7] == {"name": "Glucose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"}


def test_upload_non_csv_file():
    """Negative test: Attempt to upload a non-CSV file."""
    response = client.post("/upload", files={"file": ("molecules.txt", "Some text content", "text/plain")})
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid file format. Only CSV is accepted!"


def test_upload_duplicate_identifier_existing_storage():
    """Negative test: Attempt to upload a CSV file where a molecule identifier already exists in storage."""
    csv_content = "identifier,name,smiles\n1,DuplicateWater,O\n"  # ID 1 already exists
    response = client.post("/upload", files={"file": ("molecules.csv", csv_content, "text/csv")})
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid data: Molecule with identifier 1 already exists."


def test_upload_duplicate_identifier_in_file():
    """Negative test: Attempt to upload a CSV file with duplicate identifiers within the same file."""
    csv_content = "identifier,name,smiles\n7,Glucose,C(C1C(C(C(C(O1)O)O)O)O)O\n7," \
                  "DuplicateGlucose,C(C1C(C(C(C(O1)O)O)O)O)O\n"
    response = client.post("/upload", files={"file": ("molecules.csv", csv_content, "text/csv")})
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid data: File should not contain molecules with same identifier - 7"
