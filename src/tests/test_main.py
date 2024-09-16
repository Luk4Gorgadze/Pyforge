import pytest
from fastapi.testclient import TestClient
from main.main import app, molecules_db

client = TestClient(app)


@pytest.fixture(autouse=True)
def setup_and_teardown():
    # Setup: Clear the molecules_db before each test
    molecules_db.clear()
    yield
    # Teardown: Clear the molecules_db after each test
    molecules_db.clear()


def test_add_molecule():
    response = client.post(
        "/molecule/",
        json={
            "identifier": "mol1",
            "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule added successfully"}


def test_get_molecule():
    # Add a molecule to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})

    # Retrieve the molecule
    response = client.get("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}


def test_update_molecule():
    # Add a molecule to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})

    # Update the molecule
    response = client.put(
        "/molecule/mol1", json={"identifier": "mol1", "smiles": "CCN"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}

    # Retrieve the updated molecule
    response = client.get("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCN"}


def test_delete_molecule():
    # Add a molecule to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})

    # Delete the molecule
    response = client.delete("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully"}

    # Try to retrieve the deleted molecule
    response = client.get("/molecule/mol1")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_list_molecules():
    # Add molecules to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    client.post("/molecule/", json={"identifier": "mol2", "smiles": "CCN"})

    # List all molecules
    response = client.get("/molecules/")
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 2
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result


def test_substructure_search():
    # Add molecules to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    client.post("/molecule/", json={"identifier": "mol2", "smiles": "CCN"})
    client.post("/molecule/", json={"identifier": "mol3", "smiles": "CCC"})

    # Perform substructure search
    response = client.post(
        "/substructure_search/",
        json={
            "substructure": "CC"})
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 3
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result
    assert {"identifier": "mol3", "smiles": "CCC"} in result


def test_substructure_search_no_match():
    # Add molecules to the database
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    client.post("/molecule/", json={"identifier": "mol2", "smiles": "CCN"})

    # Perform substructure search with no match
    response = client.post(
        "/substructure_search/",
        json={
            "substructure": "NNN"})
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 0


def test_upload_molecules():
    # Create a mock file with molecule data
    file_content = "mol1 CCO\nmol2 CCN\n"
    file = {"file": ("molecules.txt", file_content)}

    # Upload the molecules
    response = client.post("/upload_molecules/", files=file)
    assert response.status_code == 200
    assert response.json() == {"message": "Molecules uploaded successfully"}

    # Verify the molecules were added
    response = client.get("/molecules/")
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 2
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result
