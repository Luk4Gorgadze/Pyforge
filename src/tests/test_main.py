import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.main.main import app, get_db
from src.mydb.database import Base, DATABASE_URL
from src.mydb.models import Molecule

SQLALCHEMY_DATABASE_URL = DATABASE_URL

engine = create_engine(SQLALCHEMY_DATABASE_URL)
TestingSessionLocal = sessionmaker(
    autocommit=False, autoflush=False, bind=engine)

Base.metadata.create_all(bind=engine)

client = TestClient(app)


@pytest.fixture(scope="module")
def db():
    Base.metadata.create_all(bind=engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()
        Base.metadata.drop_all(bind=engine)


@pytest.fixture(autouse=True)
def setup_and_teardown(db):
    yield
    db.query(Molecule).delete()
    db.commit()


def override_get_db():
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db


def test_add_molecule():
    response = client.post(
        "/molecule/",
        json={
            "identifier": "mol1",
            "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule added successfully"}


def test_get_molecule():
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    response = client.get("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}


def test_update_molecule():
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    response = client.put(
        "/molecule/mol1", json={"identifier": "mol1", "smiles": "CCN"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}
    response = client.get("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCN"}


def test_delete_molecule():
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    response = client.delete("/molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully"}
    response = client.get("/molecule/mol1")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_list_molecules():
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    client.post("/molecule/", json={"identifier": "mol2", "smiles": "CCN"})
    response = client.get("/molecules/")
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 2
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result


def test_substructure_search():
    client.post("/molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    client.post("/molecule/", json={"identifier": "mol2", "smiles": "CCN"})
    response = client.post("/substructure_search/",
                           json={"substructure": "CC"})
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 2
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result


def test_upload_molecules():
    file_content = "mol1 CCO\nmol2 CCN"
    response = client.post(
        "/upload_molecules/",
        files={"file": ("molecules.txt", file_content, "text/plain")}
    )
    assert response.status_code == 200
    assert response.json() == {"message": "Molecules uploaded successfully"}
    response = client.get("/molecules/")
    assert response.status_code == 200
    result = response.json()
    assert len(result) == 2
    assert {"identifier": "mol1", "smiles": "CCO"} in result
    assert {"identifier": "mol2", "smiles": "CCN"} in result


def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()
