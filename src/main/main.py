from fastapi import FastAPI, Depends, HTTPException, UploadFile, File
from sqlalchemy.orm import Session
from src.mydb.database import Base, SessionLocal, engine, database
from src.mydb.models import Molecule
from pydantic import BaseModel
from rdkit import Chem
from os import getenv

app = FastAPI()

# Create the database tables
Base.metadata.create_all(bind=engine)

# Dependency to get the database session


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


class MoleculeCreate(BaseModel):
    identifier: str
    smiles: str


class MoleculeUpdate(BaseModel):
    smiles: str


class SubstructureSearchRequest(BaseModel):
    substructure: str


def validate_smiles(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    return molecule


@app.on_event("startup")
async def startup():
    await database.connect()


@app.on_event("shutdown")
async def shutdown():
    await database.disconnect()


@app.post("/molecule/", response_model=dict)
def add_molecule(molecule: MoleculeCreate, db: Session = Depends(get_db)):
    if db.query(Molecule).filter(Molecule.identifier == molecule.identifier).first():
        raise HTTPException(
            status_code=400, detail="Identifier already exists")
    validate_smiles(molecule.smiles)
    db_molecule = Molecule(
        identifier=molecule.identifier, smiles=molecule.smiles)
    db.add(db_molecule)
    db.commit()
    return {"message": "Molecule added successfully"}


@app.get("/molecule/{identifier}", response_model=dict)
def get_molecule(identifier: str, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": molecule.identifier, "smiles": molecule.smiles}


@app.put("/molecule/{identifier}", response_model=dict)
def update_molecule(identifier: str, molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    validate_smiles(molecule.smiles)
    db_molecule.smiles = molecule.smiles
    db.commit()
    return {"message": "Molecule updated successfully"}


@app.delete("/molecule/{identifier}", response_model=dict)
def delete_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(db_molecule)
    db.commit()
    return {"message": "Molecule deleted successfully"}


@app.get("/molecules/", response_model=list)
def list_molecules(db: Session = Depends(get_db)):
    molecules = db.query(Molecule).all()
    return [{"identifier": molecule.identifier, "smiles": molecule.smiles} for molecule in molecules]


@app.post("/substructure_search/", response_model=list)
def substructure_search(request: SubstructureSearchRequest, db: Session = Depends(get_db)):
    substructure = validate_smiles(request.substructure)
    result = []
    molecules = db.query(Molecule).all()
    for molecule in molecules:
        mol = Chem.MolFromSmiles(molecule.smiles)
        if mol.HasSubstructMatch(substructure):
            result.append({"identifier": molecule.identifier,
                          "smiles": molecule.smiles})
    return result


@app.post("/upload_molecules/", response_model=dict)
async def upload_molecules(file: UploadFile = File(...), db: Session = Depends(get_db)):
    content = await file.read()
    content_str = content.decode("utf-8")
    lines = content_str.splitlines()
    for line in lines:
        identifier, smiles = line.split()
        if db.query(Molecule).filter(Molecule.identifier == identifier).first():
            raise HTTPException(
                status_code=400,
                detail=f"Identifier {identifier} already exists"
            )
        validate_smiles(smiles)
        db_molecule = Molecule(identifier=identifier, smiles=smiles)
        db.add(db_molecule)
    db.commit()
    return {"message": "Molecules uploaded successfully"}


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
