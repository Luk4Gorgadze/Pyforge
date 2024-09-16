from os import getenv
from fastapi import FastAPI, HTTPException, UploadFile, File
from pydantic import BaseModel
from typing import Dict
from rdkit import Chem

app = FastAPI()

molecules_db: Dict[str, str] = {}


class Molecule(BaseModel):
    identifier: str
    smiles: str


class SubstructureSearchRequest(BaseModel):
    substructure: str


def validate_smiles(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    return molecule


@app.post("/molecule/")
def add_molecule(molecule: Molecule):
    if molecule.identifier in molecules_db:
        raise HTTPException(
            status_code=400,
            detail="Identifier already exists")
    validate_smiles(molecule.smiles)
    molecules_db[molecule.identifier] = molecule.smiles
    return {"message": "Molecule added successfully"}


@app.get("/molecule/{identifier}")
def get_molecule(identifier: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": identifier, "smiles": molecules_db[identifier]}


@app.put("/molecule/{identifier}")
def update_molecule(identifier: str, molecule: Molecule):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found")
    validate_smiles(molecule.smiles)
    molecules_db[identifier] = molecule.smiles
    return {"message": "Molecule updated successfully"}


@app.delete("/molecule/{identifier}")
def delete_molecule(identifier: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules_db[identifier]
    return {"message": "Molecule deleted successfully"}


@app.get("/molecules/")
def list_molecules():
    return [{"identifier": identifier, "smiles": smiles}
            for identifier, smiles in molecules_db.items()]


@app.post("/substructure_search/")
def substructure_search(request: SubstructureSearchRequest):
    substructure = validate_smiles(request.substructure)
    result = []
    for identifier, smiles in molecules_db.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule.HasSubstructMatch(substructure):
            result.append({"identifier": identifier, "smiles": smiles})
    return result


@app.post("/upload_molecules/")
async def upload_molecules(file: UploadFile = File(...)):
    content = await file.read()
    content_str = content.decode("utf-8")
    lines = content_str.splitlines()
    for line in lines:
        identifier, smiles = line.split()
        if identifier in molecules_db:
            raise HTTPException(
                status_code=400,
                detail=f"Identifier {identifier} already exists"
            )
        validate_smiles(smiles)
        molecules_db[identifier] = smiles
    return {"message": "Molecules uploaded successfully"}


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
