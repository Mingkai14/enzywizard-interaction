from __future__ import annotations


from ..utils.logging_utils import Logger


from rdkit import Chem, DataStructs

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")



import math

from ..resources.substrate_resources import *



# Validation
def is_valid_smiles(smiles: str) -> bool:
    try:
        if not isinstance(smiles, str):
            return False
        smiles = smiles.strip()
        if not smiles:
            return False
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False

def is_valid_mol_2d(mol: Chem.Mol, logger: Logger) -> bool:
    if mol is None:
        logger.print("[ERROR] Input Mol(2D) is None.")
        return False

    try:
        if not isinstance(mol, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(2D).")
            return False

        if mol.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(2D) contains no atoms.")
            return False

        Chem.SanitizeMol(mol)
        return True
    except Exception:
        logger.print("[ERROR] Input Mol(2D) is invalid or failed sanitization.")
        return False


def is_valid_mol_h(mol_h: Chem.Mol, logger: Logger) -> bool:
    if mol_h is None:
        logger.print("[ERROR] Input Mol(H) is None.")
        return False

    try:
        if not isinstance(mol_h, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(H).")
            return False

        if mol_h.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(H) contains no atoms.")
            return False

        Chem.SanitizeMol(mol_h)

        has_h = any(atom.GetAtomicNum() == 1 for atom in mol_h.GetAtoms())
        if not has_h:
            logger.print("[ERROR] Input Mol(H) does not contain explicit hydrogen atoms.")
            return False

        return True
    except Exception:
        logger.print("[ERROR] Input Mol(H) is invalid or failed sanitization.")
        return False


def is_valid_conf_3d(conf: Chem.Conformer, logger: Logger) -> bool:
    if conf is None:
        logger.print("[ERROR] Input conformer is None.")
        return False

    try:
        if not isinstance(conf, Chem.Conformer):
            logger.print("[ERROR] Input object is not an RDKit Conformer.")
            return False

        if not conf.Is3D():
            logger.print("[ERROR] Input conformer is not 3D.")
            return False

        if conf.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input conformer contains no atoms.")
            return False

        for atom_idx in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            if any(math.isnan(v) or math.isinf(v) for v in [pos.x, pos.y, pos.z]):
                logger.print("[ERROR] Input conformer contains invalid 3D coordinates.")
                return False

        return True
    except Exception:
        logger.print("[ERROR] Input conformer(3D) is invalid.")
        return False


def is_valid_mol_3d(mol_3d: Chem.Mol, logger: Logger) -> bool:
    if mol_3d is None:
        logger.print("[ERROR] Input Mol(3D) is None.")
        return False

    try:
        if not isinstance(mol_3d, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(3D).")
            return False

        if mol_3d.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(3D) contains no atoms.")
            return False

        Chem.SanitizeMol(mol_3d)

        if mol_3d.GetNumConformers() <= 0:
            logger.print("[ERROR] Input Mol(3D) contains no conformer.")
            return False

        conf = mol_3d.GetConformer()
        if not is_valid_conf_3d(conf, logger):
            return False

        if conf.GetNumAtoms() != mol_3d.GetNumAtoms():
            logger.print("[ERROR] Atom count mismatch between Mol(3D) and conformer.")
            return False

        return True
    except Exception:
        logger.print("[ERROR] Input Mol(3D) is invalid.")
        return False
































