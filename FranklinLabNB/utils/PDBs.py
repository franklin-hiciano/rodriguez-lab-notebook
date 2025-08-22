import os
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder

class PDBs:
  def __init__(self, pdbs):
    if not isinstance(pdbs, list): pdbs = [pdbs]
    self.pdbs = [Path(f) for f in pdbs]
    self.dir = os.path.dirname(os.path.abspath(self.pdbs[0]))

def __init__(self, pdbs):
  if not isinstance(pdbs, list): pdbs = [pdbs]
  self.pdbs = [Path(f) for f in pdbs]
  self.dir = os.path.dirname(os.path.abspath(self.pdbs[0]))

PDBs.to_fasta = to_fasta
