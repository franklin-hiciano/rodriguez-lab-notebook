import os
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder

class PDBs:
  def __init__(self, pdbs):
    if not isinstance(pdbs, list): pdbs = [pdbs]
    self.pdbs = [Path(f) for f in pdbs]
    self.dir = os.path.dirname(os.path.abspath(self.pdbs[0]))

def to_fasta(self):
  pdbprs = PDBParser()
  ppb = PPBuilder()
  fastas = []
  for pdb in self.pdbs:
    chains = pdbprs.get_structure("p", pdb).get_chains()
    fasta = os.path.join(pdb.parent, pdb.stem+".fasta")
    with open(fasta, 'w') as f:
      for c in chains:
        amino_acids = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(c))
        f.write(f""">{pdb.stem}|{c.id}\n{amino_acids}\n""")
    fastas.append(fasta)
  return fastas

PDBs.to_fasta = to_fasta
