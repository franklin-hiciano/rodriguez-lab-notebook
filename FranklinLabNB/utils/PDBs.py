import os
from pathlib import Path
from Bio.Seq import Seq
from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class PDBs:
  def __init__(self, pdbs: list[str]):

    if isinstance(pdbs, list):
      try:
        self.pdbs = [Path(pdb) for pdb in pdbs]
      except TypeError as e:
        raise TypeError("The input must be a list of paths.")
    else:
      raise TypeError("The input must be a list of paths.")

def to_fasta(self, fastas_dir="fastas"):
  os.makedirs(fastas_dir, exist_ok=True)

  pdbParser = PDBParser()
  ppBuilder = PPBuilder()

  fastas = []

  for pdb in self.pdbs:
    structure_chains = pdbParser.get_structure("struc", pdb).get_chains()
    records_of_chains_in_pdb = []

    for chain in structure_chains:
      polypeptides = ppBuilder.build_peptides(chain)
      polypeptide_sequences = [str(polypeptide.get_sequence()) for polypeptide in polypeptides]
      chain_sequence = "".join(polypeptide_sequences)

      chain_record = SeqRecord(
          Seq(chain_sequence),
          id=chain.id
      )

      records_of_chains_in_pdb.append(chain_record)

    fasta_outfile = os.path.join(fastas_dir, f"{Path(pdb).stem}.fasta")

    SeqIO.write(records_of_chains_in_pdb, fasta_outfile, 'fasta')
    fastas.append(fasta_outfile)

  return fastas

PDBs.to_fasta = to_fasta
