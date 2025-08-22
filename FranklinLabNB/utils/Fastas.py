import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Fastas:
  def __init__(self, fastas: dict[os.PathLike]):
    self.fastas = [Path(f) for f in fastas]
    self.dir = os.path.dirname(os.path.abspath(self.fastas[0]))

def translate(self):
  for nt_fa in self.fastas:
    try:
      with open(os.path.join(self.dir, nt_fa.stem+'_aa'+nt_fa.suffix) ,'w+') as aa_fa:
        for rec in SeqIO.parse(nt_fa, 'fasta'):
          SeqIO.write(SeqRecord(seq=rec.seq.translate(to_stop=True), id=rec.id), aa_fa, format='fasta')
    except Exception as e:
      print(f"Could not translate {nt_fa}: {e}")

Fastas.translate = translate
