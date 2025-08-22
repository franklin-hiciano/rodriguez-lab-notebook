import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Fastas:
  def __init__(self, fastas: dict[os.PathLike]):
    self.fastas = [Path(f) for f in fastas]
    self.dir = os.path.dirname(os.path.abspath(self.fastas[0]))

def translate(self):
  for i, nt_fa in enumerate(self.fastas):
    try:
      aa_fa = os.path.join(self.dir, nt_fa.stem+'_aa'+nt_fa.suffix)
      with open(aa_fa, 'w+') as f:
        for rec in SeqIO.parse(nt_fa, 'fasta'):
          SeqIO.write(SeqRecord(seq=rec.seq.translate(to_stop=True), id=rec.id), f, format='fasta')
        self.fastas[i] = aa_fa
    except Exception as e:
      print(f"Could not translate {nt_fa}: {e}")
  return self

def merge(self, out: os.PathLike):
  with open(out,'w') as out:
    for fa in self.fastas:
      out.write(open(fa).read())
  return self

Fastas.translate = translate
Fastas.merge = merge
