import os
import sys
from datetime import datetime
import contextlib
import subprocess
from urllib.request import urlretrieve
from urllib.request import HTTPError
from pathlib import Path
import glob
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder

class SabDabRetriever:
  def __init__(self, summary_file_download_link):
    self.timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    self.summary_file = f"sabdab_summary_{self.timestamp}.tsv"
    urlretrieve(summary_file_download_link, self.summary_file)

def get_structures(self):
  urlretrieve("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/", "sabdab_downloader.py")
  os.makedirs((outpath := f"sabdab_structures_{self.timestamp}"))
  subprocess.run([sys.executable, "sabdab_downloader.py", "--summary_file", self.summary_file, "--chothia_pdb", "--output_path", outpath], check=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
  return outpath

def get_sequences_and_structures(self, structures_dir=None):
  df = pd.read_csv(self.summary_file, sep="\t")
  with contextlib.chdir(structures_dir or self.get_structures()):
    print(os.getcwd())
    pdbprs = PDBParser()
    ppb = PPBuilder()
    def pdbs():
      pdb_id = lambda pdb_file: Path(pdb_file).stem
      chains = lambda pdb_file: pdbprs.get_structure("ab-ag_complex", pdb_file).get_chains()
      amino_acids = lambda chn: ("".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chn)), print(chn.id))[0]

      return {pdb_id(f): {chn.id: amino_acids(chn) for chn in chains(f)} for f in glob.glob("*/**/*.pdb", recursive=True)}

    def insert_seqs_into_table():
      for pdb_id, chns in pdbs().items():
        for chn_id, chn_seq in chns.items():
          pdb_rows = df[df['pdb']==pdb_id]
          chn_row = pdb_rows[(pdb_rows==chn_id[0]).any(axis=1)]
          if chn_row.empty: continue
          chn_col_name = pdb_rows.columns[(pdb_rows==chn_id).any(axis=0)][0]+'_seq'
          df.loc[chn_row.index[0], chn_col_name] = chn_seq
    insert_seqs_into_table()

  df.to_csv((out := Path(self.summary_file).stem+'_with_seqs.tsv'), sep='\t', index=False)
  return out

SabDabRetriever.get_structures = get_structures
SabDabRetriever.get_sequences_and_structures = get_sequences_and_structures
