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
    urlretrieve("download_link", self.summary_file)

def get_structures(self):
  urlretrieve("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/", "sabdab_downloader.py")
  subprocess.run([sys.executable, "sabdab_downloader.py", "--summary_file", self.summary_file, "--chothia_pdb", "--output_path", (dir := f"sabdab_structures_{self.timestamp}")], check=True)
  return dir

def get_sequences_and_structures(self):
  os.chdir(self.get_structures())
  df = pd.read_csv(self.summary_file, sep="\t")
  def pdbs():
    pdb_id = lambda pdb_file: Path(pdb_file).stem
    chains = lambda pdb_file: PDBParser().get_structure("ab-ag_complex", pdb_file).get_chains()
    amino_acids = lambda chn: "".join(str(pp.get_sequence) for pp in PPBuilder().build_peptides(chn))
    return {pdb_id(f): {chn.id: amino_acids(chn) for chn in chains(f)} for f in glob.glob("*/**/*.pdb", recursive=True)}

  def insert_seqs_into_table():
    for pdb_id, chns in pdbs().items():
      for chn_id, chn_seq in chns.items:
        pdb_rows = df[df['pdb']==pdb_id]
        chn_row = pdb_rows[(pdb_rows==chn_id).any(axis=1)]
        chn_col_name = pdb_rows.columns[(pdb_rows==chn_id).any(axis=0)][0]+'_seq'
        df.loc[chn_row.index[0], chn_col_name] = chn_seq
  insert_seqs_into_table()

  df.to_csv((out := Path(self.summary_file).stem+'_with_seqs.tsv'), sep='\t', index=False)
  return out

SabDabRetriever.get_structures = get_structures
SabDabRetriever.get_sequences_and_structures = get_sequences_and_structures
