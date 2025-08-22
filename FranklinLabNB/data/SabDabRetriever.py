import os
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
  def __init__(self):
    self.summary_file_path = None
    self.structures_dir = None

def get_summary_file(self, download_link, destination=None):
  destination = destination or f"sabdab_summary_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.tsv"
  try:
    urlretrieve(download_link, destination)
    self.summary_file_path = destination
  except Exception as e:
    print(e)
    if "404" in str(e): sys.exit("[ERROR] Try a new SabDab download link, yours may be expired.")
    else: sys.exit()

def get_structures(self, destination_dir=None):
  destination_dir = destination_dir or f"sabdab_complexes_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}"
  urlretrieve("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/", "sabdab_downloader.py")
  os.makedirs(destination_dir, exist_ok=True)
  cmd = [
    sys.executable, "-u", "sabdab_downloader.py",
    "--summary_file", self.summary_file_path,
    "--chothia_pdb",
    "--output_path", destination_dir
  ]
  try:
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    subprocess.run(cmd, env=env, check=True)
    self.structures_dir = destination_dir
  except Exception as e:
    print(f"[ERROR] Error while downloading structures from SabDab: {e}")
    sys.exit(1)

def extract_sequences_from_structures(self, destination=None):
  destination = destination or Path(self.summary_file_path).stem+".tsv"
  df_summary = pd.read_csv(self.summary_file_path, sep="\t")

  with contextlib.chdir(self.structures_dir):
    pdb_files = glob.glob("*/**/*.pdb", recursive=True)

    pdbprs = PDBParser(QUIET=True)
    ppb = PPBuilder()

    def get_pdb_id(pdb_file):
      return Path(pdb_file).stem
    def get_chains(pdb_file):
      return pdbprs.get_structure("ab-ag_complex", pdb_file).get_chains()
    def get_amino_acids(chain):
      return "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))

    pdbs = {get_pdb_id(f): {chain.id: get_amino_acids(chain) for chain in get_chains(f)}
            for f in pdb_files}

    def insert_into_table(pdb_id, chain_id, chain_seq):

      pdb_rows = df_summary[df_summary['pdb'] == pdb_id]
      chain_row = pdb_rows[(pdb_rows == chain_id).any(axis=1)]
      chain_col_name = pdb_rows.columns[(pdb_rows == chain_id).any(axis=0)][0] + '_sequence'
      df_summary.loc[chain_row.index[0], chain_col_name] = chain_seq

    [insert_into_table(pdb_id, chain_id, chain_seq)
    for pdb_id, chains in pdbs.items()
    for chain_id, chain_seq in chains.items()]

  df_summary.to_csv(destination, sep='\t', index=False)

SabDabRetriever.get_summary_file = get_summary_file
SabDabRetriever.get_structures = get_structures
SabDabRetriever.extract_sequences_from_structures = extract_sequences_from_structures
