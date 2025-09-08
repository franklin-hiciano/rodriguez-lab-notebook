import sys
import os
from pathlib import Path
import subprocess
import glob
import pandas as pd
import shutil
import requests
import errno
from Bio import SeqIO
from urllib.request import urlretrieve
from urllib.error import HTTPError, URLError
from FranklinLabNB.utils.PDBs import PDBs

class SabDabRetriever:
  def __init__(self, sabdab_summary_file_download_link, sabdab_summary_file_path):
    self.sabdab_downloader_python_script = self._get_sabdab_downloader_python_script()
    self.sabdab_summary_file = self._download_sabdab_summary_file(link=sabdab_summary_file_download_link, outpath=sabdab_summary_file_path)

def _ensure_path_ends_with(self, path: str, end: str):
  if path.endswith(end):
    pass
  else:
    path += end

  return path

def _ensure_path_exists(self, path: str):
  if os.path.exists(path):
    pass
  else:
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

  return path

def _protect_path_from_being_overwritten_if_exists_already(self, path: str, overwrite_anyway=False):
  if os.path.exists(path):
    if overwrite_anyway:
      if os.path.isdir(path):
        shutil.rmtree(path)
        os.makedirs(path)
      else:
        os.remove(path)
    else:
      raise FileExistsError(f"{path} already exists. Please delete or enable overwrite_anyway.")
  else:
    if len(path.split(Path(path).stem)[-1]) < 2:
      os.makedirs(path)
    else:
      pass

  return path

def _get_sabdab_downloader_python_script(self, outpath: str=None):
  DOWNLOAD_LINK_FOR_SABDAB_DOWNLOADER_PYTHON_SCRIPT = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/"

  if outpath is not None:
    outpath = self._ensure_path_ends_with(outpath, '.py')
  else: 
    outpath ='sabdab_downloader.py'

  print("Getting SabDab downloader python script...")
  urlretrieve(DOWNLOAD_LINK_FOR_SABDAB_DOWNLOADER_PYTHON_SCRIPT, outpath)
  print(f"Successfully downloaded SabDab downloader python script to {outpath}.")

  return outpath

def _download_sabdab_summary_file(self, link: str, outpath: str):
    outpath = self._protect_path_from_being_overwritten_if_exists_already(outpath, overwrite_anyway=True)
    outpath = self._ensure_path_ends_with(outpath, '.tsv')
    try:
        print("Downloading SabDab summary file...")
        urlretrieve(link, outpath)
        print(f"Successfully downloaded SabDab summary file to {outpath}.")
    except HTTPError as err:
        if err.code == 404:
            print(f"Your SabDab download link is invalid or expired. Error: {err}")
        else:
            print(f"HTTP error while downloading SabDab summary file: {err}")
    except URLError as err:
        print(f"Network error while downloading SabDab summary file: {err}")

    return outpath

def get_structures(self, outpath: str):
  outpath = self._protect_path_from_being_overwritten_if_exists_already(outpath, overwrite_anyway=True)
  print(f"Downloadiing SabDab structures to {outpath}...")
  command = [
    sys.executable,
    self.sabdab_downloader_python_script,
    '--summary_file', self.sabdab_summary_file,
    '--output_path', outpath,
    '--chothia_pdb'
  ]
  subprocess.run(command, check=True)
  print(f"Successfully downloaded SabDab structures to {outpath}.")

def structures_to_sequences(self, structures_dir: str, outpath: str):
    import os, glob, re
    from pathlib import Path
    import pandas as pd
    from Bio import SeqIO

    self._ensure_path_ends_with(outpath, ".tsv")
    self._protect_path_from_being_overwritten_if_exists_already(outpath, overwrite_anyway=True)

    # emit FASTAs
    fasta_dir = "fastas"
    pdbs = glob.glob(os.path.join(structures_dir, "**", "*.pdb"), recursive=True)
    PDBs(pdbs).to_fasta(fasta_dir)

    # (pdb, chain) -> sequence
    seq_map = {}
    for fa in glob.glob(os.path.join(fasta_dir, "**", "*.fasta"), recursive=True):
        pdb = Path(fa).stem[:4].lower()
        for rec in SeqIO.parse(fa, "fasta"):
            m = re.search(r"([A-Za-z0-9])$", rec.id)
            if m:
                seq_map[(pdb, m.group(1).upper())] = str(rec.seq)

    # load and replace
    df = pd.read_csv(self.sabdab_summary_file, sep="\t")
    if "pdb" not in df.columns:
        raise KeyError("Table needs a 'pdb' column.")
    pdb_col = df["pdb"].astype(str).str.lower()

    targets = [c for c in ("Hchain", "Lchain", "antigen_chain") if c in df.columns]
    for c in targets:
        df[c + "_id"] = df[c].astype(str)

    def replace_cell(pdb, cell):
        if pd.isna(cell): return cell
        toks = re.split(r"[,\s]+", str(cell).strip())
        toks = [t for t in toks if t]
        return ",".join([seq_map.get((pdb, t.upper()), t) for t in toks]) if toks else ""

    for c in targets:
        df[c] = [replace_cell(p, v) for p, v in zip(pdb_col, df[c])]

    df.to_csv(outpath, sep="\t", index=False)
    print(f"Wrote {outpath}")

SabDabRetriever._ensure_path_ends_with = _ensure_path_ends_with
SabDabRetriever._ensure_path_exists = _ensure_path_exists
SabDabRetriever._protect_path_from_being_overwritten_if_exists_already = _protect_path_from_being_overwritten_if_exists_already
SabDabRetriever._get_sabdab_downloader_python_script = _get_sabdab_downloader_python_script
SabDabRetriever._download_sabdab_summary_file = _download_sabdab_summary_file
SabDabRetriever.get_structures = get_structures
SabDabRetriever.structures_to_sequences = structures_to_sequences
