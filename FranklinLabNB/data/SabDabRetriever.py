from pathlib import Path
import subprocess
import glob
import pandas as pd
import shutil
from urllib.request import urlretrieve
from FranklinLabNB.utils.PDBs import PDBs
from Bio import SeqIO

class SabDabRetriever:
  def __init__(self, summary_file_download_link, summary_table_outpath="sabdab_summary.tsv"):
    urlretrieve(summary_file_download_link, summary_table_outpath)
    self.summary_table = summary_table_outpath

def get_structures(self, outdir="sabdab_structures"):
  shutil.rmtree(outdir)
  os.makedirs(outdir, exist_ok=True)
  urlretrieve("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/", "sabdab_downloader.py")
  subprocess.run([
      sys.executable,
      "sabdab_downloader.py",
      '--summary_file', self.summary_table,
      '--output_path', outdir,
      '--chothia_pdb'
  ], check=True)

def convert_antibody_structures_to_fastas(self, structures_dir, outdir):
  antibody_structures = glob.glob(os.path.join(structures_dir, "**/*.pdb"), recursive=True)
  print(antibody_structures)
  PDBs(antibody_structures).to_fasta(fastas_dir=outdir)

def insert_antibody_sequences_into_sabdab_summary_table(self, fasta_dir, filled_summary_table_outpath):
    summary_table = pd.read_csv(self.summary_table, sep="\t")
    print("Columns:", summary_table.columns.tolist())

    # collect all FASTA files
    fastas_with_chain_sequences = glob.glob(os.path.join(fasta_dir, "**/*.fasta"), recursive=True)
    fastas_with_chain_sequences += glob.glob(os.path.join(fasta_dir, "**/*.fa"), recursive=True)

    def get_sequence_of_structure_chain(pdb_id, chain_id):
        fasta_file = next((f for f in fastas_with_chain_sequences if Path(f).stem == pdb_id), None)
        if fasta_file is None:
            print(f"No fasta found for {pdb_id}")
            return None

        for record in SeqIO.parse(fasta_file, "fasta"):
            print(str(record.seq))
            if record.id == chain_id:
              return str(record.seq)

        print(f"Sequence not found for PDB {pdb_id}, chain {chain_id}")
        return None

    def insert_antibody_chain_sequence_into_row(row):
        pdb_id = row["pdb"]
        row["Hchain_seq"] = get_sequence_of_structure_chain(pdb_id, row["Hchain"])
        row["Lchain_seq"] = get_sequence_of_structure_chain(pdb_id, row["Lchain"])
        row["antigen_chain_seq"] = get_sequence_of_structure_chain(pdb_id, row["antigen_chain"])
        return row

    summary_table = summary_table.apply(insert_antibody_chain_sequence_into_row, axis=1)
    summary_table.to_csv(filled_summary_table_outpath, sep="\t", index=False)

SabDabRetriever.get_structures = get_structures
SabDabRetriever.convert_antibody_structures_to_fastas = convert_antibody_structures_to_fastas
SabDabRetriever.insert_antibody_sequences_into_sabdab_summary_table = insert_antibody_sequences_into_sabdab_summary_table
