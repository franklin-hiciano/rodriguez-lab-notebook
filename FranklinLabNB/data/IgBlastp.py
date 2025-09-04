import os
import subprocess
import glob
import tarfile
from urllib.request import urlretrieve
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from FranklinLabNB.utils.Fastas import Fastas

class IgBlastp:
  def __init__(self):
    pass

def install(self):
  print("Installing IGBlastp...")

  def download_binaries():
    urlretrieve("https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz", "igblast.tar.gz")
    tarfile.open("igblast.tar.gz", "r:gz").extractall()
  download_binaries()

  def download_human_germline_V_genes_database():

    with open("human_germline_V_genes.fasta", 'w') as nucleotide_germline_V_genes_fasta:

      for chain in ['H', 'K', 'L']:
        try:
          process = subprocess.run(['download_germline_set',
                          'Homo sapiens',
                          f'IG{chain}',
                          '-f', 'MULTI-IGBLAST',
                          '-p', f'human_IG{chain}',],
                          check=True,
                          text=True,
                          stderr=subprocess.STDOUT,
                          stdout=subprocess.PIPE,
                          bufsize=1,
                          env=os.environ.copy())

          for line in process.stdout:
            print(line, end="")

        except subprocess.CalledProcessError as e:
          print("stdout:\n", e.stdout)
          print("stderr:\n", e.stderr)
          raise 


        with open(f"human_IG{chain}V.fasta") as f:
          current_chain_nucleotide_germline_V_genes_fasta = f.read()

        nucleotide_germline_V_genes_fasta.write(
            current_chain_nucleotide_germline_V_genes_fasta
        )

    Fastas(["human_germline_V_genes.fasta"]).translate()

    try:
      makeblastdb_binary = glob.glob("*/**/makeblastdb")[0]
      process = subprocess.run([
        makeblastdb_binary,
        "-parse_seqids", 
        "-dbtype", "prot",
        "-in", "human_germline_V_genes_aa.fasta",
        "-out", "human_gl_V"], check=True,
                               text=True,
                               stderr=subprocess.STDOUT,
                               stdout=subprocess.PIPE,
                               bufsize=1,
                               env=os.environ.copy())
      for line in process.stdout:
        print(line, end="")
    except subprocess.CalledProcessError as e:
      print("stdout:\n", e.stdout)
      print("stderr:\n", e.stderr)
      raise 


  download_human_germline_V_genes_database()

  print("Installed IGBlastp successfully.")
  return glob.glob("../**/ncbi-igblast*")[0]

def run(self, fastas: list[os.PathLike]):
  def run():
    IGBlast_dir = glob.glob("../**/ncbi-igblast*")[0]
    os.environ["IGDATA"] = IGBlast_dir

    IGBlastp_binary = glob.glob("*/**/igblastp")[0]

    for fasta in fastas:
      try:
        process = subprocess.run([
          IGBlastp_binary,
          '-query', fasta,
          '-germline_db_V', 'human_gl_V',
          '-organism', 'human'],
                                check=True,
                                text=True,
                                stderr=subprocess.STDOUT,
                                stdout=subprocess.PIPE,
                                bufsize=1,
                                env=os.environ.copy())

        for line in process.stdout:
          print(line, end="")

      except subprocess.CalledProcessError as e:
        print("stdout:\n", e.stdout)
        print("stderr:\n", e.stderr)
        raise 


  try:
    run()
  except:
    print("IGBlastp is not installed, installing and trying to run again...")
    self.install()
    run()

IgBlastp.install = install
IgBlastp.run = run
