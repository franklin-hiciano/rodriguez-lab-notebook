import subprocess
import textwrap
import os
import sys
import gdown
from urllib.request import urlretrieve
import pandas as pd

class FrankenReferenceGenome:
  """
  Desc: Downloads the Franken Reference genome
  Platforms: Minerva High-Performance Compute at Icahn School of Medicine at Mount Sinai
  """
  def __init__(self):
    self.outdir = "/sc/arion/scratch/hiciaf01/imputation/franken_reference_genome"
    self.fasta_file = os.path.join(self.outdir, "reference.fasta")
    self.fasta_url = "http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta"
    self.fasta_index_file = os.path.join(self.outdir, "reference.fasta.fai")
    self.fasta_index_url = "https://github.com/Watson-IG/immune_receptor_genomics/blob/main/240520/reference.fasta.fai"

  def run(self):
    self._ensure_dir(directory=self.outdir)
    self.download()
    # subprocess.run(["ls", "-l", self.outdir], check=True)


  def _ensure_dir(self, directory):
    if not os.path.exists(directory):
      os.makedirs(directory, exist_ok=True)

def download_fasta_index(self):
  print("Downloading reference fasta index...")
  urlretrieve(self.fasta_index_url, self.fasta_index_file)
  print("Downloaded reference fasta index successfully...")

def download_fasta(self):
  print("Downloading reference fasta...")
  urlretrieve(self.fasta_url, self.fasta_file)
  print("Downloaded reference fasta successfully...")

def download(self):
  self.download_fasta_index()
  self.download_fasta()

FrankenReferenceGenome.download_fasta_index = download_fasta_index
FrankenReferenceGenome.download_fasta = download_fasta
FrankenReferenceGenome.download = download
