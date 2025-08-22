from urllib.request import urlretrieve
import tarfile
import subprocess
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

class IgBlastp:
  def __init__(self):
    pass

def istl(self):
  def _dl():
    urlretrieve("https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz", "igblast.tar.gz")
    tarfile.open("igblast.tar.gz", "r:gz").extractall()
  _dl()

  def _dl_db():
    with open("nt.fa",'w+') as ntides, open("aa.fa",'w') as aas:
      for chn in ["H", "K", "L"]:
        subprocess.run(["download_germline_set", "Homo sapiens", f"IG{chn}", "-f", "MULTI-IGBLAST", "-p", f"human_IG{chn}"], check=True)
        ntides.write(open(f"human_IG{chn}V.fasta").read())
      for rec in SeqIO.parse('nt.fa', 'fasta'):
        SeqIO.write(SeqRecord(seq=rec.seq.translate(to_stop=True), id=rec.id), aas, format='fasta')
      subprocess.run([glob.glob("*/**/makeblastdb")[0], "-parse_seqids", "-dbtype", "prot", "-in", "aa.fa", "-out", "human_gl_V"], check=True)
  _dl_db()
  return glob.glob("../**/ncbi-igblast*")[0]

def run(self, fastas: list[os.PathLike]):
  def run():
    os.environ["IGDATA"] = glob.glob("../**/ncbi-igblast*")[0]
    p = subprocess.run([glob.glob("*/**/igblastp")[0], "-query", "test.fasta", "-germline_db_V", "human_gl_V", "-organism", "human"], check=True, text=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, bufsize=1, env=os.environ.copy())
    for line in p.stdout:
      print(line, end="")
  try:
    run()
  except:
    self.install()
    run()

IgBlastp.install = install
IgBlastp.run = run
