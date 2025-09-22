import subprocess
import textwrap
import os
import sys
import gdown
from urllib.request import urlretrieve
import pandas as pd

class ShortReads:
  """
  Desc: Downloads CRAM files of DNA short-reads.
  Platforms: Minerva High-Performance Compute at Icahn School of Medicine at Mount Sinai
  """

  def run(self):
    self.download_sample_names()
    self.download_index_files()
    self.remove_comments_from_index_files()
    self.extract_cram_files_download_links_from_index_files()
    self.download_cram_files()

def download_sample_names(self):
  file_id = "1ArzWBQqjSpl_KNQFgPW-at0l0BTZnvvu"
  url = f"https://drive.google.com/uc?id={file_id}"
  output = "sample_names.txt"  # name to save as
  gdown.download(url, output, quiet=False)

def get_short_read_urls(self):
  outpath = "urls.tsv"
  urlretrieve("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index")
  df = pd.read_csv('urls.tsv', sep='\t', comment='#')

def download_index_files(self,
                         index_file_1_path="1000G_2504_high_coverage.sequence.index",
                         index_file_2_path="1000G_698_related_high_coverage.sequence.index"
                         ):
  index_file_1_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
  index_file_2_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index"

  urlretrieve(index_file_1_url, index_file_1_path)
  urlretrieve(index_file_2_url, index_file_2_path)

def remove_comments_from_index_files(self,
    index_file_1_path="1000G_2504_high_coverage.sequence.index",
    index_file_2_path="1000G_698_related_high_coverage.sequence.index"
):
    for path in [index_file_1_path, index_file_2_path]:
        cleaned_lines = []
        with open(path, "r") as f:
            for line in f:
                if line.startswith("##"):
                    continue  # skip metadata
                if line.startswith("#"):
                    line = line[1:]  # strip leading "#" from header
                cleaned_lines.append(line)

        # Write out a new file alongside the original
        out_path = os.path.splitext(path)[0] + ".cleaned.index"
        with open(out_path, "w") as out:
            out.writelines(cleaned_lines)

        print(f"Cleaned file written to: {out_path}")

def extract_cram_files_download_links_from_index_files(self,
                         sample_names_file_path="sample_names.txt",
                         index_file_1_path="1000G_2504_high_coverage.sequence.cleaned.index",
                         index_file_2_path="1000G_698_related_high_coverage.sequence.cleaned.index",
                         download_links_table_outpath="cram_download_links.tsv"
                         ):
  sample_names = open(sample_names_file_path, "r").read().splitlines()
  df_filtered_all = pd.DataFrame()
  for df_index_file in [index_file_1_path, index_file_2_path]:
    df = pd.read_csv(df_index_file, sep="\t", comment="#")
    df_filtered = df[df["SAMPLE_NAME"].isin(sample_names)][["SAMPLE_NAME", "ENA_FILE_PATH"]]
    print(f"Found {len(df_filtered)} sample names: {df_filtered["SAMPLE_NAME"]}")
    df_filtered_all = pd.concat([df_filtered_all, df_filtered], ignore_index=True)
  df_filtered_all.to_csv(download_links_table_outpath, sep='\t', index=False)

def download_cram_files(self, outdir="/sc/arion/scratch/hiciaf01/short_reads",download_links_table_path="cram_download_links.tsv", paths_outpath="paths.txt"):
  download_links = pd.read_csv(download_links_table_path, sep="\t")
  download_links["LOCAL_PATH"] = download_links["ENA_FILE_PATH"].str.replace(
        r"^ftp://ftp\.sra\.ebi\.ac\.uk", "", regex=True
    )
  download_links["LOCAL_PATH"].to_csv(paths_outpath, index=False, header=False)

  body = textwrap.dedent("""
#!/bin/bash
# Usage: ./aspera_download.sh /sc/arion/scratch/hiciaf01/short_reads

set -euo pipefail

# Paths (from the modulefile you showed and your working example)
SIF="/hpc/packages/minerva-centos7/aspera-connect/4.2.4/src/ubuntu_latest.sif"
ASPERA_KEY="/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh"
ASPERA_USER="era-fasp@fasp.sra.ebi.ac.uk"

# Input arg
URL_PATH="/vol1/run/ERR323/ERR3239480/NA12718.final.cram"

# Local output directory (you can adjust this)
OUTDIR="$1"
mkdir -p "$OUTDIR"

# Extract filename from URL_PATH
FN=$(basename "$URL_PATH")

# Run inside container with Singularity exec
function download {
singularity exec -B /hpc/packages:/hpc/packages "$SIF" \
  ascp \
    -i "$ASPERA_KEY" \
    -P33001 -O33001 \
    -T -L- \
    "${ASPERA_USER}:${URL_PATH}" \
    "${OUTDIR}/${FN}"
}
module load singularity/3.6.4   # or just 'module load singularity' if that’s your default

# Make sure these are defined once:
SIF="/hpc/packages/minerva-centos7/aspera-connect/4.2.4/src/ubuntu_latest.sif"
ASPERA_KEY="/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh"

# Export env into the container (this is the key bit you were missing)
export SINGULARITYENV_PATH="/hpc/packages/minerva-centos7/aspera-connect/4.2.4/aspera/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
export SINGULARITYENV_LD_LIBRARY_PATH="/hpc/packages/minerva-centos7/aspera-connect/4.2.4/aspera/lib:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/lib"

function download_batch {
  xargs -a paths.txt -P4 -I{} \
    singularity exec -B /hpc/packages:/hpc/packages "$SIF" \
      ascp -i "$ASPERA_KEY" -P33001 -O33001 -k2 -QT --overwrite=diff \
      "era-fasp@fasp.sra.ebi.ac.uk:{}" "$OUTDIR/"
}

download_batch
""")

  bsub_cmd = (
    f"bsub -P acc_oscarlr -q interactive -W 12:00 -n 6 "
    f"-J genotype_SVs -Is /bin/bash -lc "
    f"\"nohup bash aspera_download.sh {shlex.quote(outdir)} > genotype_SVs.$LSB_JOBID.out 2>&1 & exit\""
)


  out = subprocess.run(bsub_cmd, shell=True,
                         text=True, capture_output=True, check=False)
  print(out.stdout)
  print(out.stderr)

ShortReads.download_sample_names = download_sample_names
ShortReads.get_short_read_urls = get_short_read_urls
ShortReads.download_index_files = download_index_files
ShortReads.remove_comments_from_index_files = remove_comments_from_index_files
ShortReads.extract_cram_files_download_links_from_index_files = extract_cram_files_download_links_from_index_files
ShortReads.download_cram_files = download_cram_files
