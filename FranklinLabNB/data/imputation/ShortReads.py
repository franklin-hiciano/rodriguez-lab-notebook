import subprocess
import textwrap
import os
import sys
import gdown
from urllib.request import urlretrieve
import pandas as pd
import re
from pathlib import Path
import os, shlex, subprocess, textwrap, pandas as pd

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
    self.launch_aspera_batch_to_download_cram_files()


  def __init__(self):
    # ---INTERMEDIARY DATA, i.e. handling download links---

    self.SAMPLES_FILE = "/sc/arion/work/hiciaf01/short_reads/HGSVC_2024_sample_ids.txt"
    self.INDEX1_FILE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
    self.INDEX2_FILE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index"
    self.INDEX1_FILE = "1000G_2504_high_coverage.sequence.index"
    self.INDEX2_FILE = "1000G_698_related_high_coverage.sequence.index"

    # ---DOWNLOADING THE SHORT READS---
    # inputs / outputs
    self.DOWNLOAD_LINKS_TABLE = "/sc/arion/work/hiciaf01/short_reads/cram_download_links.tsv"  # must contain ENA_FILE_PATH column
    self.PATHS_TXT           = "/sc/arion/work/hiciaf01/short_reads/aspera_paths.txt"
    self.ASPERA_SCRIPT       = "/sc/arion/work/hiciaf01/short_reads/aspera_download.sh"
    self.OUTDIR              = "/sc/arion/scratch/hiciaf01/short_reads"

    # singularity / aspera env
    self.MODULE_SINGULARITY        = "singularity/3.6.4"  # or "singularity"
    self.SIF                       = "/hpc/packages/minerva-centos7/aspera-connect/4.2.4/src/ubuntu_latest.sif"
    self.ASPERA_KEY                = "/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh"
    self.SINGULARITYENV_PATH       = "/hpc/packages/minerva-centos7/aspera-connect/4.2.4/aspera/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
    self.SINGULARITYENV_LD_LIBRARY_PATH = "/hpc/packages/minerva-centos7/aspera-connect/4.2.4/aspera/lib:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/lib"

    # batching
    self.XARGS_PROCS = 4  # parallel downloads

    # LSF
    self.LSF_ACCOUNT  = "acc_oscarlr"
    self.LSF_QUEUE    = "interactive"
    self.LSF_WALLTIME = "12:00"
    self.LSF_CORES    = 6
    self.LSF_JOB_NAME = "genotype_SVs"

def download_sample_names(self):
  file_id = "1ArzWBQqjSpl_KNQFgPW-at0l0BTZnvvu"
  url = f"https://drive.google.com/uc?id={file_id}"
  output = self.SAMPLES_FILE  # name to save as
  gdown.download(url, output=output, quiet=False)

def get_short_read_urls(self):
  outpath = "urls.tsv"
  urlretrieve("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index")
  df = pd.read_csv('urls.tsv', sep='\t', comment='#')

def download_index_files(self):

  urlretrieve(self.INDEX1_FILE_URL, self.INDEX1_FILE)
  urlretrieve(self.INDEX2_FILE_URL, self.INDEX2_FILE)

def remove_comments_from_index_files(self):
    for path in [self.INDEX1_FILE, self.INDEX2_FILE]:
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

def _load_samples(self):
    """Reads newline-separated sample names from self.SAMPLES_FILE."""
    return [
        ln.strip()
        for ln in Path(self.SAMPLES_FILE).read_text(errors="ignore").splitlines()
        if ln.strip()
    ]

def _scan_index_text_for_links(self, index_text, wanted):
    """
    Scan raw index text for ftp URLs and derive sample name from the filename.
    Returns dict {sample: url} for samples in 'wanted'.
    """
    url_re = re.compile(r'ftp://[^\s\t]+', re.IGNORECASE)
    found = {}
    for line in index_text.splitlines():
        m = url_re.search(line)
        if not m:
            continue
        url = m.group(0).rstrip("/")
        fname = url.split("/")[-1]       # e.g., HG00405.final.cram
        sample = fname.split(".")[0]     # -> HG00405
        if sample in wanted and sample not in found:
            found[sample] = url
    return found

def _ensure_tsv_header(self, out_path):
    """
    Ensure a TSV exists with the header row:
    SAMPLE_NAME\tENA_FILE_PATH
    """
    p = Path(out_path)
    if not p.exists() or p.stat().st_size == 0:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("SAMPLE_NAME\tENA_FILE_PATH\n")
        print(f"[init] created TSV with header → {out_path}")

def _already_written_samples(self, out_path):
    """
    Return a set of sample names already present in the output TSV (SAMPLE<TAB>URL).
    If file doesn't exist, returns empty set.
    """
    p = Path(out_path)
    if not p.exists():
        return set()
    existing = set()
    for line in p.read_text(errors="ignore").splitlines():
        if not line.strip():
            continue
        sample = line.split("\t", 1)[0].strip()
        if sample:
            existing.add(sample)
    return existing

def extract_links_from_1000G_2504_high_coverage(self):
    """
    Uses:
      self.SAMPLES_FILE
      self.INDEX1_FILE
      self.DOWNLOAD_LINKS_TABLE  (headered TSV with columns: SAMPLE_NAME, ENA_FILE_PATH)

    Appends new rows and reports how many were added.
    """
    # read sample names and scan index text
    samples = self._load_samples()
    wanted = set(samples)
    index_text = Path(self.INDEX1_FILE).read_text(errors="ignore")
    mapping = self._scan_index_text_for_links(index_text, wanted)

    # ensure TSV exists with header, then append only new samples
    self._ensure_tsv_header(self.DOWNLOAD_LINKS_TABLE)
    already = self._already_written_samples(self.DOWNLOAD_LINKS_TABLE)

    added = 0
    with open(self.DOWNLOAD_LINKS_TABLE, "a") as f:
        for s in samples:
            if s in mapping and s not in already:
                f.write(f"{s}\t{mapping[s]}\n")
                added += 1

    print(f"[2504_high_cov] added={added} → {self.DOWNLOAD_LINKS_TABLE}")
    return added

def extract_links_from_1000G_698_related_high_coverage(self):
    """
    Uses:
      self.SAMPLES_FILE
      self.INDEX2_FILE
      self.DOWNLOAD_LINKS_TABLE  (headered TSV with columns: SAMPLE_NAME, ENA_FILE_PATH)

    Appends new rows and reports how many were added.
    """
    # read sample names and scan index text
    samples = self._load_samples()
    wanted = set(samples)
    index_text = Path(self.INDEX2_FILE).read_text(errors="ignore")
    mapping = self._scan_index_text_for_links(index_text, wanted)

    # ensure TSV exists with header, then append only new samples
    self._ensure_tsv_header(self.DOWNLOAD_LINKS_TABLE)
    already = self._already_written_samples(self.DOWNLOAD_LINKS_TABLE)

    added = 0
    with open(self.DOWNLOAD_LINKS_TABLE, "a") as f:
        for s in samples:
            if s in mapping and s not in already:
                f.write(f"{s}\t{mapping[s]}\n")
                added += 1

    print(f"[698_related_high_cov] added={added} → {self.DOWNLOAD_LINKS_TABLE}")
    return added

def extract_cram_files_download_links_from_index_files(self):
  self.extract_links_from_1000G_2504_high_coverage()
  self.extract_links_from_1000G_698_related_high_coverage()

def _build_paths_txt(self):
    """
    Read self.DOWNLOAD_LINKS_TABLE (TSV with ENA_FILE_PATH),
    strip the ftp host to relative ENA paths, and write to self.PATHS_TXT.
    """
    df = pd.read_csv(self.DOWNLOAD_LINKS_TABLE, sep="\t", dtype=str)
    # build LOCAL_PATH by removing ftp host prefix (works for ftp://ftp.sra.ebi.ac.uk)
    df["LOCAL_PATH"] = df["ENA_FILE_PATH"].str.replace(
        r"^ftp://ftp\.sra\.ebi\.ac\.uk", "", regex=True
    )
    Path(self.PATHS_TXT).parent.mkdir(parents=True, exist_ok=True)
    df["LOCAL_PATH"].to_csv(self.PATHS_TXT, index=False, header=False)
    print(f"[paths] wrote {len(df)} lines → {self.PATHS_TXT}")

def _write_aspera_script(self):
    """
    Write an ascp batch script to self.ASPERA_SCRIPT.
    Uses Singularity and environment from self.* constants.
    """
    body = textwrap.dedent(f"""\
    #!/bin/bash
    # Usage: ./$(basename "$0") "{shlex.quote(self.OUTDIR)}"
    set -euo pipefail

    # Modules / containers
    module load {self.MODULE_SINGULARITY}

    # Container + Aspera bits
    SIF="{self.SIF}"
    ASPERA_KEY="{self.ASPERA_KEY}"

    # Export env into the container so ascp is on PATH and libs are found
    export SINGULARITYENV_PATH="{self.SINGULARITYENV_PATH}"
    export SINGULARITYENV_LD_LIBRARY_PATH="{self.SINGULARITYENV_LD_LIBRARY_PATH}"

    OUTDIR="${{1:-{shlex.quote(self.OUTDIR)}}}"
    mkdir -p "$OUTDIR"

    # Parallel batch download with xargs
    # Reads relative ENA paths from: {self.PATHS_TXT}
    function download_batch {{
      xargs -a "{self.PATHS_TXT}" -P {int(self.XARGS_PROCS)} -I{{}} \\
        singularity exec -B /hpc/packages:/hpc/packages "$SIF" \\
          ascp -i "$ASPERA_KEY" -P33001 -O33001 -k2 -QT --overwrite=diff \\
          "era-fasp@fasp.sra.ebi.ac.uk:{{}}" "$OUTDIR/"
    }}

    download_batch
    """)
    Path(self.ASPERA_SCRIPT).parent.mkdir(parents=True, exist_ok=True)
    Path(self.ASPERA_SCRIPT).write_text(body)
    os.chmod(self.ASPERA_SCRIPT, 0o755)
    print(f"[script] wrote aspera batch script → {self.ASPERA_SCRIPT}")

def launch_aspera_batch_to_download_cram_files(self):
    """
    Build paths file, write the script, then submit via bsub.
    Prints LSF submission output; returns subprocess returncode.
    """
    # Prep inputs
    self._build_paths_txt()
    self._write_aspera_script()

    # LSF command (interactive shell wrapper that backgrounds the job)
    bsub_cmd = (
        f'bsub -P {shlex.quote(self.LSF_ACCOUNT)} '
        f'-q {shlex.quote(self.LSF_QUEUE)} '
        f'-W {shlex.quote(self.LSF_WALLTIME)} '
        f'-n {int(self.LSF_CORES)} '
        f'-J {shlex.quote(self.LSF_JOB_NAME)} '
        f'-Is /bin/bash -lc '
        f'"nohup bash {shlex.quote(self.ASPERA_SCRIPT)} {shlex.quote(self.OUTDIR)} '
        f'> {shlex.quote(self.LSF_JOB_NAME)}.$LSB_JOBID.out 2>&1 & exit"'
    )

    print(f"[bsub] {bsub_cmd}")
    out = subprocess.run(bsub_cmd, shell=True, text=True, capture_output=True)
    if out.stdout.strip():
        print(out.stdout.strip())
    if out.stderr.strip():
        print(out.stderr.strip())
    return out.returncode

ShortReads.download_sample_names = download_sample_names
ShortReads.get_short_read_urls = get_short_read_urls
ShortReads.download_index_files = download_index_files
ShortReads.remove_comments_from_index_files = remove_comments_from_index_files
ShortReads._load_samples = _load_samples
ShortReads._scan_index_text_for_links = _scan_index_text_for_links
ShortReads._ensure_tsv_header = _ensure_tsv_header
ShortReads._already_written_samples = _already_written_samples
ShortReads.extract_links_from_1000G_2504_high_coverage = extract_links_from_1000G_2504_high_coverage
ShortReads.extract_links_from_1000G_698_related_high_coverage = extract_links_from_1000G_698_related_high_coverage
ShortReads.extract_cram_files_download_links_from_index_files = extract_cram_files_download_links_from_index_files
ShortReads._build_paths_txt = _build_paths_txt
ShortReads._write_aspera_script = _write_aspera_script
ShortReads.launch_aspera_batch_to_download_cram_files = launch_aspera_batch_to_download_cram_files
