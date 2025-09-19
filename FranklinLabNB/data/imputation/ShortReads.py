import subprocess
import textwrap
import os
import sys
import gdown
from urllib.request import urlretrieve

class ShortReads:
  def run(self):
    # self.download_sample_names()
    # self.test_short_read_downloads_without_live_logging()
    self.submit_job(ftp_url="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram", outdir="/sc/arion/scratch/hiciaf01", run_live=True)

def download_sample_names(self):
  file_id = "1ArzWBQqjSpl_KNQFgPW-at0l0BTZnvvu"
  url = f"https://drive.google.com/uc?id={file_id}"
  output = "sample_names.txt"  # name to save as
  gdown.download(url, output, quiet=False)

def get_short_read_urls(self):
  outpath = "urls.tsv"
  urlretrieve("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index")
  df = pd.read_csv('urls.tsv', sep='\t', comment='#')

def test_short_read_downloads_without_live_logging(self):
    os.makedirs("downloads", exist_ok=True)

    bash_script = textwrap.dedent("""\
    #!/bin/bash
    set -euo pipefail

    ASPERA_PROGRAM=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/bin/ascp
    ASPERA_KEY=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh
    ftp_url=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram
    udp_url=vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz
    echo ${ASPERA_PROGRAM}

    echo ${ASPERA_KEY}

    function download {
      wget ${ftp_url} ./testing
      # ${ASPERA_PROGRAM} \
        # -i ${ASPERA_KEY} \
        # -P33001 -O33001 -T -L- \
        # era-fasp@fasp.sra.ebi.ac.uk:${udp_url} \
        # /sc/arion/scratch/hiciaf01/test.cram/
    }

    download
    """)



    result = subprocess.run(
        bash_script,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("❌ Script failed with exit code", result.returncode)
        print("\n--- STDOUT ---\n", result.stdout)
        print("\n--- STDERR ---\n", result.stderr)
    else:
        print("✅ Script succeeded!")
        print("\n--- STDOUT ---\n", result.stdout)

def submit_job(self,
               ftp_url=None,                 # e.g. ftp://ftp.sra.ebi.ac.uk/vol1/run/...
               udp_url=None,                 # e.g. vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz
               outdir="/sc/",
               job_name="aspera_test",
               allocation_account="acc_oscarlr",
               queue="express",
               run_live=False,
               use_ascp_first=False,
               proxy_url=None,               # e.g. http://proxy.chimera.mssm.edu:3128 (example)
               extra_wget_args=None):
    """
    Submit an LSF job that (optionally) tries Aspera first, then falls back to HTTPS wget.
    - Automatically converts EBI ftp:// to https:// (safer on locked-down queues).
    - Supports HTTP(S) proxy via env vars and wget --execute settings.
    - Avoids non-ASCII emojis in the script to prevent mojibake.
    """
    import os, subprocess, textwrap, shlex

    os.makedirs(outdir, exist_ok=True)

    # Convert ftp:// to https:// when it points to EBI; keeps filename identical.
    https_url = None
    if ftp_url:
        if ftp_url.startswith("ftp://ftp.sra.ebi.ac.uk/"):
            https_url = "https://" + ftp_url[len("ftp://"):]
        else:
            # generic fallback: try https mirror if someone passed an ftp URL
            https_url = ftp_url.replace("ftp://", "https://", 1)
    # Prefer https_url if we created one
    dl_url = https_url or ftp_url

    # Assemble proxy env exports if provided
    proxy_exports = ""
    wget_proxy_exec = ""
    if proxy_url:
        # export both lowercase and uppercase envs for maximal compatibility
        proxy_exports = textwrap.dedent(f"""\
        export http_proxy={shlex.quote(proxy_url)}
        export https_proxy={shlex.quote(proxy_url)}
        export HTTP_PROXY={shlex.quote(proxy_url)}
        export HTTPS_PROXY={shlex.quote(proxy_url)}
        """)
        # also tell wget explicitly
        wget_proxy_exec = f"--execute http_proxy={shlex.quote(proxy_url)} --execute https_proxy={shlex.quote(proxy_url)}"

    # Extra wget args (retries, continue, timeouts)
    if extra_wget_args is None:
        extra_wget_args = ["--tries=3", "--timeout=60", "--continue", "--no-verbose"]
    wget_args_str = " ".join(extra_wget_args)

    # Aspera command (optional first try)
    # Note: EBI Aspera host: fasp.sra.ebi.ac.uk (or hl-xfer-fasp.ebi.ac.uk); path is usually same as ftp after host.
    ascp_lines = ""
    if use_ascp_first and udp_url:
        ascp_lines = textwrap.dedent(f"""\
        echo "[INFO] Trying Aspera (ascp) first..."
        ASPERA_PROGRAM=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/bin/ascp
        ASPERA_KEY=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh

        if [[ -x "$ASPERA_PROGRAM" ]]; then
          set +e
          "$ASPERA_PROGRAM" -QT -l 100m -i "$ASPERA_KEY" \
            fasp.sra.ebi.ac.uk:/{shlex.quote(udp_url)} {shlex.quote(outdir)}
          rc=$?
          set -e
          if [[ $rc -eq 0 ]]; then
            echo "[INFO] Aspera download succeeded."
          else
            echo "[WARN] Aspera failed with code $rc; falling back to HTTPS."
          fi
        else
          echo "[WARN] Aspera binary not found; skipping Aspera and using HTTPS."
        fi

        """)

    # Wget fallback (or primary if no Aspera)
    wget_lines = ""
    if dl_url:
        wget_lines = textwrap.dedent(f"""\
        echo "[INFO] Downloading via HTTPS with wget..."
        filename=$(basename {shlex.quote(dl_url)})
        wget {wget_args_str} {wget_proxy_exec} -O {shlex.quote(outdir)}/"$filename" {shlex.quote(dl_url)}
        echo "[INFO] Saved to: {shlex.quote(outdir)}/$filename"
        """)

    script_content = textwrap.dedent(f"""\
    #!/bin/bash
    #BSUB -P {allocation_account}
    #BSUB -J {job_name}
    #BSUB -o {job_name}.out
    #BSUB -e {job_name}.err
    #BSUB -W 0:30
    #BSUB -M 2000
    #BSUB -n 1
    #BSUB -R "rusage[mem=2000]"
    #BSUB -q {queue}

    set -euo pipefail

    # Ensure predictable locale; avoid UTF-8 emojis that can mojibake
    export LC_ALL=C.UTF-8 || true
    export LANG=C.UTF-8 || true

    mkdir -p {shlex.quote(outdir)}

    {proxy_exports.strip()}

    {ascp_lines}{wget_lines}

    echo "[INFO] Done."
    """)

    print("=== SCRIPT PREVIEW ===")
    print(script_content)
    print("======================")

    script_path = f"./{job_name}.sh"
    with open(script_path, "w", newline="\n") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)

    if run_live:
        subprocess.run(["bash", script_path], check=True)
    else:
        result = subprocess.run(f"bsub < {shlex.quote(script_path)}",
                                capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            print("❌ bsub submission failed")
            print(result.stderr)
        else:
            print("✅ Job submitted:", result.stdout.strip())
            return result.stdout.strip()

ShortReads.download_sample_names = download_sample_names
ShortReads.get_short_read_urls = get_short_read_urls
ShortReads.test_short_read_downloads_without_live_logging = test_short_read_downloads_without_live_logging
ShortReads.submit_job = submit_job
