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

def submit_job(self, ftp_url=None, udp_url=None, outdir="/sc/", 
               job_name="aspera_test", allocation_account="acc_oscarlr",
               queue="express", run_live=False):
    import os, subprocess, textwrap
    os.makedirs(outdir, exist_ok=True)

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

ASPERA_PROGRAM=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/bin/ascp
ASPERA_KEY=/hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh

mkdir -p {outdir}

echo "Testing wget fallback..."
wget -O {outdir}/$(basename {ftp_url}) {ftp_url}

echo "✅ Done"
""")


    print("=== SCRIPT PREVIEW ===")
    print(script_content)
    print("======================")

    script_path = f"./{job_name}.sh"
    with open(script_path, "w", newline="\n") as f:
      f.write(script_content)

    os.chmod(script_path, 0o755)

    if run_live:
        # Just execute script directly
        subprocess.run([script_path], check=True)
    else:
        # Submit with bsub
        result = subprocess.run(f"bsub < {script_path}",
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
