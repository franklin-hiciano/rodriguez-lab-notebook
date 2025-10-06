class ShortReadsDownloadLinks:
    """
    Desc: Creates a .tsv file relating sample names and download links for short-read data from 1000 Genomes project.
    Platforms: Minerva HPC at Mount Sinai
    """

    def __init__(self):
      self.sample_names_txt_file = "/sc/arion/work/hiciaf01/imputation/sample_names.txt"
      self.outdir = "/sc/arion/work/hiciaf01/imputation/short_read_download_links/"

      self.index_1_file_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
      self.index_2_file_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index"

      self.index_1_file_name = "1000G_2504_high_coverage.sequence.index"
      self.index_2_file_name = "1000G_698_related_high_coverage.sequence.index"
      self.index_1_file = os.path.join(self.outdir, self.index_1_file_name)
      self.index_2_file = os.path.join(self.outdir, self.index_2_file_name)

      self.download_links_tsv_file_name = "cram_download_links.tsv"
      self.download_links_tsv_file = os.path.join(self.outdir, self.download_links_tsv_file_name)

    def run(self):
        """
        Download the index files, extract the download links from them, and 
        a .tsv file relating each sample name with its short-reads download link.
        """
        print("Running ShortReadsDownloadLinks...")
        self.download_index_files_containing_download_links()
        self.write_download_links_tsv_file()
        print("ShortReadsDownloadLinks finished.")

def download_index_files_containing_download_links(self):
  """
  Downloads the two index files containing download links for the short read samples.
  """
  print(f"Downloading the two index files. Destinations: {self.index_1_file}, {self.index_2_file}. Sources: {self.index_1_file_url}, {self.index_2_file_url}.")
  urlretrieve(self.index_1_file_url, self.index_1_file)
  urlretrieve(self.index_2_file_url, self.index_2_file)
  print(f"Downloaded the two index files. Destinations: {self.index_1_file}, {self.index_2_file}. Sources: {self.index_1_file_url}, {self.index_2_file_url}.")

def write_download_links_tsv_file(self):
  """
  Create a TSV with two columns: sample <tab> url
  - Reads samples from self.sample_names_txt_file (one per line; blanks ignored).
  - Searches BOTH index files for the FIRST matching link for each sample.
  - If no link is found, leaves the url cell blank.
  """

  # -------------------------------------------------
  # 1) GET THE SAMPLE LIST
  # -------------------------------------------------
  print("[1] loading sample list…")
  with open(self.sample_names_txt_file) as f:
    samples = [s.strip() for s in f if s.strip()]
  print(f"[1] {len(samples)} samples loaded")

  # -------------------------------------------------
  # 2) HELPER: get_single_download_link(sample)
  #    Searches BOTH index files in one go,
  #    clearly indicating which file is being searched.
  # -------------------------------------------------
  def get_single_download_link(sample):
    def find_in_index(index_path):
      print(f"    [search] {sample} in {index_path}")
      with open(index_path) as fh:
        # find header
        sample_col = path_col = None
        for line in fh:
          if line.startswith("#ENA_FILE_PATH"):
            header = line[1:].rstrip("\n").split("\t")
            sample_col = header.index("SAMPLE_NAME")
            path_col = header.index("ENA_FILE_PATH")
            break
        # scan rows
        for line in fh:
          if not line or line.startswith("#"):
            continue
          parts = line.split("\t")
          if len(parts) > sample_col and parts[sample_col] == sample:
            return parts[path_col].strip()
      return None

    # search order: index_1_file → index_2_file
    for path in [self.index_1_file, self.index_2_file]:
      url = find_in_index(path)
      if url:
        print(f"    [hit] {sample} → {url}")
        return url
    print(f"    [miss] {sample} not found")
    return ""

  # -------------------------------------------------
  # 3) WRITE THE RESULTS TO TSV
  # -------------------------------------------------
  print(f"[3] writing results → {self.download_links_tsv_file}")
  with open(self.download_links_tsv_file, "w") as out:
    out.write("sample\turl\n")
    for s in samples:
      out.write(f"{s}\t{get_single_download_link(s)}\n")
  print("[done] TSV written")

ShortReadsDownloadLinks.download_index_files_containing_download_links = download_index_files_containing_download_links
ShortReadsDownloadLinks.write_download_links_tsv_file = write_download_links_tsv_file
