class Samples66:
    """
    Desc: This class downloads the 66 samples that Nefte found matched between short-read and long-read data.
    Platforms: Minerva HPC at Mount Sinai School of Medicine
    """

    def __init__(self):
        self.output_file_name = "sample_names.txt"
        self.output_dir = "/sc/arion/work/hiciaf01/imputation/"
        self.output_file = os.path.join(self.output_dir, self.output_file_name)

        self.google_drive_file_id = "1ArzWBQqjSpl_KNQFgPW-at0l0BTZnvvu"
        self.google_drive_file_url = f"https://drive.google.com/uc?id={self.google_drive_file_id}"

    def run(self):
        """
        Downloads the sample names.
        """
        print("Running Samples66...")
        self.download()
        print("Successfully ran Samples66.")

def download(self):
    """
    This function downloads the 66 samples that Nefte found matched between short-read and long-read data.
    """
    print(f"Downloading 66 sample names. Destination: {self.output_file}. Source: {self.google_drive_file_url}.")
    gdown.download(self.google_drive_file_url, output=self.output_file, quiet=False)
    print(f"66 sample names successfully downloaded. Destination: {self.output_file}. Source: {self.google_drive_file_url}.")

Samples66.download = download
