import os
import subprocess
import sys
import time
from Bio import Entrez

# Set email for NCBI API
Entrez.email = "your_email@example.com"


def display_progress(step_name, duration=2, sleep=0.1, _break_=False):
    """
    Displays a progress arrow with a status message.

    Args:
        step_name (str): Description of the current step.
        duration (int): Duration to simulate progress (seconds).
    """
    progress_bar = "["
    for _ in range(duration * 5):  # Updates approximately every 0.2 seconds
        progress_bar += "-"
        print(f"\r{progress_bar}>] {step_name}", end="", flush=True)
        time.sleep(sleep)
    print(f" [done].")

    if _break_:
        time.sleep(0.5)


def fetch_geo_samples(geo_id):
    """
    Fetches sample information for a given GEO ID.

    Args:
        geo_id (str): GEO ID (e.g., "GSE123456").

    Returns:
        list: A list of samples with Accession and Title.
    """
    try:
        display_progress(f"Searching for internal UID for GEO ID: {geo_id}", _break_=True)
        search_handle = Entrez.esearch(db="gds", term=geo_id)
        search_result = Entrez.read(search_handle)
        search_handle.close()

        if not search_result["IdList"]:
            display_progress("No results found for the given GEO ID.")
            return None

        uid = search_result["IdList"][0]
        display_progress(f"Found internal UID: {uid}", _break_=True)

        display_progress(f"Fetching detailed information for UID: {uid}", _break_=True)
        summary_handle = Entrez.esummary(db="gds", id=uid)
        summary_result = Entrez.read(summary_handle)
        summary_handle.close()

        samples = summary_result[0].get("Samples", [])
        if not samples:
            print("No samples found for the given GEO ID.")
            return None

        return samples
    except Exception as e:
        print(f"Error fetching GEO samples: {e}")
        return None


def download_with_prefetch(accession, output_dir, max_retries=3):
    """
    Downloads SRA files using prefetch.

    Args:
        accession (str): The sample accession number.
        output_dir (str): Directory to save the SRA files.
        max_retries (int): Maximum number of retries for the download.

    Returns:
        bool: True if the download succeeds, False otherwise.
    """
    attempt = 0
    while attempt < max_retries:
        try:
            print(f"Attempt {attempt + 1}/{max_retries} for downloading {accession} using prefetch...")
            prefetch_command = ["prefetch", accession, "--output-directory", output_dir]
            subprocess.run(prefetch_command, check=True)
            print(f"Successfully downloaded {accession} to {output_dir}.")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Download failed for {accession}: {e}")
            attempt += 1
            time.sleep(5)
    print(f"Exceeded maximum retries for {accession}. Skipping.")
    return False


def convert_all_sra_to_fastq(output_dir):
    """
    Converts all .sra files in the output directory to FASTQ using fasterq-dump.

    Args:
        output_dir (str): Directory containing the SRA files.

    Returns:
        None
    """
    print(f"Converting all SRA files in {output_dir} to FASTQ...")
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".sra"):
                sra_file_path = os.path.join(root, file)
                print(f"Found SRA file: {sra_file_path}. Starting conversion...")
                try:
                    fasterq_command = [
                        "fasterq-dump",
                        "--split-files",
                        "--gzip",
                        "-O",
                        output_dir,
                        sra_file_path,
                    ]
                    subprocess.run(fasterq_command, check=True)
                    print(f"Successfully converted {sra_file_path} to FASTQ.")
                except subprocess.CalledProcessError as e:
                    print(f"Failed to convert {sra_file_path}: {e}")


def validate_args(args):
    """
    Validates command-line arguments.

    Args:
        args (list): Command-line arguments.

    Returns:
        tuple: GEO ID and output directory if valid, otherwise exits the program.
    """
    if len(args) != 3:
        print("\nInvalid Arguments")
        print("\nUsage: python script.py <GEO_ID> <OUTPUT_DIR>")
        sys.exit(1)

    geo_id = args[1]
    output_dir = args[2]

    if not os.path.isdir(output_dir):
        print(f"Error: Output directory '{output_dir}' does not exist or is not valid.")
        sys.exit(1)

    return geo_id, output_dir


def main():
    geo_id, output_dir = validate_args(sys.argv)
    samples = fetch_geo_samples(geo_id)

    if samples:
        display_progress(f"Samples found for {geo_id}:", _break_=True)
        for sample in samples:
            display_progress(f"Accession: {sample['Accession']}, Title: {sample['Title']}", sleep=0.005)
        print(f"\nFound {len(samples)} total samples")

        # Download all SRA files first
        for sample in samples:
            accession = sample["Accession"]
            print(f"\nDownloading sample: {accession}")
            download_with_prefetch(accession, output_dir)

        # Convert all downloaded SRA files to FASTQ
        convert_all_sra_to_fastq(output_dir)
    else:
        print("No samples found for the given GEO ID.")


if __name__ == "__main__":
    main()
