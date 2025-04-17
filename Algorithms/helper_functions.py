### Imports
import os
import requests
import tarfile
import zipfile
import gzip
import shutil
import pandas as pd
from typing import Dict

### Get 10x data from URL
def download_and_extract_10x_data(data_path: str, folder_path: str, urls: Dict[str, str]) -> None:
    """
    Downloads files from the given URLs and extracts any .tar.gz archives.

    Args:
        data_path (str): Base directory path where data will be stored.
        folder_path (str): Specific subfolder path where files will be downloaded and extracted.
        urls (Dict[str, str]): A dictionary mapping filenames to their download URLs.

    Example:
        urls = {
            "data.h5": "https://example.com/data.h5",
            "archive1.tar.gz": "https://example.com/archive1.tar.gz",
            "archive2.tar.gz": "https://example.com/archive2.tar.gz"
        }
        download_and_extract_10x_data("./data/", "./data/my_folder/", urls)

    Notes:
        - The function creates `folder_path` if it doesn't exist.
        - Any files ending with `.tar.gz` will be automatically extracted after download.
        - HTTP status codes are checked to ensure successful downloads.
    """
    # Create the folder path if it doesn't exist
    try:
        os.makedirs(folder_path, exist_ok=True)
        print(f"Folder path '{folder_path}' created successfully.")
    except Exception as e:
        print(f"An error occurred while creating the folder path: {e}")

    # Download files from the URLs
    for filename, url in urls.items():
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"{filename} successfully downloaded")

            # If the file is a .tar.gz, extract it
            if filename.endswith(".tar.gz"):
                try:
                    with tarfile.open(file_path, "r:gz") as tar:
                        tar.extractall(path=folder_path)
                    print(f"Extracted {filename} to: {folder_path}")
                except Exception as e:
                    print(f"An error occurred while extracting {filename}: {e}")
            
            # If the file is a .zip, extract it
            elif filename.endswith(".zip"):
                try:
                    with zipfile.ZipFile(file_path, 'r') as zip_ref:
                        zip_ref.extractall(folder_path)
                    print(f"Extracted {filename} to: {folder_path}")
                except Exception as e:
                    print(f"An error occurred while extracting {filename}: {e}")

        else:
            print(f"Failed to download {filename}, status code: {response.status_code}")

    # List final contents
    print("Final folder contents:", os.listdir(folder_path))


### Decompress .csv.gz File
def decompress_file(gz_path: str) -> None:
    """
    Decompresses a .gz file to its original CSV form.

    Args:
        gz_path (str): Path to the .gz compressed file.
    """
    csv_path = gz_path.rstrip('.gz')
    with gzip.open(gz_path, 'rb') as f_in:
        with open(csv_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed {gz_path} to {csv_path}")
    
### Load and Inspect 10x Data (Generic)
def load_and_inspect_10x_csv_data(sample_path: str) -> Dict[str, pd.DataFrame]:
    """
    Loads and inspects all CSV/CSV.GZ files in a given 10x data sample directory.
    Automatically decompresses any .csv.gz files if needed.

    Args:
        sample_path (str): Path to the unzipped 10x data sample folder.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary where keys are file names (without extension)
                                 and values are corresponding pandas DataFrames.
    """
    data_frames = {}

    for file_name in os.listdir(sample_path):
        file_path = os.path.join(sample_path, file_name)

        if file_name.endswith(".csv.gz"):
            csv_path = file_path.rstrip(".gz")

            # Decompress if .csv file doesn't exist
            if not os.path.exists(csv_path):
                decompress_file(file_path)

            # Load DataFrame
            df = pd.read_csv(csv_path)
            key_name = os.path.splitext(os.path.basename(csv_path))[0]
            data_frames[key_name] = df
            print(f"Loaded {key_name} (from .csv.gz) with shape {df.shape}")

        elif file_name.endswith(".csv"):
            df = pd.read_csv(file_path)
            key_name = os.path.splitext(os.path.basename(file_path))[0]
            data_frames[key_name] = df
            print(f"Loaded {key_name} (from .csv) with shape {df.shape}")

    # Quick inspection — show heads
    for name, df in data_frames.items():
        print(f"\n{name} preview:")
        print(df.head())

    return data_frames