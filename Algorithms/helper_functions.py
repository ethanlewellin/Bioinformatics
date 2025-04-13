### Imports
import os
import requests
import tarfile
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
        else:
            print(f"Failed to download {filename}, status code: {response.status_code}")

    # List final contents
    print("Final folder contents:", os.listdir(folder_path))
