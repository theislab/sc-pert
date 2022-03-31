import requests
import os

def download_binary_file(
    file_url: str, output_path: str, overwrite: bool = False
) -> None:
    """
    Download binary data file from a URL.

    Args:
    ----
        file_url: URL where the file is hosted.
        output_path: Output path for the downloaded file.
        overwrite: Whether to overwrite existing downloaded file.

    Returns
    -------
        None.
    """
    file_exists = os.path.exists(output_path)
    if (not file_exists) or (file_exists and overwrite):
        request = requests.get(file_url)
        with open(output_path, "wb") as f:
            f.write(request.content)
        print(f"Downloaded data from {file_url} at {output_path}")
    else:
        print(
            f"File {output_path} already exists. "
            "No files downloaded to overwrite the existing file."
        )