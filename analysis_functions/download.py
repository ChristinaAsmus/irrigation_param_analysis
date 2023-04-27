from os import path as op

import requests
from tqdm import tqdm


def download_file(url):
    local_filename = op.basename(url)
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size_in_bytes = int(r.headers.get("content-length", 0))
        chunk_size = 8192
        progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
        with open(local_filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                progress_bar.update(len(chunk))
                # if chunk:
                f.write(chunk)
    return local_filename


def download_file_in_dir(url, data_dir):
    local_filename = op.basename(url)
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size_in_bytes = int(r.headers.get("content-length", 0))
        chunk_size = 8192
        progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
        with open(str(data_dir) + "/" + str(local_filename), "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                progress_bar.update(len(chunk))
                # if chunk:
                f.write(chunk)
    return local_filename
