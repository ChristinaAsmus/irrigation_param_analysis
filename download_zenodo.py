import os
import tarfile
from os import path as op

import xarray as xr

from analysis_functions.download import download_file_in_dir

files = ["067017.tar.gz", "067019.tar.gz", "067020.tar.gz"]  # , "067015", "067016"]
url = "https://zenodo.org/record/7867329/files/"

# create folder for data
dir_working = os.getcwd()
# creates dir in parent directory
# dir_out = os.path.join(os.pardir, "data")
dir_out = "data"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
print("Data directory is: ", dir_out)


def download(singlefile):
    print("Downloading", str(singlefile))
    filename = download_file_in_dir(url + str(singlefile), dir_out)
    assert filename == op.basename(url + str(singlefile))
    return filename


for singlefile in files:
    download(singlefile)
    # unpack tarball
    file = tarfile.open(os.path.join(dir_out, singlefile))
    file.extractall(str(dir_out))
