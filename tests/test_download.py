from os import path as op

import xarray as xr

from analysis_functions.download import download_file

url = "https://sandbox.zenodo.org/record/1186200/files/e067015e_c743_201706.nc"


def test_download():
    filename = download_file(url)
    assert filename == op.basename(url)
    ds = xr.open_dataset(filename)
    print(ds)
    return ds


if __name__ == "__main__":
    test_download()
