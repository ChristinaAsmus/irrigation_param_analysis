import xarray as xr

from analysis_functions.download import download_file

url = "https://sandbox.zenodo.org/record/1185969/files/e065000e_c167_198001.nc"


filename = download_file(url)
ds = xr.open_dataset(filename)
print(ds)
