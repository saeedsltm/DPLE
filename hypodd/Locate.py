import os
from glob import glob
from pathlib import Path
from shutil import copy

from tqdm import tqdm

from hypodd.Extra import writexyzm
from hypodd.Input import prepareHypoddInputs
from hypodd.Extra import readHypoddConfig
from core.Extra import divide_chunks
from obspy import read_events
from pandas import concat, read_csv


def readxyzm(inp):
    return read_csv(inp, delim_whitespace=True)


def locateHypoDD(config):
    locationPath = os.path.join("results", "location", "hypoDD")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    catalogs = glob(os.path.join("results", "location", "hypocenter", "select.out"))
    hypoddConfig = readHypoddConfig()
    for catalogFile in catalogs:
        copy(catalogFile, locationPath)
    root = os.getcwd()
    os.chdir(locationPath)
    print("+++ Loading catalog ...")
    catalog = read_events("select.out")
    desc = "+++ Locating catalog using 'HypoDD' ... "
    for c, chunck_catalog in tqdm(enumerate(divide_chunks(catalog, 1e4)), desc=desc):
        print(f"+++ Working on chunck {c+1} ...")
        outName = f"hypoDD_{c}"
        prepareHypoddInputs(config,
                            hypoddConfig,
                            chunck_catalog,
                            locationPath)
        cmd = "ph2dt ph2dt.inp >/dev/null 2>/dev/null"
        os.system(cmd)
        cmd = "hypoDD hypoDD.inp >/dev/null 2>/dev/null"
        os.system(cmd)
        writexyzm(outName)
        for f in glob("hypoDD.reloc*"):
            os.remove(f)
    xyzmfiles = glob("xyzm_hypoDD_*.dat")
    df = concat(map(readxyzm, xyzmfiles))
    columns = ["ORT", "Lon", "Lat", "Dep", "Mag",
               "Nus", "NuP", "NuS", "ADS", "MDS", "GAP", "RMS", "ERH", "ERZ"]
    with open("xyzm_hypoDD.dat", "w") as f:
        df.to_string(f, columns=columns, index=False, float_format="%7.3f")
    for f in xyzmfiles:
        os.remove(f)
    os.chdir(root)
