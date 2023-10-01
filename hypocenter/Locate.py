import os
from glob import glob
from pathlib import Path
from shutil import copy

from core.Catalog import catalog2xyzm
from core.Magnitude import setMagnitude
from tqdm import tqdm
from obspy import read_events

from hypocenter.Station import toSTATION0HYP


def locateHypocenter(config):
    locationPath = os.path.join("results", "location", "hypocenter")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    cmd = "cat results/*_*.out > results/all.out"
    os.system(cmd)
    catalog = read_events(os.path.join("results", "all.out"))
    toSTATION0HYP(config, catalog)
    catalogs = glob(os.path.join("results", "all.out"))
    desc = "+++ Locate catalog using 'Hypocenter' ..."
    for catalogFile in tqdm(catalogs, desc=desc):
        copy(catalogFile, locationPath)
    copy(os.path.join("files", "select.inp"), locationPath)
    copy(os.path.join("results", "STATION0.HYP"), locationPath)
    root = os.getcwd()
    os.chdir(locationPath)
    with open("hyp.inp", "w") as f:
        f.write("all.out\nn\n")
    cmd = "hyp < hyp.inp >/dev/null 2>/dev/null"
    os.system(cmd)
    setMagnitude("hyp.out", config)
    catalog2xyzm(config, "hyp.out", "hypocenter")
    cmd = "select select.inp >/dev/null 2>/dev/null"
    os.system(cmd)
    os.chdir(root)
