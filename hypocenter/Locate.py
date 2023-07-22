import os
from glob import glob
from pathlib import Path
from shutil import copy

from core.Catalog import catalog2xyzm
from tqdm import tqdm

from hypocenter.Station import toSTATION0HYP


def locateHypocenter(config):
    toSTATION0HYP(config)
    locationPath = os.path.join("results", "location", "hypocenter")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    cmd = "cat results/*_*.out > results/all.out"
    os.system(cmd)
    catalogs = glob(os.path.join("results", "all.out"))
    desc = "+++ Locate catalog using 'Hypocenter' ..."
    for catalogFile in tqdm(catalogs, desc=desc):
        copy(catalogFile, locationPath)
    copy(os.path.join("results", "STATION0.HYP"), locationPath)
    root = os.getcwd()
    os.chdir(locationPath)
    with open("hyp.inp", "w") as f:
        f.write("all.out\nn\n")
    cmd = "hyp < hyp.inp >/dev/null 2>/dev/null"
    os.system(cmd)
    catalog2xyzm(config, "hyp.out", "hypocenter")
    os.chdir(root)
