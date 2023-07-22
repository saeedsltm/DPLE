import os
from glob import glob
from pathlib import Path
from shutil import copy

from core.Catalog import catalog2xyzm
from tqdm import tqdm

from hypodd.Extra import writexyzm
from hypodd.Input import prepareHypoddInputs


def locateHypoDD(config):
    locationPath = os.path.join("results", "location", "hypoDD")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    catalogs = glob(os.path.join("results", "location", "hypocenter", "hyp.out"))
    for catalogFile in catalogs:
        copy(catalogFile, locationPath)
    root = os.getcwd()
    os.chdir(locationPath)
    desc = "+++ Locate catalog using 'HypoDD' ..."
    for catalogFile in tqdm(glob("hyp.out"), desc=desc):
        outName = catalogFile.split(os.sep)[-1].split(".")[0]
        outName = "hypoDD"
        prepareHypoddInputs(config,
                            catalogFile,
                            locationPath)
        cmd = "ph2dt ph2dt.inp >/dev/null 2>/dev/null"
        os.system(cmd)
        cmd = "hypoDD hypoDD.inp >/dev/null 2>/dev/null"
        os.system(cmd)
        writexyzm(outName)
        for f in glob("hypoDD.reloc*"):
            os.remove(f)
    os.chdir(root)
