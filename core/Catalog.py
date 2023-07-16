from core.Extra import handleNone, getRMS, getHer, getZer
from numpy import nan, round_, mean
from pandas import DataFrame
from obspy import read_events
from tqdm import tqdm


import warnings
warnings.filterwarnings("ignore")


def catalog2xyzm(config, hypInp, catalogFileName):
    """Convert catalog to xyzm file format

    Args:
        hypInp (str): file name of NORDIC file
        catalogFileName (str): file name of xyzm.dat file
    """
    print("+++ Loading catalog ...")
    cat = read_events(hypInp)
    outputFile = "xyzm_{catalogFileName:s}.dat".format(
        catalogFileName=catalogFileName.split(".")[0])
    catDict = {}
    desc = "Convert to xyzm ..."
    for i, event in tqdm(enumerate(cat), desc=desc):
        preferred_origin = event.preferred_origin()
        preferred_magnitude = event.preferred_magnitude()
        arrivals = preferred_origin.arrivals
        ort = preferred_origin.time
        lat = preferred_origin.latitude
        lon = preferred_origin.longitude
        mag = preferred_magnitude.mag if preferred_magnitude else nan
        if len(arrivals) == 0:
            continue
        try:
            dep = preferred_origin.depth*0.001
        except TypeError:
            dep = nan
        try:
            nus = handleNone(
                preferred_origin.quality.used_station_count, dtype="int")
        except AttributeError:
            nus = nan
        nuP = len(
            [arrival.phase for arrival in arrivals if "P" in arrival.phase.upper()])
        nuS = len(
            [arrival.phase for arrival in arrivals if "S" in arrival.phase.upper()])
        mds = handleNone(
            min([handleNone(arrival.distance) for arrival in preferred_origin.arrivals]), degree=True)
        ads = round_(handleNone(
            mean([handleNone(arrival.distance) for arrival in preferred_origin.arrivals]), degree=True), 2)
        try:
            gap = handleNone(
                preferred_origin.quality.azimuthal_gap, dtype="int")
        except AttributeError:
            gap = nan
        rms = getRMS(preferred_origin.arrivals)
        erh = getHer(event)
        erz = getZer(event)
        catDict[i] = {
            "ORT": ort,
            "Lon": lon,
            "Lat": lat,
            "Dep": dep,
            "Mag": mag,
            "Nus": nus,
            "NuP": nuP,
            "NuS": nuS,
            "ADS": ads,
            "MDS": mds,
            "GAP": gap,
            "RMS": rms,
            "ERH": erh,
            "ERZ": erz,
        }
    df = DataFrame(catDict).T
    df = df.replace({"None": nan})
    with open(outputFile, "w") as f:
        df.to_string(f, index=False, float_format="%7.3f")