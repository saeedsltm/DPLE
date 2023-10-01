#!/usr/bin/env python3

from obspy import read, read_inventory
from obspy import UTCDateTime as utc
import os
from pathlib import Path
from tqdm import tqdm
from pandas import DataFrame, Series
from numpy import array, nan
from glob import glob
from core.Extra import handle_masked_arr
from obspy.core.inventory.inventory import Inventory
from gamma.utils import estimate_station_spacing


def prepareWaveforms(starttime, endtime, config):
    path = Path("tmp")
    path.mkdir(parents=True, exist_ok=True)

    if config["preprocess_data"]:
        for f in glob(os.path.join(path, "*")):
            os.remove(f)

    # Read Station Data
    invFile = os.path.join(
        "DB",
        f"{starttime.strftime('%Y%m%d')}_{endtime.strftime('%Y%m%d')}",
        "stations",
        "*.xml")

    try:
        inventory = read_inventory(invFile)
        inv = Inventory()
        for network in config["networks"]:
            inv += inventory.select(network)
    except Exception:
        return

    inv.write(os.path.join("tmp", "stations.xml"), format="STATIONXML")

    # Get Station Codes
    stations = sorted(set([
        s.split(".")[1] for s in inv.get_contents()["channels"]
    ]))

    with open(os.path.join("tmp", "mseed.csv"), "w") as fp:
        fp.write("fname,E,N,Z\n")
        for station in tqdm(stations, desc="+++ Preparing raw data"):
            st = read(os.path.join(
                "DB",
                f"{starttime.strftime('%Y%m%d')}_{endtime.strftime('%Y%m%d')}",
                "waveforms",
                f"??.{station}.*.???__{starttime.strftime('%Y%m%d')}T000000Z__{endtime.strftime('%Y%m%d')}T000000Z.mseed"),
                format="MSEED", check_compression=False)
            st.merge(fill_value=None)
            st = handle_masked_arr(st)
            if len(config["searchWindow"]):
                s, t = config["searchWindow"]
                st = st.slice(utc(s), utc(t))
            if len(st):
                st = st.slice(utc(starttime), utc(endtime))
                st.write(os.path.join("tmp", f"{station}.mseed"))
                sta = st[0].stats.station
                chn = st[0].stats.channel[:-1]
                fp.write(
                    f"{sta}.mseed,{chn}E,{chn}N,{chn}Z\n")
    return True


def prepareInventory(config, proj, st, et):
    stationxml = os.path.join(
        "DB",
        f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}",
        "stations",
        "*.xml")
    inv = read_inventory(stationxml)
    station_df = []
    for net in inv:
        for sta in net:
            station_df.append({
                "id": f"{net.code}.{sta.code}.",
                "longitude": sta.longitude,
                "latitude": sta.latitude,
                "elevation(m)": sta.elevation,
                "unit": "m/s",
                "component": "E,N,Z",
            })
            break
    station_df = DataFrame(station_df)
    cx1 = station_df["longitude"] >= config["xlim_degree"][0]
    cx2 = station_df["longitude"] < config["xlim_degree"][1]
    cy1 = station_df["latitude"] >= config["ylim_degree"][0]
    cy2 = station_df["latitude"] < config["ylim_degree"][1]
    c = (cx1) & (cx2) & (cy1) & (cy2)
    station_df = station_df[c]
    station_df.reset_index(inplace=True, drop=True)
    station_df[["x(km)", "y(km)"]] = station_df.apply(
        lambda x: Series(
            proj(longitude=x.longitude, latitude=x.latitude)),
        axis=1)
    station_df["z(km)"] = station_df["elevation(m)"].apply(lambda x: -x*1e-3)
    station_dict = {station: (x, y) for station, x, y in zip(
        station_df["id"], station_df["x(km)"], station_df["y(km)"])}
    return station_df, station_dict


def picks2DF(picks, pick_outfile):
    outFile = os.path.join("results", pick_outfile)
    pick_df = []
    for pick in picks:
        pick_df.append(
            {"id": pick.trace_id,
             "timestamp": pick.peak_time.datetime,
             "prob": pick.peak_value,
             "type": pick.phase.lower(),
             "amp": nan,  # pick.phase_amplitude,
             "phase_amp": nan  # pick.phase_amplitude
             })
    pick_df = DataFrame(pick_df)
    pick_df.to_csv(outFile, index=False, float_format="%0.3f")


def applyGaMMaConfig(config, stations):
    method = {
        "BGMM": 8,
        "GMM": 1}
    config["oversample_factor"] = method[config["method"]]

    config["vel"] = {"p": 6.0, "s": 6.0 / 1.75}
    config["dims"] = ["x(km)", "y(km)", "z(km)"]
    clon = config["center"][0]
    clat = config["center"][1]
    config["x(km)"] = (array(config["xlim_degree"]) -
                       array(clon))*config["degree2km"]
    config["y(km)"] = (array(config["ylim_degree"]) -
                       array(clat))*config["degree2km"]
    config["z(km)"] = (config["zlim_degree"][0], config["zlim_degree"][1])
    config["bfgs_bounds"] = (
        (config["x(km)"][0] - 1, config["x(km)"][1] + 1),
        (config["y(km)"][0] - 1, config["y(km)"][1] + 1),
        (config["z(km)"][0], config["z(km)"][1] + 1),
        (None, None),
    )

    if config["useEikonal"]:
        zz = config["zz"]
        vp = config["vp"]
        vp_vs_ratio = config["vp_vs_ratio"]
        vs = [v / vp_vs_ratio for v in vp]
        h = config["h"]
        vel = {"z": zz, "p": vp, "s": vs}
        config["eikonal"] = {"vel": vel,
                             "h": h,
                             "xlim": config["x(km)"],
                             "ylim": config["y(km)"],
                             "zlim": config["z(km)"]}
    if config["dbscan_eps"] == "A":
        config["dbscan_eps"] = 2.5 * estimate_station_spacing(stations) / config["vel"]["p"]


    return config
