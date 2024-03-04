#!/usr/bin/env python3

from obspy import read, read_inventory
from obspy import UTCDateTime as utc
import os
from pathlib import Path
from tqdm import tqdm
from pandas import DataFrame, Series
from glob import glob
from core.Extra import handle_masked_arr
from obspy.core.inventory.inventory import Inventory
import joblib


def saveModel(model, pick_outname):
    joblib.dump(model, os.path.join("results", f"{pick_outname}.jlib"))


def loadModel(pick_filename):
    model = joblib.load(pick_filename)
    return model


def prepareWaveforms(starttime, endtime, config):
    path = Path("tmp")
    path.mkdir(parents=True, exist_ok=True)

    if config["preprocess_data"]:
        for f in glob(os.path.join(path, "*")):
            os.remove(f)
    else:
        return True

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
        desc = "+++ Preparing raw data"
        for station in tqdm(stations, desc=desc, unit="station"):
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


def prepareInventory(config, proj, st, et, onsite=False):
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
    if onsite:
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
