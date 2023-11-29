import os
from datetime import timedelta as td
from pathlib import Path
from glob import glob

import seisbench.models as sbm
from gamma.utils import association
from obspy import read
from pandas import DataFrame, Series, date_range, concat, read_csv
from pyproj import Proj
from tqdm import tqdm
from obspy.core.stream import Stream

from core.PrepareData import (applyGaMMaConfig, picks2DF, prepareInventory,
                              prepareWaveforms)
from core.Extra import divide_chunks


def runSeisBench(config):
    path = Path("results")
    path.mkdir(parents=True, exist_ok=True)

    startTime = config["starttime"]
    endTime = config["endtime"]

    startDateRange = date_range(startTime, endTime-td(days=1), freq="1D")
    endDateRange = date_range(startTime+td(days=1), endTime, freq="1D")

    proj = Proj(f"+proj=sterea\
                +lon_0={config['center'][0]}\
                +lat_0={config['center'][1]}\
                +units=km")

    # Loop over one day data
    for st, et in zip(startDateRange, endDateRange):
        print(f"\n+++ Working on period: {st} - {et}")

        # Prepare One-day-length data
        dataExists = prepareWaveforms(st, et, config)
        if not dataExists:
            continue

        chunksData = glob(os.path.join("tmp", "*.mseed"))
        for c, chunkData in enumerate(divide_chunks(chunksData, 10)):

            pick_outfile = f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}_{c}.csv"

            if not config["repick_data"] and os.path.exists(
                    os.path.join("results", pick_outfile)):
                continue

            stream = Stream()
            for s in chunkData:
                stream += read(s)

            min_p_prob = config["min_p_prob"]
            min_s_prob = config["min_s_prob"]

            # Apply PhaseNet predict method
            print("+++ Applying SeisBench ...")
            ncpu = os.cpu_count() - 2
            picker = sbm.PhaseNet.from_pretrained(config["model"])
            picks = picker.classify(stream,
                                    batch_size=64,
                                    P_threshold=min_p_prob,
                                    S_threshold=min_s_prob,
                                    parallelism=ncpu).picks
            if config["repick_data"] and os.path.exists(
                    os.path.join("results", f"{pick_outfile}.csv")):
                os.remove(os.path.join("results", f"{pick_outfile}.csv"))

            if len(picks):
                picks2DF(picks, pick_outfile)

        # Create DataFrame for stations and picks
        picks_df_list = glob(os.path.join(
            "results", f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}_*.csv"))
        pick_df = concat(map(read_csv, picks_df_list))
        pick_df.reset_index(inplace=True, drop=True)
        pick_df.to_csv(os.path.join(
            "results", f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv"))
        station_df, station_dict = prepareInventory(config, proj, st, et)

        # Apply GaMMa configuration
        config = applyGaMMaConfig(config, station_df)

        # Removes picks without amplitude if amplitude flag is set to True
        if config["use_amplitude"]:
            pick_df = pick_df[pick_df["phase_amplitude"] != -1]

        # Rum GaMMa associator
        event_index0 = 0
        assignments = []
        pbar = tqdm(1)
        for f in [f"assignments_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv",
                  f"picks_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv",
                  f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv",
                  f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.out"]:
            if os.path.exists(os.path.join("results", f)):
                os.remove(os.path.join("results", f))
        catalogs, assignments = association(
            pick_df,
            station_df,
            config,
            event_index0,
            config["method"],
            pbar=pbar)
        event_index0 += len(catalogs)

        if event_index0 == 0:
            continue

        # Create catalog
        catalog_csv = os.path.join(
            "results",
            f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        catalogs = DataFrame(
            catalogs,
            columns=["time"]+config["dims"]+[
                "magnitude",
                "sigma_time",
                "sigma_amp",
                "cov_time_amp",
                "event_index",
                "gamma_score"])
        catalogs[[
            "longitude",
            "latitude"]] = catalogs.apply(lambda x: Series(
                proj(longitude=x["x(km)"],
                     latitude=x["y(km)"],
                     inverse=True)),
            axis=1)
        catalogs["depth(m)"] = catalogs["z(km)"].apply(lambda x: x*1e3)
        catalogs.replace({"magnitude": 999}, 99, inplace=True)
        with open(catalog_csv, "w") as fp:
            catalogs.to_csv(
                fp,
                sep="\t",
                index=False,
                float_format="%.3f",
                date_format='%Y-%m-%dT%H:%M:%S.%f',
                columns=[
                    "time",
                    "magnitude",
                    "longitude",
                    "latitude",
                    "depth(m)",
                    "sigma_time",
                    "sigma_amp",
                    "cov_time_amp",
                    "event_index",
                    "gamma_score"])

        # Add assignment to picks
        assignments_csv = os.path.join(
            "results",
            f"assignments_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        assignments = DataFrame(
            assignments,
            columns=[
                "pick_index",
                "event_index",
                "gamma_score"])
        with open(assignments_csv, "w") as fp:
            assignments.to_csv(
                fp,
                sep="\t",
                index=False,
                columns=[
                    "pick_index",
                    "event_index",
                    "gamma_score"])
        picks_csv = os.path.join(
            "results",
            f"picks_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        pick_df = pick_df.join(
            assignments.set_index("pick_index")
        ).fillna(-1).astype({'event_index': int})
        with open(picks_csv, "w") as fp:
            pick_df.to_csv(
                fp,
                sep="\t",
                index=False,
                date_format='%Y-%m-%dT%H:%M:%S.%f',
                columns=[
                    "id",
                    "timestamp",
                    "type",
                    "prob",
                    "phase_amp",
                    "event_index",
                    "gamma_score"])
