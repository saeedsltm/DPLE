import os
from datetime import timedelta as td

import proplot as plt
from numpy import (abs, arange, array, histogram, max, mean, min, random, sqrt,
                   unique)
from obspy import Stream
from obspy import UTCDateTime as utc
from obspy import read, read_events
from obspy.geodetics.base import degrees2kilometers as d2k
from pandas import DataFrame, Series, date_range, read_csv
from pyproj import Proj
from tqdm import tqdm

from core.Extra import handle_masked_arr, weighted_avg_and_std, weightMapper
from core.PrepareData import prepareInventory
from pathlib import Path


def plotSeismicity(config):
    path = Path(os.path.join("results", "figures"))
    path.mkdir(parents=True, exist_ok=True)
    startTime = config["starttime"]
    endTime = config["endtime"]

    startDateRange = date_range(startTime, endTime-td(days=1), freq="1D")
    endDateRange = date_range(startTime+td(days=1), endTime, freq="1D")

    if config["plotResults"]:
        proj = Proj(f"+proj=sterea\
                    +lon_0={config['center'][0]}\
                    +lat_0={config['center'][1]}\
                    +units=km")
        for st, et in tqdm(
                zip(startDateRange, endDateRange),
                desc="+++ Plotting seismicity maps"):
            catalog = os.path.join(
                "results",
                f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
            if not os.path.exists(catalog):
                continue
            catalog = read_csv(catalog, sep="\t")
            if len(catalog) <= 1:
                continue
            station_df, station_dict = prepareInventory(config, proj, st, et)
            catalog[[
                "x(km)",
                "y(km)"]] = catalog.apply(lambda x: Series(
                    proj(longitude=x["longitude"],
                         latitude=x["latitude"],
                         inverse=False)),
                axis=1)
            catalog["z(km)"] = catalog["depth(m)"]*1e-3

            fig, axs = plt.subplots()
            [ax.grid(ls=":") for ax in axs]
            ax = axs[0]
            ax.format(
                ultitle=f"{len(catalog)} events from {st.strftime('%Y-%m-%d')} to {et.strftime('%Y-%m-%d')}",
                fontsize=5)
            ax.set_aspect("equal")
            cb = ax.scatter(
                catalog["x(km)"],
                catalog["y(km)"],
                c=catalog["z(km)"],
                s=catalog["magnitude"] if catalog["magnitude"].mean() != 99 else 10,
                cmap="inferno_r",
                mew=0.25, mec="k", mfc="r",
                vmin=0)
            cbar = fig.colorbar(cb)
            cbar.ax.set_ylim(cbar.ax.get_ylim()[::-1])
            cbar.set_label("Depth[km]")
            ax.scatter(
                station_df["x(km)"],
                station_df["y(km)"],
                m="^", s=10, mew=0.25, mec="k", mfc="gray", alpha=0.5)
            for x, y, s in zip(
                    station_df["x(km)"],
                    station_df["y(km)"],
                    station_df["id"]):
                ax.text(x, y, s.split(".")[1], fontsize=3, alpha=0.5)
            ax.set_xlabel("Easting [km]")
            ax.set_ylabel("Northing [km]")
            fig.save(os.path.join(
                path,
                f"seismicity_{st.strftime('%Y%m%d')}_{st.strftime('%Y%m%d')}.png"))


def pickerTest(config):
    path = Path(os.path.join("results", "figures"))
    path.mkdir(parents=True, exist_ok=True)
    startTime = config["starttime"]
    endTime = config["endtime"]

    startDateRange = date_range(startTime, endTime-td(days=1), freq="1D")
    endDateRange = date_range(startTime+td(days=1), endTime, freq="1D")

    for st, et in tqdm(
            zip(startDateRange, endDateRange),
            desc="+++ Plotting picker test samples"):

        catalog = os.path.join(
            "results",
            f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        pick = os.path.join(
            "results",
            f"picks_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        assignment = os.path.join(
            "results",
            f"assignments_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
        proj = Proj(f"+proj=sterea\
                    +lon_0={config['center'][0]}\
                    +lat_0={config['center'][1]}\
                    +units=km")

        if not os.path.exists(catalog) and not os.path.exists(pick):
            continue

        catalog_df = read_csv(catalog, sep="\t")
        catalog_df.sort_values(by=["time"], inplace=True)
        pick_df = read_csv(pick, sep="\t")
        assignment_df = read_csv(assignment, sep="\t")
        station_df, station_dict = prepareInventory(config, proj, st, et)
        catalog_df[[
            "x(km)",
            "y(km)"]] = catalog_df.apply(lambda x: Series(
                proj(longitude=x["longitude"],
                     latitude=x["latitude"],
                     inverse=False)),
            axis=1)
        catalog_df["z(km)"] = catalog_df["depth(m)"]*1e-3

        for n in range(config["nTests"]):

            event_index = random.randint(len(catalog_df))
            event_picks = [pick_df.iloc[i]
                           for i in assignment_df[assignment_df["event_index"] == event_index]["pick_index"]]
            event = catalog_df.iloc[event_index]

            times = [utc(pick.timestamp) for pick in event_picks]

            if len(times) == 0:
                continue

            first, last = min(times), max(times)

            sub = Stream()

            for station in unique(
                    [pick.id.split(".")[1] for pick in event_picks]):
                streamFile = os.path.join(
                    "DB",
                    f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}",
                    "waveforms",
                    f"??.{station}.*.???__{st.strftime('%Y%m%d')}T000000Z__{et.strftime('%Y%m%d')}T000000Z.mseed")
                try:
                    stream = read(streamFile)
                except Exception:
                    continue
                stream.merge(fill_value=None)
                stream = handle_masked_arr(stream)
                streamZ = stream.select(station=station, channel="??Z")
                if len(streamZ):
                    sub.append(streamZ[0])
            if len(sub) == 0:
                continue
            sub = sub.slice(first - 5, last + 5)

            sub = sub.copy()
            sub.detrend()
            sub.filter("highpass", freq=1.0)

            fig, axs = plt.subplots()
            ax = axs[0]
            [ax.grid(ls=":") for ax in axs]
            ax.format(
                ultitle=f"Ort={utc(event['time']).strftime('%Y-%m-%dT%H:%M:%S')}, Lon={event['longitude']:0.3f}, Lat={event['latitude']:0.3f}, Dep={event['depth(m)']*1e-3:0.3f}, Mag={event['magnitude']:0.1f}",
                fontsize=4)

            for i, trace in enumerate(sub):
                normed = trace.data - mean(trace.data)
                normed = normed / max(abs(normed))
                station_x, station_y = station_dict[trace.id[:-4]]
                y = sqrt((station_x - event["x(km)"]) ** 2 +
                         (station_y - event["y(km)"]) ** 2 +
                         event["z(km)"] ** 2)
                ax.plot(trace.times(), 5 * normed + y, lw=0.25)

            for pick in event_picks:
                station_x, station_y = station_dict[pick.id]
                y = sqrt((station_x - event["x(km)"]) ** 2 +
                         (station_y - event["y(km)"]) ** 2 +
                         event["z(km)"] ** 2)
                x = utc(pick.timestamp) - trace.stats.starttime
                if pick.type.upper() == "P":
                    ls = '-'
                else:
                    ls = '--'
                ax.plot([x, x], [y - 3, y + 3], 'k', ls=ls, lw=0.5)
                ax.text(x, y+3, f"{pick.prob:0.1f}", fontsize=3)

            ax.set_ylim(0)
            ax.set_xlim(0, max(trace.times()))
            ax.set_ylabel("Hypocentral distance [km]")
            ax.set_xlabel("Time [s]")

            fig.save(os.path.join(
                path,
                f"pickerTest_{n}_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.png"))


def pickerStats(config):
    path = Path(os.path.join("results", "figures"))
    path.mkdir(parents=True, exist_ok=True)
    for hypFile in ["hyp", "select"]:
        catFile = os.path.join("results", "location", "hypocenter", f"{hypFile}.out")
        if not os.path.exists(catFile):
            return

        print(f"+++ loading catalog from {hypFile} file ...")
        catalog = read_events(catFile)
        data = {}

        # Statistical metrics
        dist_P, ttim_P, resi_P, wegt_P = [], [], [], []
        dist_S, ttim_S, resi_S, wegt_S = [], [], [], []

        # Loop over catalog
        desc = "+++ Exporting catalog to DataFrame"
        for event in tqdm(catalog, desc=desc):
            preferred_origin = event.preferred_origin()
            picks = event.picks
            arrivals = preferred_origin.arrivals
            for arrival in arrivals:
                pha = arrival.phase
                sta = [
                    pick.waveform_id.station_code
                    for pick in picks if pick.resource_id ==
                    arrival.pick_id][0]
                art = [
                    pick.time for pick in picks if pick.resource_id ==
                    arrival.pick_id][0]
                wet = [
                    int(pick.extra.get("nordic_pick_weight")["value"])
                    for pick in picks if pick.resource_id ==
                    arrival.pick_id][0]
                ttm = art - preferred_origin.time
                if abs(ttm) > 1e3:
                    continue
                if sta not in data:
                    data[sta] = {
                        "DIST_P": [],
                        "TTIM_P": [],
                        "RESI_P": [],
                        "WEGT_P": [],
                        "DIST_S": [],
                        "TTIM_S": [],
                        "RESI_S": [],
                        "WEGT_S": [],
                    }
                if "P" in pha.upper():
                    dist_P.append(arrival.distance)
                    ttim_P.append(ttm)
                    resi_P.append(arrival.time_residual)
                    wegt_P.append(wet)
                    data[sta]["DIST_P"].append(arrival.distance)
                    data[sta]["TTIM_P"].append(ttm)
                    data[sta]["RESI_P"].append(arrival.time_residual)
                    data[sta]["WEGT_P"].append(wet)
                if "S" in pha.upper():
                    dist_S.append(arrival.distance)
                    ttim_S.append(ttm)
                    resi_S.append(arrival.time_residual)
                    wegt_S.append(wet)
                    data[sta]["DIST_S"].append(arrival.distance)
                    data[sta]["TTIM_S"].append(ttm)
                    data[sta]["RESI_S"].append(arrival.time_residual)
                    data[sta]["WEGT_S"].append(wet)

        df_P = DataFrame(
            {
                "DIST_P": dist_P, "TTIM_P": ttim_P,
                "RESI_P": resi_P, "WEGT_P": wegt_P
            }
        )
        df_S = DataFrame(
            {
                "DIST_S": dist_S, "TTIM_S": ttim_S,
                "RESI_S": resi_S, "WEGT_S": wegt_S
            }
        )

        df_P["DIST_P"] = d2k(df_P["DIST_P"])
        df_S["DIST_S"] = d2k(df_S["DIST_S"])

        # Travel times curve
        print(f"+++ Plot travel time curve using {hypFile} file ...")
        fig, axs = plt.subplots()
        [ax.grid(ls=":") for ax in axs]
        ax = axs[0]
        ax.format(
            xlabel="Distance (km)",
            ylabel="Travel time (s)",
            xlim=(0, config["maxDist"]))

        p = ax.scatter(
            df_P["DIST_P"], df_P["TTIM_P"],
            m="^", c=df_P["RESI_P"], s=5, cmap="rdylbu_r",
            vmin=-config["minmaxRes"], vmax=config["minmaxRes"],
            ec="k", ew=0.2)
        s = ax.scatter(
            df_S["DIST_S"], df_S["TTIM_S"],
            m="s", c=df_S["RESI_S"], s=5, cmap="rdylbu_r",
            vmin=-config["minmaxRes"], vmax=config["minmaxRes"],
            mec="k", mew=0.2)

        ax.colorbar(p, loc="r", label="Residuals (s)")
        fig.save(os.path.join(path, f"traveltime_{hypFile}.png"))

        # Statistics figures
        desc = f"+++ Plotting statistics using {hypFile} file ..."
        for station in tqdm(data.keys(), desc=desc):

            db_P = DataFrame(
                {
                    "DIST_P": data[station]["DIST_P"],
                    "TTIM_P": data[station]["TTIM_P"],
                    "RESI_P": data[station]["RESI_P"],
                    "WEGT_P": data[station]["WEGT_P"]
                }
            )
            db_S = DataFrame(
                {
                    "DIST_S": data[station]["DIST_S"],
                    "TTIM_S": data[station]["TTIM_S"],
                    "RESI_S": data[station]["RESI_S"],
                    "WEGT_S": data[station]["WEGT_S"]
                }
            )

            db_P.dropna(inplace=True)
            db_S.dropna(inplace=True)

            W_M_RESI_P, W_STD_RESI_P = weighted_avg_and_std(
                db_P["RESI_P"], weightMapper(db_P["WEGT_P"], reverse=True))
            W_M_RESI_S, W_STD_RESI_S = weighted_avg_and_std(
                db_S["RESI_S"], weightMapper(db_S["WEGT_S"], reverse=True))

            axShape = array(
                [[1, 2]]
            )
            fig, axs = plt.subplots(axShape, sharey=False, sharex=True)
            [ax.grid(ls=":") for ax in axs]
            axs.format(
                suptitle=station,
                xlabel="Time residues (s)",
                ylabel="Number of picks (#)",
                xlocator=("maxn", 5),
                xlim=(-config["minmaxRes"], config["minmaxRes"])
            )

            W = [0, 1, 2, 3, 4]
            C = ["gray2", "gray4", "gray6", "gray8", "red9"]
            dr = 0.1
            bins = arange(-config["minmaxRes"], config["minmaxRes"] + dr, dr)
            areas = []
            for d, l, ax in zip([db_P, db_S], ["P", "S"], axs):
                if l == "P":
                    s = "$\overline{m}_w=$"+f"{W_M_RESI_P:0.2f}" + \
                        ", $\mu_w=$"+f"{W_STD_RESI_P:0.2f}"
                if l == "S":
                    s = "$\overline{m}_w=$"+f"{W_M_RESI_S:0.2f}" + \
                        ", $\mu_w=$"+f"{W_STD_RESI_S:0.2f}"
                ax.format(
                    ultitle=s,
                    urtitle=l)
                for w, c in zip(W, C):
                    df = d[d[f"WEGT_{l}"] == w]
                    h, edges = histogram(df[f"RESI_{l}"], bins=bins)
                    x = [mean([i, j]) for i, j in zip(edges[:-1], edges[1:])]
                    area = ax.area(x, h, color=c, lw=0.5, ec="k", label=str(w))
                    if l == "P":
                        areas.append(area)
            fig.legend(areas, loc="r", title="Weights", ncols=1)
            fig.save(os.path.join(path, f"{station}_stat_{hypFile}.png"))
