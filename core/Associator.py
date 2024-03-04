from gamma.utils import association
from pandas import DataFrame, Series, date_range
from tqdm import tqdm
import os
import pyocto
from obspy import read_inventory
from core.PrepareData import prepareInventory, loadModel
from numpy import array, nan
from pathlib import Path
from datetime import timedelta as td
from pyproj import Proj
from gamma.utils import estimate_eps
from obspy.geodetics.base import degrees2kilometers as d2k


def applyGaMMaConfig(config, stations):
    method = {
        "BGMM": 5,
        "GMM": 1}
    config["oversample_factor"] = method[config["method"]]
    config["ncpu"] = os.cpu_count() - 2
    config["degree2km"] = 111.19492474777779
    config["h"] = 3.0
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
        config["dbscan_eps"] = estimate_eps(stations, config["vel"]["p"])
    return config


def assocciateWithGamma(config, st, et, proj):
    picks = loadModel(os.path.join(
        "results", f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.jlib"))
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
    catalogs, assignments = association(
        pick_df,
        station_df,
        config,
        event_index0,
        config["method"],
        pbar=pbar)
    event_index0 += len(catalogs)

    if event_index0 == 0:
        return False

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
    catalogs["depth"] = catalogs["z(km)"]
    catalogs.replace({"magnitude": 999}, 99, inplace=True)
    catalogs["magnitude"] = nan
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
                "depth",
                "sigma_time",
                "sigma_amp",
                "cov_time_amp",
                "event_index",
                "gamma_score"])

    # Add assignment to picks
    assignments = DataFrame(
        assignments,
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
    pick_df["time"] = pick_df["timestamp"]
    pick_df["phase"] = pick_df["type"]
    pick_df["probability"] = pick_df["prob"]
    pick_df["station"] = pick_df["id"]
    pick_df["event_idx"] = pick_df["event_index"]
    with open(picks_csv, "w") as fp:
        pick_df.to_csv(
            fp,
            sep="\t",
            index=False,
            date_format='%Y-%m-%dT%H:%M:%S.%f',
            columns=[
                "station",
                "time",
                "phase",
                "probability",
                "phase_amp",
                "event_idx",
                "gamma_score"])
    return True


def assocciateWithPyocto(config, st, et, proj):
    picks = loadModel(os.path.join(
        "results", f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.jlib"))
    stationxml = os.path.join(
        "DB",
        f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}",
        "stations",
        "*.xml")
    inv = read_inventory(stationxml)
    association_cutoff_distance = config["association_cutoff_distance"]
    time_before = config["time_before"]
    n_picks = config["n_picks"]
    n_p_picks = config["n_p_picks"]
    n_s_picks = config["n_s_picks"]
    n_p_and_s_picks = config["n_p_and_s_picks"]
    vp = array(config["vp"])
    z = array(config["zz"])
    vs = vp/config["vp_vs_ratio"]
    model = {"vp": vp,
             "vs": vs,
             "depth": z}
    model = DataFrame(model)
    path = os.path.join("results", "vm.dat")
    maxradius = d2k(config["maxradius"])
    maxdepth = config["zlim_degree"][1]
    pyocto.VelocityModel1D.create_model(model, 1.0, maxradius, maxdepth, path)
    tolerance = config["tolerance"]
    association_cutoff_distance = config["association_cutoff_distance"]
    velocity_model = pyocto.VelocityModel1D(path,
                                            tolerance,
                                            association_cutoff_distance,
                                            surface_p_velocity=vp[0],
                                            surface_s_velocity=vs[0])
    associator = pyocto.OctoAssociator.from_area(
        lat=config["ylim_degree"],
        lon=config["xlim_degree"],
        zlim=config["zlim_degree"],
        time_before=time_before,
        velocity_model=velocity_model,
        n_picks=n_picks,
        n_p_picks=n_p_picks,
        n_s_picks=n_s_picks,
        n_p_and_s_picks=n_p_and_s_picks,
    )
    stations = associator.inventory_to_df(inv)
    events, assignments = associator.associate_seisbench(picks, stations)
    associator.transform_events(events)
    if len(events) == 0:
        return False
    # save catalog
    catalog_csv = os.path.join(
        "results",
        f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
    events["magnitude"] = nan
    with open(catalog_csv, "w") as fp:
        events.to_csv(
            fp,
            sep="\t",
            index=False,
            float_format="%.3f",
            date_format='%Y-%m-%dT%H:%M:%S.%f',
            columns=[
                "time",
                "longitude",
                "latitude",
                "depth",
                "magnitude",
                "picks",
                "x",
                "y",
                "z",
                "idx"])
    # save picks
    picks_csv = os.path.join(
        "results",
        f"picks_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv")
    assignments["amplitude"] = nan
    with open(picks_csv, "w") as fp:
        assignments.to_csv(
            fp,
            sep="\t",
            index=False,
            date_format='%Y-%m-%dT%H:%M:%S.%f',
            columns=[
                "event_idx",
                "pick_idx",
                "residual",
                "station",
                "time",
                "probability",
                "phase",
                "amplitude"])
    return True


def runAssociator(config):
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
    for st, et in zip(startDateRange, endDateRange):
        print(f"+++ Run {config['associator']} Associator on period: {st} - {et}")
        if config["associator"] == "GaMMA":
            assocciateWithGamma(config, st, et, proj)
        elif config["associator"] == "PyOcto":
            assocciateWithPyocto(config, st, et, proj)
