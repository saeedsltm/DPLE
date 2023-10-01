import os
import warnings

from obspy import UTCDateTime as utc
from obspy import read_events
from obspy.clients.fdsn import Client
from pandas import DataFrame, concat
from tqdm import tqdm

warnings.filterwarnings("ignore")


client = Client("http://localhost:8080")
default_weight = {
    "nordic_pick_weight": {
        "value": "0",
        "namespace": "https://seis.geus.net/software/seisan/node239.html"}}

scdb = "mysql://sysop:sysop@localhost/seiscomp"


def setMagnitude(catalogPath, config):
    starttime = utc(config["starttime"])
    endtime = utc(config["endtime"])
    inventory = client.get_stations(network="*",
                                    station="*",
                                    starttime=starttime,
                                    endtime=endtime,
                                    level="channel")
    weight_df = DataFrame([{"pick_id": None, "extra": default_weight}])
    catalog = read_events(catalogPath)
    catalog = catalog.filter(f"time > {starttime.isoformat()}",
                             f"time < {endtime.isoformat()}")
    desc = "+++ Editing catalog ..."
    for event in tqdm(catalog, desc=desc):
        preferred_origin = event.preferred_origin()
        arrivals = preferred_origin.arrivals
        for a, arrival in enumerate(arrivals):
            if not arrival.distance:
                arrivals.pop(a)
        for pick in event.picks:
            station = pick.waveform_id.station_code
            inv = inventory.select(
                station=station,
                starttime=preferred_origin.time,
                endtime=preferred_origin.time+1)
            network = inv[0].code
            pick.waveform_id.network_code = network
            pick.waveform_id.channel_code = inv[0].stations[0].channels[0].code[:2] + \
                pick.waveform_id.channel_code[-1]
            try:
                data = DataFrame(
                    [{"pick_id": pick.resource_id, "extra": pick.extra}])
                weight_df = concat([weight_df, data])
            except AttributeError:
                data = DataFrame([{"pick_id": None, "extra": default_weight}])
                weight_df = concat([weight_df, data])
    catalog.write("tmp.xml", format="sc3ml")
    cmd = f"scamp --ep tmp.xml -d {scdb} --reprocess --force > amp.xml"
    os.system(cmd)
    cmd = f"scmag --ep amp.xml -d {scdb} --reprocess > tmp.xml"
    os.system(cmd)
    catalog = read_events("tmp.xml")
    catalog = catalog.filter(f"time > {starttime.isoformat()}",
                             f"time < {endtime.isoformat()}")
    for event in catalog:
        preferred_origin = event.preferred_origin()
        for pick in event.picks:
            extra = weight_df[weight_df.pick_id == pick.resource_id].extra.values[0]
            pick.update({
                "extra": extra})
    catalog.write("fin.out", format="NORDIC", high_accuracy=False)
    os.rename("fin.out", catalogPath)
    os.remove("tmp.xml")
    os.remove("amp.xml")
