import os
import warnings

from obspy import UTCDateTime as utc
from obspy import read_events
from obspy.core.event.origin import Pick
from obspy.clients.fdsn import Client
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
    weight_dict = {}
    print("+++ Loading selected catalog ...")
    catalog = read_events(catalogPath)
    catalog = catalog.filter(f"time > {starttime.isoformat()}",
                             f"time < {endtime.isoformat()}")
    desc = "+++ Calculating amplitude and magnitude using SeisComP ..."
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
            pick.waveform_id.channel_code = \
                inv[0].stations[0].channels[0].code[:2] + \
                pick.waveform_id.channel_code[-1]
            weight_dict[pick.resource_id] = pick.extra
    catalog.write("tmp.xml", format="sc3ml")
    cmd = f"scamp --ep tmp.xml -d {scdb} --reprocess --force > amp.xml"
    os.system(cmd)
    cmd = f"scmag --ep amp.xml -d {scdb} --reprocess > tmp.xml"
    os.system(cmd)
    print("+++ Loading QuickML catalog ...")
    catalog = read_events("tmp.xml")
    catalog = catalog.filter(f"time > {starttime.isoformat()}",
                             f"time < {endtime.isoformat()}")
    desc = "+++ Editing catalog for setting weights ..."
    for event in tqdm(catalog, desc=desc):
        preferred_origin = event.preferred_origin()
        picks = event.picks
        amplitudes = event.amplitudes
        for pick in picks:
            extra = weight_dict[pick.resource_id]
            pick.update({
                "extra": extra})
        for amplitude in amplitudes:
            amplitude.generic_amplitude *= (1e6/2080)  # count to Wood-Anderson nm
            amp_pick_id = amplitude.pick_id
            pick = [pick for pick in picks if pick.resource_id == amp_pick_id][0]
            new_amp_pick = Pick()
            new_amp_pick.time = pick.time
            new_amp_pick.waveform_id = pick.waveform_id
            new_amp_pick.onset = "impulsive"
            new_amp_pick.phase_hint = "AML"
            new_amp_pick.evaluation_mode = pick.evaluation_mode
            picks.append(new_amp_pick)
            amplitude.pick_id = new_amp_pick.resource_id
        event.picks = picks
        event.amplitudes = amplitudes

    catalog.write("fin.out", format="NORDIC", high_accuracy=False)
    os.rename("fin.out", catalogPath)
    os.remove("tmp.xml")
    os.remove("amp.xml")
