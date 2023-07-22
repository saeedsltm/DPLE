from datetime import timedelta as td

from obspy.clients.fdsn.mass_downloader import (CircularDomain, MassDownloader,
                                                Restrictions)
from pandas import date_range
from tqdm import tqdm


def fetchRawWaveforms(config):

    # Data Time Span
    startTime = config["starttime"]
    endTime = config["endtime"]

    startDateRange = date_range(startTime, endTime-td(days=1), freq="1D")
    endDateRange = date_range(startTime+td(days=1), endTime, freq="1D")

    # Data Region
    domain = CircularDomain(
        longitude=config["center"][0],
        latitude=config["center"][1],
        minradius=config["minradius"],
        maxradius=config["maxradius"])

    for st, et in tqdm(
            zip(startDateRange, endDateRange),
            desc="+++ Downloading waveforms"):

        # Data Restrictions
        restrictions = Restrictions(
            starttime=st,
            endtime=et,
            reject_channels_with_gaps=False,
            chunklength_in_sec=86400,
            minimum_length=0.0,
            network=config["network"],
            channel_priorities=["HH[ENZ]", "BH[ENZ]", "SH[ENZ]", "HN[ENZ]"],
            sanitize=True,
        )

        # Download Data
        mdl = MassDownloader(config["fdsn_urls"], configure_logging=False)
        mdl.download(
            domain,
            restrictions,
            mseed_storage=f"DB/{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}/waveforms",
            stationxml_storage=f"DB/{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}/stations")
