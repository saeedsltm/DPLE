import os

import latlon as ll
from core.Extra import getStationMetaData, roundTo
from numpy import mean, sqrt
from obspy.geodetics.base import degrees2kilometers as d2k
from pandas import Series


def toSTATION0HYP(config, catalog):
    print("+++ Generating STATION0.HYP file ...")
    resetsPath = os.path.join("files", "resets.dat")
    station0hypPath = os.path.join("results", "STATION0.HYP")
    station_db = getStationMetaData(config["networks"],
                                    "*",
                                    config["starttime"],
                                    config["endtime"],
                                    *config["center"],
                                    config["maxradius"],
                                    catalog)
    velocities = config["vp"]
    depths = config["zz"]
    VpVs = config["vp_vs_ratio"]
    trialDepth = 10
    xNear = mean(station_db.apply(lambda x: mean(
        Series(sqrt((x.lon-station_db.lon)**2 + (x.lat-station_db.lat)**2))),
        axis=1))
    xNear = roundTo(d2k(xNear), base=5)
    xFar = 2.5*xNear
    stationLine = "  {code:4s}{latDeg:2.0f}{latMin:05.2f}N {lonDeg:2.0f}\
{lonMin:05.2f}E{elv:4.0f}\n"
    modelLine = " {v:5.2f}  {z:6.3f}             \n"
    controlLine = "{trialDepth:4.0f}.{xNear:4.0f}.{xFar:4.0f}. {VpVs:4.2f}"
    with open(resetsPath) as f, open(station0hypPath, "w") as g:
        for line in f:
            g.write(line)
        g.write("\n\n")
        for r, row in station_db.iterrows():
            code = row.code
            lat = ll.Latitude(row.lat)
            lon = ll.Longitude(row.lon)
            elv = row.elv
            g.write(stationLine.format(
                code=code,
                latDeg=lat.degree, latMin=lat.decimal_minute,
                lonDeg=lon.degree, lonMin=lon.decimal_minute,
                elv=elv
            ))
        g.write("\n")
        for v, z in zip(velocities, depths):
            g.write(modelLine.format(
                v=v, z=z
            ))
        g.write("\n")
        g.write(controlLine.format(
            trialDepth=trialDepth, xNear=xNear, xFar=xFar, VpVs=VpVs
        ))
        g.write("\nNew")
