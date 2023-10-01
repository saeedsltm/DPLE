from numpy import sqrt
from pandas import read_csv, to_datetime


def loadHypoDDRelocFile():
    names = ["ID",  "LAT",  "LON",  "DEPTH",
             "X",  "Y",  "Z",
             "EX",  "EY",  "EZ",
             "YR",  "MO",  "DY",  "HR",  "MI",  "SC",
             "MAG",
             "NCCP",  "NCCS",
             "NCTP",  "NCTS",
             "RCC",  "RCT",
             "CID "]
    hypodd_df = read_csv("hypoDD.reloc",
                         delim_whitespace=True,
                         names=names,
                         na_values="********")
    return hypodd_df


def writexyzm(outName):
    hypodd_df = loadHypoDDRelocFile()
    outputFile = f"xyzm_{outName}.dat"
    hypodd_df["year"] = hypodd_df.YR
    hypodd_df["month"] = hypodd_df.MO
    hypodd_df["day"] = hypodd_df.DY
    hypodd_df["hour"] = hypodd_df.HR
    hypodd_df["minute"] = hypodd_df.MI
    hypodd_df["second"] = hypodd_df.SC
    hypodd_df["ORT"] = to_datetime(hypodd_df[["year",
                                              "month",
                                              "day",
                                              "hour",
                                              "minute",
                                              "second"]])
    hypodd_df["ORT"] = hypodd_df["ORT"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    hypodd_df["Lon"] = hypodd_df.LON
    hypodd_df["Lat"] = hypodd_df.LAT
    hypodd_df["Dep"] = hypodd_df.DEPTH
    hypodd_df["Mag"] = hypodd_df.MAG
    hypodd_df["Nus"] = hypodd_df.NCTP
    hypodd_df["NuP"] = hypodd_df.NCTP
    hypodd_df["NuS"] = hypodd_df.NCTS
    hypodd_df["ADS"] = 0
    hypodd_df["MDS"] = 0
    hypodd_df["GAP"] = 0
    hypodd_df["RMS"] = hypodd_df.RCT
    hypodd_df["ERH"] = sqrt(hypodd_df.EX**2 + hypodd_df.EY**2)*1e-3
    hypodd_df["ERZ"] = hypodd_df.EZ*1e-3
    columns = ["ORT", "Lon", "Lat", "Dep", "Mag",
               "Nus", "NuP", "NuS", "ADS", "MDS", "GAP", "RMS", "ERH", "ERZ"]
    with open(outputFile, "w") as f:
        hypodd_df.to_string(f, columns=columns, index=False, float_format="%7.3f")
