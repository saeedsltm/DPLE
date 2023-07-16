from pandas import read_fwf


def getPick(picks, requestedPickID):
    for pick in picks:
        if pick.resource_id == requestedPickID:
            return pick


def roundTo(x, base=5):
    return base * round(x/base)


def loadhypo71Out():
    names = ["yy", "mo", "dd", "A", "hh", "mm", "B", "sssss", "C",
             "yd", "D", "ymmmm", "E", "xdd", "F", "xmmmm", "G",
             "depth_", "HHHH", "mag", "L", "ns", "M", "gap", "N", "dmin", "O",
             "rms_", "erh__", "erz__", "P", "qm"]
    widths = [len(name) for name in names]
    bulletin_df = read_fwf("results/hyp71.out", names=names,
                           widths=widths, header=0)
    return bulletin_df
