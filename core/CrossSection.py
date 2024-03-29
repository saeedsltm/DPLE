import math
import os
from string import ascii_uppercase

import numpy as np
import proplot as plt
from pandas import DataFrame, Series, read_csv
from pyproj import Proj


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Point(%g, %g, %g)" % (self.x, self.y, self.z)

    def __sub__(self, other):
        """P - Q"""
        if isinstance(other, Vector):
            return Point(self.x - other.x,
                         self.y - other.y,
                         self.z - other.z)

        return Vector(self.x - other.x,
                      self.y - other.y,
                      self.z - other.z)


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Vector(%g, %g, %g)" % (self.x, self.y, self.z)

    def norm(self):
        """u / |u|"""
        d = math.sqrt(self.x**2 + self.y**2 + self.z**2)

        return Vector(self.x / d, self.y / d, self.z / d)

    def __mul__(self, other):
        """dot product u · v or scaling x · u"""
        if isinstance(other, Vector):
            return (self.x * other.x
                    + self.y * other.y
                    + self.z * other.z)

        return Vector(self.x * other,
                      self.y * other,
                      self.z * other)


def cross(a, b):
    return Vector(a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x)


def getProfileDimension(x, y, z, length, theta):
    x_str, x_stp = x + 0.5*length * \
        np.cos(np.deg2rad(-theta+270)), x + 0.5*length*np.cos(np.deg2rad(-theta+90))
    y_str, y_stp = y + 0.5*length * \
        np.sin(np.deg2rad(-theta+270)), y + 0.5*length*np.sin(np.deg2rad(-theta+90))
    z_str, z_stp = 0, z
    if theta[0] > -90.0 and theta[0] < 0.0:
        x_str, x_stp = x_stp, x_str
        y_str, y_stp = y_stp, y_str
    return x_str[0], x_stp[0], y_str[0], y_stp[0], z_str, z_stp[0]


def distance(pointsA, pointB):
    distances = []
    for pointA in pointsA:
        distance = np.sqrt((pointA.x-pointB.x)**2 + (pointA.y-pointB.y)**2)
        distances.append(distance)
    return distances


def highlightArea(xmin, xmax, ymin, profile_w, profile_theta):
    x = np.linspace(xmin, xmax, 10)
    y = ymin + np.tan(np.deg2rad(profile_theta)) * x
    y_err = profile_w
    return x, y - y_err, y + y_err


def getEndSegments(xmin, xmax, ymin, ymax, width, theta):
    xminr = xmin + width*np.cos(np.deg2rad(-theta))
    yminr = ymin + width*np.sin(np.deg2rad(-theta))
    xminl = xmin + width*np.cos(np.deg2rad(-theta+180))
    yminl = ymin + width*np.sin(np.deg2rad(-theta+180))
    xmaxr = xmax + width*np.cos(np.deg2rad(-theta))
    ymaxr = ymax + width*np.sin(np.deg2rad(-theta))
    xmaxl = xmax + width*np.cos(np.deg2rad(-theta+180))
    ymaxl = ymax + width*np.sin(np.deg2rad(-theta+180))
    return xminr, yminr, xminl, yminl, xmaxr, ymaxr, xmaxl, ymaxl


def load_xyzm(locator, proj):
    xyzmFile = os.path.join("results", "location", locator, f"xyzm_{locator}.dat")
    xyzm_df = read_csv(xyzmFile, delim_whitespace=True)
    xyzm_df[[
        "x(km)",
        "y(km)"]] = xyzm_df.apply(lambda x: Series(
            proj(longitude=x["Lon"],
                 latitude=x["Lat"],
                 inverse=False)),
        axis=1)
    xyzm_df["z(km)"] = xyzm_df.Dep
    return xyzm_df


def make_profile(config, proj, dx, dy, length, width, z, theta):
    profile = DataFrame([{
        "Lon": config["center"][0],
        "Lat": config["center"][1],
        "Length": length,
        "Width": width,
        "Z": z,
        "Orientation": theta}])
    profile[[
        "X",
        "Y"]] = profile.apply(lambda x: Series(
            proj(longitude=x["Lon"],
                 latitude=x["Lat"],
                 inverse=False)),
        axis=1)
    profile.X += dx
    profile.Y += dy
    return profile


def getProfileGeometry(config, proj, Profile):
    dx = Profile["dx"]
    dy = Profile["dy"]
    length = Profile["length"]
    width = Profile["width"]
    z = Profile["z"]
    theta = Profile["theta"]
    profile = make_profile(config, proj, dx, dy, length, width, z, theta)
    xmin, xmax, ymin, ymax, zmin, zmax = getProfileDimension(profile.X,
                                                             profile.Y,
                                                             profile.Z,
                                                             profile.Length,
                                                             profile.Orientation)
    (xminr, yminr,
     xminl, yminl,
     xmaxr, ymaxr,
     xmaxl, ymaxl) = getEndSegments(xmin, xmax,
                                    ymin, ymax,
                                    width, theta)
    return (profile, xmin, xmax, ymin, ymax, zmin, zmax,
            xminr, yminr, xminl, yminl, xmaxr, ymaxr, xmaxl, ymaxl)


def plotCrossSection(config, locator):
    print(f"+++ Plotting cross sections for {locator} ...")
    center_lon, center_lat = config["center"]
    proj = Proj(f"+proj=sterea +lon_0={center_lon} +lat_0={center_lat} +units=km")
    xyzm_df = load_xyzm(locator, proj)

    for P, Profile in enumerate(config["profiles"]):
        cross_name = ascii_uppercase[P]
        outName = os.path.join("results", "location", locator,
                               f"cross_{cross_name}.png")
        (profile,
         xmin, xmax, ymin, ymax, zmin, zmax,
         xminr, yminr,
         xminl, yminl,
         xmaxr, ymaxr,
         xmaxl, ymaxl) = getProfileGeometry(config, proj, Profile)
        plane = (
            Point(xmin, ymin, zmin),
            Point(xmax, ymax, zmax),
            Point(profile.X[0], profile.Y[0], profile.Z[0])
        )

        points = [Point(x, y, z) for x, y, z in zip(xyzm_df["x(km)"],
                                                    xyzm_df["y(km)"],
                                                    xyzm_df["z(km)"])]

        x = plane[1] - plane[0]
        y = plane[2] - plane[0]
        u = cross(x, y).norm()

        newPoints = []
        for p in points:
            d = (p - plane[0]) * u
            pp = p - u * d
            if abs(d) > profile.Width[0]:
                continue
            newPoints.append(pp)

        plt.rc.update(
            {"fontsize": 9,
             "legend.fontsize": 6,
             "label.weight": "bold"})
        axShape = np.array([
            [0, 1, 0],
            [2, 2, 2]
        ])
        fig, axs = plt.subplots(axShape, share=False)
        axs.format(xlocator=("maxn", 4),
                   ylocator=("maxn", 4))
        [ax.set_aspect("equal") for ax in axs]
        [ax.grid(ls=":") for ax in axs]
        ax1 = axs[0]
        ax2 = axs[1]
        # shadows of other profiles
        for P_, Profile_ in enumerate(config["profiles"]):
            cross_name_ = ascii_uppercase[P_]
            (profile_,
             xmin_, xmax_, ymin_, ymax_, zmin_, zmax_,
             xminr_, yminr_,
             xminl_, yminl_,
             xmaxr_, ymaxr_,
             xmaxl_, ymaxl_) = getProfileGeometry(config, proj, Profile_)
            ax1.plot([xmin_, xmax_], [ymin_, ymax_], color="gray", ls="-", alpha=0.3)
            ax1.plot([xminr_, xminl_], [yminr_, yminl_],
                     color="gray", ls="-", alpha=0.3)
            ax1.plot([xmaxr_, xmaxl_], [ymaxr_, ymaxl_],
                     color="gray", ls="-", alpha=0.3)
            ax1.text(xmin_, ymin_, f"{cross_name_}",
                     rotation=profile_.Orientation[0],
                     border=True,
                     borderinvert=True,
                     borderwidth=1,
                     color="k",
                     alpha=0.3,
                     **{"weight": "bold", "size": "small", "ha": "center"})
            ax1.text(xmax_, ymax_, f"{cross_name_}'",
                     rotation=profile.Orientation[0],
                     border=True,
                     borderinvert=True,
                     borderwidth=1,
                     color="k",
                     alpha=0.3,
                     **{"weight": "bold", "size": "small", "ha": "center"})

        ax1.plot([xmin, xmax], [ymin, ymax], "r-")
        ax1.plot([xminr, xminl], [yminr, yminl], "r-")
        ax1.plot([xmaxr, xmaxl], [ymaxr, ymaxl], "r-")
        ax1.text(xmin, ymin, f"{cross_name}",
                 rotation=profile.Orientation[0],
                 border=True,
                 borderinvert=True,
                 borderwidth=1,
                 **{"weight": "bold", "size": "small", "ha": "center"})
        ax1.text(xmax, ymax, f"{cross_name}'",
                 rotation=profile.Orientation[0],
                 border=True,
                 borderinvert=True,
                 borderwidth=1,
                 **{"weight": "bold", "size": "small", "ha": "center"})
        x, y, z = [p.x for p in points], [p.y for p in points], [p.z for p in points]
        ax1.scatter(x, y, s=10, c=z,
                    mec="k",
                    mew=0.2,
                    cmap="Plasma_r")
        x, y, z = [p.x for p in newPoints], [
            p.y for p in newPoints], [p.z for p in newPoints]
        ax1.format(xlabel="Easting (km)",
                   ylabel="Northing (km)",
                   xreverse=False,
                   yreverse=False)

        distances = distance(
            newPoints,
            Point(xmin, ymin, zmin))
        sc = ax2.scatter(distances, z, s=10, c=z,
                         mec="k",
                         mew=0.2,
                         cmap="Plasma_r")
        dmax = distance(
            [Point(xmin, ymin, zmin)],
            Point(xmax, ymax, zmax))
        ax2.format(xlabel="Distance (km)",
                   ylabel="Depth (km)",
                   rtitle=f"N{profile.Orientation[0]:.0f} {cross_name}'",
                   ltitle=f"{cross_name}",
                   xlim=(0, dmax[0]),
                   ylim=(0, profile.Z[0]),
                   yreverse=True
                   )

        ax1.colorbar(sc, row=1, loc="r",
                     extend="both",
                     label="Depth (km)",
                     reverse=True)
        fig.save(outName, bbox_inches="tight")
