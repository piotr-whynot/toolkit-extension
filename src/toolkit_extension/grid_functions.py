# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

import geopandas as gpd
import numpy as np
from cartopy.io.shapereader import natural_earth

# functions to support grid alignment
# Domain bounds mapping
# Each bounding box is in the form: [north, west, south, east]
DOMAIN_BOUNDS = {
        "westafrica":      [20, -20, 0, 25],
        "waf":             [20, -20, 0, 25],

        "centralafrica":   [15, 5, -15, 35],
        "caf":             [15, 5, -15, 35],

        "eastafrica":      [18, 22, -15, 52],
        "eaf":             [18, 22, -15, 52],

        "southernafrica":  [0, 10, -35, 42],
        "saf":             [0, 10, -35, 42],

        "northwestafrica": [38, -20, 12, 10],
        "nwa":             [38, -20, 12, 10],

        "northeastafrica": [38, 10, 12, 52],
        "nea":             [38, 10, 12, 52],

        "madagascar":      [-10, 42, -28, 54],
        "mdg":             [-10, 42, -28, 54],
}

def match_domain_from_bbox_bounds(bbox_bounds, tol=1e-6):
    """
    Check whether the provided bbox matches one of the predefined domains.
    """
    bbox = list(map(float, bbox_bounds))
    for name, bounds in DOMAIN_BOUNDS.items():
        bounds = list(map(float, bounds))
        if all(abs(a - b) <= tol for a, b in zip(bounds, bbox)):
            return name
    return None

def get_country_bbox_bounds(country_name, resolution="110m"):
    """
    Get bounding box for a country from Natural Earth.

    Parameters
    ----------
    country_name : str
        Country name as stored in Natural Earth 'ADMIN' field.
    resolution : str
        Natural Earth resolution. '110m' is sufficient for coarse grids.

    Returns
    -------
    list[float]
        Bounding box [north, west, south, east].
    """
    shpfilename = natural_earth(
        resolution=resolution,
            category="cultural",
            name="admin_0_countries"
    )

    world = gpd.read_file(shpfilename)

    country_name_norm = country_name.strip().casefold()
    admin_norm = world["ADMIN"].astype(str).str.strip().str.casefold()

    country = world[admin_norm == country_name_norm]

    if country.empty:
        raise ValueError(
            f"Country '{country_name}' not found in Natural Earth country boundaries."
        )

    minx, miny, maxx, maxy = country.total_bounds

    return [round(float(maxy), 2),
                round(float(minx), 2),
                round(float(miny), 2),
                round(float(maxx), 2)]

def format_coord(value, lat=True):
    """
    Format coordinate values for filename construction.
    """
    hemi = "N" if lat and value >= 0 else "S" if lat else "E" if value >= 0 else "W"
    return f"{abs(value)}{hemi}"

def snap_bbox_to_grid_containing(bbox_bounds, grid="1.5x1.5"):
    """
    Snap bbox [north, west, south, east] to the nearest valid bounds on the
    global grid such that the resulting bbox fully CONTAINS the original bbox.
    """
    north, west, south, east = map(float, bbox_bounds)

    grid_clean = grid.replace("/", "x")
    dlon, dlat = map(float, grid_clean.split("x"))

    snapped_north = np.ceil(north / dlat) * dlat
    snapped_south = np.floor(south / dlat) * dlat
    snapped_west = np.floor(west / dlon) * dlon
    snapped_east = np.ceil(east / dlon) * dlon

    snapped_north = min(90.0, snapped_north)
    snapped_south = max(-90.0, snapped_south)
    snapped_west = max(-180.0, snapped_west)
    snapped_east = min(180.0, snapped_east)

    return [
            round(snapped_north, 6),
            round(snapped_west, 6),
            round(snapped_south, 6),
            round(snapped_east, 6)
        ]

