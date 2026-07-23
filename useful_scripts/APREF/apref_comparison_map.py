import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def make_map(df):
    """
    Generate a map of Asutralia showing Apref stations, highlighting stations with changes.

    :param df: pandas.DataFrame containing station information with columns 'lon', 'lat', and 'Category'.
    :return: None. Saves map as 'apref_station_map.png'.
    """

    data_crs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection=data_crs)

    # Extent for GDA2020
    ax.set_extent([99, 166, -44.5, -8], crs=data_crs)

    # Base map styling based on GA guide colours
    ocean_colour = (190/255, 232/255, 255/255)    # Ocean
    land_colour = (255/255, 218/255, 181/255)     # Landmass
    coast_colour = (0/255, 54/255, 255/255)       # Coastline
    state_colour = (78/255, 78/255, 78/255)       # State borders

    # Add basic features to map
    ax.set_facecolor(ocean_colour)

    ax.add_feature(cfeature.LAND, facecolor=land_colour, edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor=ocean_colour, edgecolor="none", zorder=0)
    ax.coastlines(resolution="10m", color=coast_colour, linewidth=0.6, zorder=3)
    ax.add_feature(cfeature.STATES, edgecolor=state_colour, linewidth=0.75, linestyle="-", facecolor="none", zorder=2)

    # Define styles for different categories of stations
    styles = {
        "Existing Station": {
            "marker": "o",
            "s": 6,
            "facecolor": "#434643",
            "edgecolor": "0.45",
            "linewidth": 0.4,
            "zorder": 3
        },
        "New Station": {
            "marker": "o",
            "s": 36,
            "facecolor": "#05B632",
            "edgecolor": "black",
            "linewidth": 0.4,
            "zorder": 5
        },
        "New Discontinuity": {
            "marker": "o",
            "s": 45,
            "facecolor": "#0077E6",
            "edgecolor": "black",
            "linewidth": 0.4,
            "zorder": 6
        },
        "Moved > 10 mm": {
            "marker": "o",
            "s": 42,
            "facecolor": "#D55E00",
            "edgecolor": "black",
            "linewidth": 0.4,
            "zorder": 7
        },
        "Removed Station": {
            "marker": "x",
            "s": 42,
            "facecolor": "#EE0B0B",
            "linewidth": 2,
            "zorder": 7
        }
    }

    # Plot stations on map based on category
    for category, style in styles.items():
        sub = df[df["Category"] == category]

        ax.scatter(
            sub["lon"],
            sub["lat"],
            transform=data_crs,
            label=category,
            **style
        )

    # Clean legend (avoid duplicates)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="lower left")

    plt.tight_layout()
    plt.savefig("apref_station_map.png", dpi=300, bbox_inches="tight")