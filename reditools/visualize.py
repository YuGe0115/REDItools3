"""Commandline tool for visualizing REDItools3 output."""

import argparse
import sys
from reditools.histogram import HistogramPlotter
from reditools.manhattan import ManhattanPlotter
from reditools.utils import load_data
from reditools.config import load_config

PLOTTERS = {
    "histogram": HistogramPlotter,
    "manhattan": ManhattanPlotter
}

def parse_options():
    """
    Parse commandline options for visualization.

    Returns:
        namespace: commandline args
    """
    parser = argparse.ArgumentParser(
        prog="reditools visualize",
        description="Visualize REDItools3 output as plots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to REDItools3 output TSV file",
    )
    parser.add_argument(
        "--plot-type",
        type=str,
        choices=["histogram", "manhattan"],
        required=True,
        help="Type of plot to generate",
    )
    parser.add_argument(
        "--output",
        default="output_plot",
        help="Output file prefix (without extension)",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to JSON config file",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=0,
        help="Minimum coverage threshold (Coverage-q30)",
    )
    parser.add_argument(
        "--subs-type",
        default=None,
        help="Substitution type to filter (e.g., AG)",
    )
    return parser.parse_args()

def main():
    """Generate visualizations for REDItools3 data."""
    options = parse_options()
    try:
        data = load_data(options.input)
        if options.min_coverage > 0:
            data = data[data["Coverage-q30"] >= options.min_coverage]
        if options.subs_type:
            data = data[data["AllSubs"].str.contains(options.subs_type, na=False)]
        config_data = load_config(options.config) if options.config else {}
        plotter_class = PLOTTERS.get(options.plot_type)
        if not plotter_class:
            sys.stderr.write(f"Error: Unsupported plot type: {options.plot_type}\n")
            sys.exit(1)
        plotter = plotter_class(data, options.output, config_data)
        plotter.plot()
    except Exception as e:
        sys.stderr.write(f"Error: Plotting failed: {e}\n")
        sys.exit(1)