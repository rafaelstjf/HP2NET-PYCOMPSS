"""PyCOMPSs entry point for the HP2NET workflow."""

import os

os.environ["HP2NET_RUNTIME"] = "pycompss"

from parsl_workflow import cli


if __name__ == "__main__":
    cli()
