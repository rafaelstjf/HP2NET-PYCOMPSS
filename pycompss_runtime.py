"""Small PyCOMPSs adapter used by HP2NET application tasks."""

from __future__ import annotations

import functools
import inspect
import os
import subprocess
from typing import Callable

try:
    from pycompss.api.parameter import COLLECTION_IN
    from pycompss.api.task import task
except ImportError as exc:  
    raise RuntimeError(
        "PyCOMPSs is required. Start this workflow with 'runcompss', not "
        "with the regular Python interpreter."
    ) from exc


AUTO_LOGNAME = None


def python_task(*_args, **_kwargs):

    def decorator(function: Callable):
        options = {"returns": 1}
        if "inputs" in inspect.signature(function).parameters:
            options["inputs"] = COLLECTION_IN
        return task(**options)(function)

    return decorator


def bash_task(*_args, **_kwargs):

    def decorator(command_builder: Callable):
        @task(returns=1, inputs=COLLECTION_IN)
        @functools.wraps(command_builder)
        def run_command(*args, **kwargs):
            stdout_path = kwargs.pop("stdout", None)
            stderr_path = kwargs.pop("stderr", None)
            command = command_builder(*args, **kwargs)
            if not command:
                return None

            stdout_handle = open(stdout_path, "a") if isinstance(stdout_path, str) else None
            stderr_handle = open(stderr_path, "a") if isinstance(stderr_path, str) else None
            try:
                try:
                    completed = subprocess.run(
                        command,
                        shell=True,
                        check=True,
                        stdout=stdout_handle,
                        stderr=(stderr_handle if stderr_handle is not None else subprocess.PIPE),
                        text=True,
                        errors="replace",
                        executable=None if os.name == "nt" else "/bin/bash",
                    )
                    return completed.returncode
                except subprocess.CalledProcessError as exc:
                    details = (exc.stderr or "").strip()
                    message = f"External command failed with exit code {exc.returncode}: {command}"
                    if details:
                        message += f"\n--- command stderr ---\n{details}"
                    raise RuntimeError(message) from exc
            finally:
                if stdout_handle is not None:
                    stdout_handle.close()
                if stderr_handle is not None:
                    stderr_handle.close()

        return run_command

    return decorator
