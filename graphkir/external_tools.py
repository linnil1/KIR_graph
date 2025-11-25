"""
External tools runner for different engines (docker, podman, singularity, local)
"""

import uuid
import subprocess
from typing import Callable
from dataclasses import dataclass
from .utils import runShell


def _convertArgsForSingularity(args: list[str]) -> list[str]:
    """Convert docker-style arguments to singularity format"""
    converted = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "-e":
            converted.extend(["--env", args[i + 1]])
            i += 2
        elif arg == "-v":
            mount = args[i + 1].replace(":/app/", f":$PWD/")
            converted.extend(["--bind", mount])
            i += 2
        else:
            converted.append(arg)
            i += 1
    return converted


@dataclass
class EngineConfig:
    name: str
    path: str
    run_args: list[str]
    image_prefix: str = ""
    name_flag: str | None = None
    convert_args_func: Callable[[list[str]], list[str]] | None = None


# Global engine configuration
_engine = "podman"
_engine_configs: dict[str, EngineConfig] = {
    "podman": EngineConfig(
        name="podman",
        path="podman",
        run_args=["run", "-t", "--rm", "-u", "root", "-w", "/app", "-v", "$PWD:/app"],
        name_flag="--name",
    ),
    "docker": EngineConfig(
        name="docker",
        path="/usr/bin/docker",
        run_args=["run", "-t", "--rm", "-u", "root", "-w", "/app", "-v", "$PWD:/app"],
        name_flag="--name",
    ),
    "singularity": EngineConfig(
        name="singularity",
        path="singularity",
        run_args=["run", "-B", "$PWD"],
        image_prefix="docker://",
        convert_args_func=_convertArgsForSingularity,
    ),
    "local": EngineConfig(
        name="local",
        path="",
        run_args=[],
    ),
}

# Tool image mappings
_images = {
    "samtools": "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    "clustalo": "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    "hisat": "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle": "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
    "bwa": "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
}


def setEngine(engine: str) -> None:
    """Set the engine to use"""
    global _engine
    if engine not in _engine_configs:
        raise ValueError(
            f"Unsupported engine: {engine}. Supported: {list(_engine_configs.keys())}"
        )
    _engine = engine


def getEngine() -> str:
    """Get current engine"""
    return _engine


def addCustomEngine(config: EngineConfig) -> None:
    """Add a custom engine configuration"""
    _engine_configs[config.name] = config


def getEngineConfig(engine: str) -> EngineConfig:
    if engine not in _engine_configs:
        raise ValueError(
            f"Unsupported engine: {engine}. Supported: {list(_engine_configs.keys())}"
        )
    return _engine_configs[engine]


def prepare_container_cmd(
    tool_name: str,
    cmd_args: list[str],
    extra_args: list[str] | None = None,
) -> list[str]:
    """Prepare container command for research scripts"""
    config = _engine_configs[_engine]

    # For local engine, return original command
    if not config.path:
        return cmd_args

    # Get image name
    image_name = config.image_prefix + _images.get(tool_name, tool_name)

    # Build container command
    cmd = [config.path] + config.run_args

    # Add container name if engine supports it
    if config.name_flag:
        random_name = str(uuid.uuid4()).split("-", 1)[0]
        cmd.extend([config.name_flag, random_name])

    # Add extra arguments
    if extra_args and config.convert_args_func:
        extra_args = config.convert_args_func(extra_args)
        cmd.extend(extra_args)

    # Add image and command
    cmd.append(image_name)
    cmd.extend(cmd_args)
    return cmd


def runTool(
    tool_name: str,
    cmd_args: list[str],
    capture_output: bool = False,
    cwd: str | None = None,
    extra_args: list[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """
    Run a tool using the configured engine

    Args:
        tool_name: Name of the tool (from _images) or full image name
        cmd_args: Command arguments as list
        capture_output: Whether to capture output
        cwd: Working directory
        extra_args: Additional container arguments
    """
    config = _engine_configs[_engine]

    # For local engine, run directly without container
    if not config.path:
        return runShell(cmd_args, capture_output=capture_output, cwd=cwd)

    # Prepare container command
    cmd = prepare_container_cmd(tool_name, cmd_args, extra_args)
    return runShell(cmd, capture_output=capture_output, cwd=cwd)
