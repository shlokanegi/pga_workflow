import os
import math

# For each Slurm partition name, what is its max wall time in minutes?
# TODO: Put this in the config
SLURM_PARTITIONS = [
    ("short", 60),
    ("medium", 12 * 60),
    ("long", 7 * 24 * 60)
]

#Different phoenix nodes seem to run at different speeds, so we can specify which node to run
#This gets added as a slurm_extra for the selected rules
SLURM_EXTRA = config.get("slurm_extra", None) or ""

# If set to True, jobs where we care about speed will demand entire nodes.
# If False, they will just use one thread per core.
EXCLUSIVE_TIMING = config.get("exclusive_timing", True)

def auto_slurm_extra(wildcards):
    """
    Determine Slurm extra arguments for a timed job.
    """
    if exclusive_timing(wildcards):
        return "--exclusive " + SLURM_EXTRA
    else:
        return SLURM_EXTRA

def auto_full_cluster_nodes(wildcards):
    """
    Determine number of full cluster nodes for a timed job.

    TODO: Is this really used by Slurm?
    """
    if exclusive_timing(wildcards):
        return 1
    else:
        return 0

def choose_partition(minutes):
    """
    Get a Slurm partition that can fit a job running for the given number of
    minutes, or raise an error.
    """
    for name, limit in SLURM_PARTITIONS:
        if minutes <= limit:
            return name
    raise ValueError(f"No Slurm partition accepts jobs that run for {minutes} minutes")



def get_mem_mb(input, buffer_factor=1.2, min_mb=500):
    """
    Calculate memory requirement in MB based on input file sizes.
    `input` is the snakemake input object.
    `buffer_factor` is the factor by which to multiply the total size of the input files. By default, adds a 20% buffer.
    `min_mb` is the minimum memory requirement in MB.
    """
    total_size_mb = 0
    
    def iter_paths(obj):
        if obj is None:
            return
        if isinstance(obj, (list, tuple, set)):
            for p in obj:
                yield from iter_paths(p)
        elif isinstance(obj, dict):
            for p in obj.values():
                yield from iter_paths(p)
        else:
            yield str(obj)

    for path in iter_paths(input):
        try:
            if os.path.isfile(path):
                total_size_mb += os.path.getsize(path) / 1_000_000.0
        except OSError:
            pass

    required_mem = math.ceil(total_size_mb * buffer_factor)
    return max(min_mb, int(required_mem))



def get_runtime_min(benchmark_path, default_min, buffer_factor=1.2):
    """
    NOTE: This function is not used in the workflow yet!
    Estimate runtime in minutes from a benchmark file.
    Falls back to a default if the benchmark file does not exist.
    """
    try:
        with open(benchmark_path, 'r') as f:
            # Skip header
            next(f)
            # Read first data line
            parts = next(f).strip().split('\t')
            wall_time_seconds = float(parts[0])
        
        # Apply buffer and convert to minutes
        runtime_min = math.ceil((wall_time_seconds / 60) * buffer_factor)
        return int(runtime_min)
    except (IOError, StopIteration):
        # File not found or empty, return default
        return default_min
