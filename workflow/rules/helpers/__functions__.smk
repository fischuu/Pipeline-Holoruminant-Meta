def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024


def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt

def get_resources(wildcards, starting_profile="small", escalation_order=None):
    """
    Return the resource set for the current attempt, with escalation.

    Args:
        wildcards: Snakemake wildcards (needed for lambda in resources)
        starting_profile: the resource profile to use for the first attempt
        escalation_order: list of profile names in escalation order; 
                          if None, defaults to ["small", "medium", "highmem", "longrun", "hm_longrun"]

    Returns:
        dict with resource allocations from config["resource_sets"]
    """
    attempt = snakemake.attempt
    if escalation_order is None:
        escalation_order = ["small", "medium", "highmem", "longrun", "hm_longrun"]

    try:
        start_idx = escalation_order.index(starting_profile)
    except ValueError:
        start_idx = 0  # fallback to first in escalation_order

    # pick the profile for this attempt
    profile_idx = min(start_idx + attempt - 1, len(escalation_order) - 1)
    profile_name = escalation_order[profile_idx]

    return config["resource_sets"][profile_name]

