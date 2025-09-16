def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024


def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt

def get_resources(wildcards, attempt, escalation_order=None):
    """
    Return the resource set for the current attempt, with escalation.

    Args:
        wildcards: Snakemake wildcards (needed for lambda in resources)
        attempt: Snakemake attempt counter (starts at 1)
        escalation_order: list of profile names in escalation order.
                          Example: ["small", "medium", "large"]

    Returns:
        dict with resource allocations from config["resource_sets"]
    """
    if escalation_order is None:
        escalation_order = ["small", "medium", "highmem", "longrun", "hm_longrun"]

    # clamp attempt to list length
    profile_idx = min(attempt - 1, len(escalation_order) - 1)
    profile_name = escalation_order[profile_idx]

    return config["resource_sets"][profile_name]
