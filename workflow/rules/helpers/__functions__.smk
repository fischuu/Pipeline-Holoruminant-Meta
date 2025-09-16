def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024


def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt

def get_resources(attempt=None, escalation_order=None):
    if attempt is None:
        attempt = 1  # default for dry-run / first attempt

    if escalation_order is None:
        escalation_order = ["small", "medium", "large"]

    profile_idx = min(attempt - 1, len(escalation_order) - 1)
    profile_name = escalation_order[profile_idx]

    resources = config["resource_sets"][profile_name].copy()

    # Ensure numeric types
    resources["runtime"] = int(resources["runtime"])
    resources["cpus"] = int(resources["cpus"])
    resources["mem_mb"] = int(resources["mem_mb"])
    resources["mem_per_cpu"] = int(resources["mem_mb"] / resources["cpus"])

    return resources
