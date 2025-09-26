def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024


def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt

# --- Escalation functions ---

def get_escalation_order(rule_name):
    return ESCALATION_CONFIG.get(rule_name, ["small", "medium", "large"])

def get_resources_old(wildcards, attempt, rule_name):
    escalation_order = get_escalation_order(rule_name)
    profile_idx = min(attempt - 1, len(escalation_order) - 1)
    profile_name = escalation_order[profile_idx]
    resources = config["resource_sets"][profile_name].copy()

    # Ensure numeric types
    resources["runtime"] = int(resources["runtime"])
    resources["cpus"] = int(resources["cpus"])
    resources["mem_mb"] = int(resources["mem_mb"])
    resources["mem_per_cpu"] = int(resources["mem_mb"] / resources["cpus"])
    resources["nvme"] = int(resources["nvme"])
    return resources
  
def get_resources(wildcards, attempt, rule_name):
    escalation_order = get_escalation_order(rule_name)
    profile_idx = min(attempt - 1, len(escalation_order) - 1)
    profile_name = escalation_order[profile_idx]
    resources = config["resource_sets"][profile_name].copy()

    # Ensure numeric types
    resources["runtime"] = int(resources["runtime"])
    resources["cpus"] = int(resources["cpus"])
    resources["mem_mb"] = int(resources["mem_mb"])
    resources["mem_per_cpu"] = int(resources["mem_mb"] / resources["cpus"])
    resources["nvme"] = resources["nvme"]
    return resources

def esc(key, rule_name):
    """Return a lambda for Snakemake resources/threads, bound to a specific rule name."""
    return lambda wc, attempt: get_resources(wc, attempt, rule_name)[key]

def esc_val(key, rule_name, attempt=1):
    return get_resources(None, attempt, rule_name)[key]
