# fastp
def _get_adapter(wildcards, end):
    """Get the adapter of the en from a file"""
    assert end in ["forward", "reverse"]
    end = "forward_adapter" if end == "forward" else "reverse_adapter"
    adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][end].tolist()[0]
    return adapter


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return _get_adapter(wildcards, end="forward")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return _get_adapter(wildcards, end="reverse")


# bowtie2
def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    lb_field = f"LB:{sample_id}_{library_id}"
    pl_field = "PL:Illumina"
    sm_field = f"SM:{sample_id}"
    return f"{lb_field}\t{pl_field}\t{sm_field}"


def _get_input_file_for_host_mapping(wildcards, end):
    """Compose the input file for host mapping"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample_id}.{wildcards.library_id}_{end}.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        PRE_BOWTIE2
        / f"non{prev_genome}"
        / f"{wildcards.sample_id}.{wildcards.library_id}_{end}.fq.gz"
    )


def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return _get_input_file_for_host_mapping(wildcards, end="forward")


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
    return _get_input_file_for_host_mapping(wildcards, end="reverse")


# finals
def _get_final_file_from_pre(wildcards, end):
    """Get the last host or FASTP forward or reverse"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    last_host = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
    return PRE_BOWTIE2 / f"non{last_host}" / f"{sample_id}.{library_id}_{end}.fq.gz"


def get_final_forward_from_pre(wildcards):
    """Get the last host forward file or the result from FASTP"""
    return _get_final_file_from_pre(wildcards, end="forward")


def get_final_reverse_from_pre(wildcards):
    """Get the last host reverse file or the result from FASTP"""
    return _get_final_file_from_pre(wildcards, end="reverse")
