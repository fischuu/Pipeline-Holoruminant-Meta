def _get_read_file(wildcards, end):
    assert end in ["forward_filename", "reverse_filename"]
    return samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ][end].values[0]


def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    return _get_read_file(wildcards, end="forward_filename")


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    return _get_read_file(wildcards, end="reverse_filename")
