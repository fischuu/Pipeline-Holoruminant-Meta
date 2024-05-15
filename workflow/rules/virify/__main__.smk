include: "virify.smk"


rule virify:
    """Run virify"""
    input:
        rules.virify__run.input,
