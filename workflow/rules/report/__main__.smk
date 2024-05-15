include: "__functions__.smk"
include: "step.smk"
include: "sample.smk"


rule report:
    """Report by step and by assembly"""
    input:
        rules.report__step.input,
        rules.report__sample.input,
