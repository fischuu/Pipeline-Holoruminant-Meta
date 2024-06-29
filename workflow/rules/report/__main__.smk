include: "__functions__.smk"
include: "preprocess.smk"
include: "step.smk"
include: "sample.smk"

rule report:
    """Report by step and by assembly"""
    input:
        rules.report__step.input,
        rules.report__sample.input,

rule preprocess:
    """Report the results from the preprocess module"""
    input:
        rules._report__preprocess_rules.output,