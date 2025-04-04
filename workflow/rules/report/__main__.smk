include: "__functions__.smk"
include: "preprocess.smk"
include: "reads.smk"
include: "step.smk"
include: "sample.smk"

rule report:
    """Report by step and by assembly"""
    input:
        rules.report__reads.output,
        rules.report__step.input,
        rules.report__sample.input,
        rules.report__preprocess.output,

rule report_preprocess:
    """Report preprocess module"""
    input:
        rules.report__preprocess.output,

rule report_reads:
    """Report reads module"""
    input:
        rules.report__reads.output,
