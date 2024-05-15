include: "hosts.smk"


rule reference:
    input:
        rules.reference__hosts.input,
