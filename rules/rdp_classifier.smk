rule taxonomic_assignment:
    input:
        config["dir"] + "/cat.denoised.nonchimeras"
    output:
        config["dir"] + "/rdp.out.tmp"
    threads: 4 
    params:
        memory = config["RDP"]["memory"],
        options = lambda wildcards: rdp_options(config)
    shell:
        """
        python3 python_scripts/parallel_rdp.py \
            --input {input} \
            --output {output} \
            --threads {threads} \
            --memory '{params.memory}' \
            --options '{params.options}'
        """
