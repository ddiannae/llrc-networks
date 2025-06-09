rule get_all_enrichments:
    input:
        expand(config["datadir"]+"/{{tissue}}/"+config["netdir"]+"_"+config["data_format"]+"/enrichments/{enrch}-{cond}-comm-all-{{cutoff}}.tsv",
        enrch=["kegg","go","ncg","onco"], cond=["cancer","normal"])
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/enrichments/all-enrichments-{cutoff}.txt"
    shell:
        "echo done > {output}"

# Gets KEGG and oncogenic gene sets (C6) enrichments for all communities in a tissue
rule get_other_enrichments:
    input:
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-{ctype}-louvain-{cutoff}.tsv",
        universe=getGeneUniverse
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/enrichments/{etype}-{cond}-comm-{ctype}-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{etype}_enrichment_{cond}_{ctype}_{cutoff}.log"
    params:
        enrich_type="{etype}"
    script:
        "../scripts/other_enrichments.R"

# Gets Gene Ontology (GO) enrichments for all communities in a tissue
rule get_go_enrichments:
    input:
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-{ctype}-louvain-{cutoff}.tsv",
        universe=getGeneUniverse
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/enrichments/go-{cond}-comm-{ctype}-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/go_enrichment_{cond}_{ctype}_{cutoff}.log"
    script:
        "../scripts/go_enrichments.R"
