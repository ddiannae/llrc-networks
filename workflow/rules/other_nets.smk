# Gets Gene Ontology (GO) enrichments for all communities in a tissue
rule get_go_enrichments_other:
    input:
        membership=config["datadir"]+"/{disease}/{cond}-comm-all-louvain.tsv",
        universe=config["datadir"]+"/{disease}/genes.tsv"
    output:
        config["datadir"]+"/{disease}/go-{cond}-comm-all.tsv",
    log:
        config["datadir"]+"/{disease}/log/go_enrichment_{cond}_all.log"
    script:
        "../scripts/go_enrichments.R"

# Get networks communities using an specified algorithm
rule get_comms_other:
    input:
        config["datadir"]+"/{disease}/log/done.txt",
        interactions=config["datadir"]+"/{disease}/{cond}-interactions.tsv",
        vertices=config["datadir"]+"/{disease}/{cond}-vertices.tsv"
    output:
        comm=config["datadir"]+"/{disease}/{cond}-comm-all-louvain.tsv",
        comm_info=config["datadir"]+"/{disease}/{cond}-comm-info-all-louvain.tsv",
    params:
        comm_type="all",
        commalg="louvain"
    log:
        config["datadir"]+"/{disease}/log/{cond}_communities_all_louvain.log"
    script:
        "../scripts/communities.R"

rule setup_log:
    output:
        config["datadir"]+"/{disease}/log/done.txt"
    shell:
        """
        touch {output}
        """