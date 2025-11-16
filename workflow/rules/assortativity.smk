rule done_assortativity_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/assortativity/cancer-comm-assort-enrich-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/assortativity/normal-comm-assort-enrich-{cutoff}.png",
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/assortativity/assort-{cutoff}.txt",
    shell:
        "echo done > {output}"

# Gets assortativity plots for communities in a tissue
rule get_assortativities_plot:
    input:
        chr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-chr-assortativity-{cutoff}.tsv",
        expr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-expr-assortativity-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-info-all-louvain-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv",
        enrich=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/enrichments/go-{cond}-comm-all-{cutoff}.tsv"
    output:
        comm_summary=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-comm-summary-{cutoff}.tsv",
        comm_assortativity_png=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/assortativity/{cond}-comm-assort-enrich-{cutoff}.png",
    params:
        tissue="{tissue}",
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_assortativity_plot_{cutoff}.log"
    script:
        "../scripts/assortativityEnrichmentPlot.R"

# Gets expression and chromosomal assortativity values for communities in a tissue
rule get_assortativities:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv",
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-all-louvain-{cutoff}.tsv",
        expression=config["datadir"]+"/{tissue}/deg/deg_results.tsv"
    output:
        chr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-chr-assortativity-{cutoff}.tsv",
        expr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-expr-assortativity-{cutoff}.tsv",
        diff_expr_summary=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/assortativity/{cond}-diff-expr-summary-{cutoff}.tsv",
    params:
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_assortativity_{cutoff}.log"
    script:
        "../scripts/communitiesAssortativity.R"
