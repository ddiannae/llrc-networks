rule get_network_plots_output:
    input:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/communities/comm-size-histogram-network-all-louvain-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/normal-comm-network-all-louvain-{cutoff}.tsv",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/cancer-comm-network-all-louvain-{cutoff}.tsv",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/mi-density-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/degree-distribution-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/shared-vertices-{cutoff}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/network-plots-{cutoff}.txt"
    shell:
        "echo done > {output}"

rule get_cancer_normal_interaction:
    input:
        inter_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/normal-interactions-{cutoff}.tsv",
        inter_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/cancer-interactions-{cutoff}.tsv",
        ver_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/normal-vertices-{cutoff}.tsv",
        ver_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/cancer-vertices-{cutoff}.tsv",
    output:
        inter_normal_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/normal_only-interactions-{cutoff}.tsv",
        inter_cancer_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/cancer_only-interactions-{cutoff}.tsv",
        inter_shared=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/shared-interactions-{cutoff}.tsv",
        ver_normal_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/normal_only-vertices-{cutoff}.tsv",
        ver_cancer_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/cancer_only-vertices-{cutoff}.tsv",
        ver_shared=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/shared-vertices-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/networks_interaction_{cutoff}.log"
    script:
        "../scripts/networksComparison.R"

rule get_comms_plots:
    input:
        comm_info_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/normal-comm-info-{ctype}-{commalg}-{cutoff}.tsv",
        comm_info_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/cancer-comm-info-{ctype}-{commalg}-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/"+config["netdir"]+"_"+config["data_format"]+"_plots/communities/comm-{plotfactor}-{plottype}-network-{{ctype}}-{{commalg}}-{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/communities_{ctype}_{commalg}_{cutoff}_plots.log"
    script:
        "../scripts/communitiesStatsPlots.R"
   
rule get_comms_net:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv",
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-{ctype}-{commalg}-{cutoff}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-network-{ctype}-{commalg}-{cutoff}.tsv",
    params:
        tissue="{tissue}",
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_communities_network_{ctype}_{commalg}_{cutoff}.log"
    script:
        "../scripts/communitiesNetwork.R"

rule get_comms:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        comm=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-{ctype}-{commalg}-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/communities/{cond}-comm-info-{ctype}-{commalg}-{cutoff}.tsv",
    params:
        comm_type="{ctype}",
        commalg="{commalg}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_communities_{ctype}_{commalg}_{cutoff}.log"
    script:
        "../scripts/communities.R"

rule get_degree_distribution_plots:
    input:
        network_normal=config["datadir"]+"/{tissue}/rdata/normal_"+config["netdir"]+"_"+config["data_format"]+"_{cutoff}.RData",
        network_cancer=config["datadir"]+"/{tissue}/rdata/cancer_"+config["netdir"]+"_"+config["data_format"]+"_{cutoff}.RData",
    params:
        tissue="{tissue}"
    output:
        dd=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/degree-distribution-{cutoff}.png",
        cdd=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/degree-distribution-cumulative-{cutoff}.png",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/degree_distribution_{cutoff}_plots.log"
    script:
        "../scripts/degreeDistributionPlots.R"

rule get_network_stats:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        network=config["datadir"]+"/{tissue}/rdata/{cond}_"+config["netdir"]+"_"+config["data_format"]+"_{cutoff}.RData",
        network_stats=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-network-stats-{cutoff}.tsv",
        node_attributes=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-node-attributes-{cutoff}.tsv",
        edge_attributes=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-edge-attributes-{cutoff}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_network_stats_{cutoff}.log"
    script:
        "../scripts/networkStats.R"

rule get_distribution_plots:
    input:
        normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/normal-interactions-{cutoff}.tsv",
        cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/cancer-interactions-{cutoff}.tsv"
    output:
        boxplot=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/mi-boxplot-network-{cutoff}.png",
        density=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"_plots/mi-density-network-{cutoff}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/network_{cutoff}_distribution_plots.log"
    script:
        "../scripts/MIDistributionPlots.R"

rule get_network_tables:
    input:
        mi_matrix=get_mi_matrix,
        done=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/done.txt"
    output:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv"
    params:
        annot_cytobands=config["biomart"],
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cutoff="{cutoff}",
        cond="{cond}"
    threads: 18 
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/{cond}_network_table_{cutoff}.log"
    script:
        "../scripts/networkTables.R"
    
