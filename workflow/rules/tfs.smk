rule get_tf_pgscore_plots:
    input:
        #interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-interactions-{cutoff}.tsv",
        #vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/{cond}-vertices-{cutoff}.tsv"
        cancer_interactions=config["datadir"]+"/cancer/{tissue}_edges.csv",
        cancer_vertices=config["datadir"]+"/cancer/{tissue}_nodes.csv",
        normal_interactions=config["datadir"]+"/normal/{tissue}_edges.csv",
        normal_vertices=config["datadir"]+"/normal/{tissue}_nodes.csv"
    output:
        "figs/{tissue}_tf_pgcores.png"
    params:
        tissue="{tissue}",
        tfs=config["datadir"]+"/TFs_Ensembl_v_1.01.txt"
    script:
        "../scripts/pgScoresTF.R"
