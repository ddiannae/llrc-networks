
rule get_intra_inter_plots:
    input:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-count-onek-bins.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-count-log-bins.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-count-onek-chunks.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-count-log-chunks.png",
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/intra-inter-ks-log-bins.tsv",
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/intra-inter-ks-onek-bins.tsv",
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-plots.txt"
    shell:
        "echo done > {output}"

# Gets Kolmogorov-Smirnov test results for intra-chromosomal interactions
rule get_intra_inter_ks:
    input:
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-intra-inter-count-{bintype}-{bc}.tsv",
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-intra-inter-count-{bintype}-{bc}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/intra-inter-ks-{bintype}-{bc}.tsv"
    threads: 38
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/intra_inter_ks_{bintype}_{bc}.log"
    script:
        "../scripts/intraInterKS.R"

# Gets plots for intra- and inter- chromosomal interactions
rule get_intra_plot:
    input:
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-intra-inter-count-{bintype}-{bc}.tsv",
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-intra-inter-count-{bintype}-{bc}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-inter-count-{bintype}-{bc}.png"
    params:
        bintype="{bintype}",
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/intra_inter_count_{bintype}_{bc}_plot.log"
    script:
        "../scripts/intraInterPlot.R"

# Gets counts for intra- and inter- chromosomal interactions
rule get_intra_inter_count:
    input:
        mi_matrix=get_mi_matrix
    output:
        onek_chunks=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-intra-inter-count-onek-chunks.tsv",
        onek_bins=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-intra-inter-count-onek-bins.tsv",
        log_bins=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-intra-inter-count-log-bins.tsv",
        log_chunks=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-intra-inter-count-log-chunks.tsv"
    params:
        cond="{cond}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 18
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_get_intra_inter_count.log"
    script:
        "../scripts/intraInterCount.R"

rule get_intra_inter_summ_cancer:
    input:
        expand(config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-intra-inter-count-log-bins-{n}.tsv", n = [x+1 for x in range(int(config["bsamples"]))], allow_missing=True)
    output:
        config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-log-bins-summ.tsv"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_bootstrap_log_bins_summ_{s}.log"
    script:
        "../scripts/bootstrapIntraInter.R"

rule get_intra_inter_summ_cancer_aracne:
    input:
        expand(config["datadir"]+"/{tissue}/correlation/bootstrap_samples/{cond}-intra-inter-count-log-bins-{n}.tsv", n = [x+1 for x in range(int(config["bsamples"]))], allow_missing=True)
    output:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/{cond}-log-bins-summ.tsv"
    log:
        config["datadir"]+"/{tissue}/log/{cond}_bootstrap_log_bins_summ.log"
    script:
        "../scripts/bootstrapIntraInter.R"

rule get_intra_inter_count_bootstrp:
    input:
        mi_matrix=get_mi_matrix_boot
    output:
        onek_chunks=config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-intra-inter-count-onek-chunks-{n}.tsv",
        onek_bins=config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-intra-inter-count-onek-bins-{n}.tsv",
        log_bins=config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-intra-inter-count-log-bins-{n}.tsv",
        log_chunks=config["datadir"]+"/{tissue}/correlation/bootstrap_spearman_{s}/{cond}-intra-inter-count-log-chunks-{n}.tsv"
    params:
        cond="{cond}",
        annot=config["datadir"]+"/{tissue}/annot.RData",
        del_input=1
    threads: 5
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_get_intra_inter_count_{n}_{s}.log"
    script:
        "../scripts/intraInterCount.R"
