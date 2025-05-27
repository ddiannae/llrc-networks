rule get_distance_plots:
    input:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-distance-"+str(config["distbin"])+"-mean_fitted.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-distance-"+str(config["distbin"])+"-mean.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-size-"+str(config["sizebin"])+"-mean.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-size-"+str(config["sizebin"])+"-median.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-size-"+str(config["sizebin"])+"-mean_fitted.png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-size-bychr-"+str(config["sizebin"])+".png",
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-distance-bychr-"+str(config["distbin"])+".png",
        # config["datadir"]+"/{tissue}/"+get_fig_dir()+"/heatmap-bins-size-all-ttests-"+str(config["sizebin"])+".png",
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/intra-plots.txt"
    shell:
        "echo done > {output}"

rule get_bin_distance_plots:
    input:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-{bintype}-{binsize}-{stat}.png"
    params:
        tissue = "{tissue}",
        stat = "{stat}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{bintype}_plot_{binsize}_{stat}.log" 
    script:
        "../scripts/binDistancePlot.R"

rule get_bin_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-bins-{bintype}-all-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/fitted_bins_{bintype}_all_{binsize}.log" 
    script:
        "../scripts/binFitting.R"

rule get_bin_chr_plots:
    input:
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-bins-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-bins-{bintype}-bychr-{binsize}.tsv",
        fitted=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/fitted-bins-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/bin-{bintype}-bychr-{binsize}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{bintype}_plot_bychr_{binsize}.log" 
    script:
        "../scripts/binDistanceByChrPlot.R"

rule get_bin_chr_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-bins-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-bins-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/fitted-bins-{bintype}-bychr-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/fitted_bins_{bintype}_bychr_{binsize}.log" 
    script:
        "../scripts/binFittingByChr.R"

rule get_ttest_heatmap:
    input:
        cancer=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/cancer-bins-{bintype}-all-tests-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/normal-bins-{bintype}-all-tests-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_fig_dir()+"/heatmap-bins-{bintype}-all-ttests-{binsize}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/heatmap_tests_{bintype}_{binsize}.log" 
    threads: 18
    script:
      "../scripts/heatmapTtest.R"

rule get_bins_tests:
    input:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-all-distance-mi.tsv"
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-bins-{bintype}-all-tests-{binsize}.tsv"
    params:
        binsize="{binsize}",
        cond="{cond}",
        bintype="{bintype}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_bin_tests_{bintype}_{binsize}.log" 
    threads: 18
    script:
      "../scripts/binTest.R"

rule get_bins:
    input:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-all-distance-mi.tsv"
    output:
        by_chr=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-bins-{bintype}-bychr-{binsize}.tsv",
        all=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-bins-{bintype}-all-{binsize}.tsv"
    params:
        binsize="{binsize}",
        cond="{cond}",
        bintype="{bintype}"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_bins_{bintype}_{binsize}.log" 
    threads: 18
    script:
      "../scripts/binStats.R"

rule get_intra_interactions:
    input:
        mi_matrix=get_mi_matrix,
        log_file=config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/done.txt"
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/{cond}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/{cond}_intra_interactions.log" 
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 18
    script:
        "../scripts/intraInteractions.R"

