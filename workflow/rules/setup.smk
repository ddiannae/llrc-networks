rule setup_distance_log:
    output:
        config["datadir"]+"/{tissue}/"+get_dist_dir()+"/log/done.txt"
    shell:
        """
        touch {output}
        """

rule setup_network_log:
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["data_format"]+"/log/done.txt"
    shell:
        """
        touch {output}
        """
