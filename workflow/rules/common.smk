import glob

def get_output_files(wildcards):
    files = []
    if config["end"] == "other":
                files.append(f'{config["datadir"]}/diabetes/go-normal-comm-all.tsv')
                files.append(f'{config["datadir"]}/diabetes/go-cases-comm-all.tsv')
                files.append(f'{config["datadir"]}/alzheimer/go-normal-comm-all.tsv')
                files.append(f'{config["datadir"]}/alzheimer/go-cases-comm-all.tsv')
    else:
        for t in config["tissues"]:
            if config["end"] == "plots":
                files.append(f'{config["datadir"]}/{t}/{config["netdir"]}_{config["data_format"]}/enrichments/all-enrichments-{config["cutoff"]}.txt')
                files.append(f'{config["datadir"]}/{t}/{config["netdir"]}_{config["data_format"]}_plots/assortativity/assort-{config["cutoff"]}.txt')
                files.append(f'{config["datadir"]}/{t}/{config["netdir"]}_{config["data_format"]}_plots/network-plots-{config["cutoff"]}.txt')
                files.append(f'{config["datadir"]}/{t}/{get_fig_dir()}/intra-plots.txt')
                files.append(f'{config["datadir"]}/{t}/{get_fig_dir()}/intra-inter-plots.txt')
            elif config["end"] == "boots":
                if t == "esophagus":
                    files.append(f'{config["datadir"]}/{t}/correlation/bootstrap_spearman_50/normal-log-bins-summ.tsv')
                    files.append(f'{config["datadir"]}/{t}/correlation/bootstrap_spearman_100/normal-log-bins-summ.tsv')
                else:
                    files.append(f'{config["datadir"]}/{t}/correlation/bootstrap_spearman_50/cancer-log-bins-summ.tsv')
                    files.append(f'{config["datadir"]}/{t}/correlation/bootstrap_spearman_100/cancer-log-bins-summ.tsv')
            elif config["end"] == "tfs":
                files.append(f'figs/{t}_tf_pgcores.png')

    return files

def get_dist_dir():
    return f'{config["distdir"]}_{config["data_format"]}'

def get_fig_dir():
    return f'{config["figdir"]}_{config["data_format"]}'

def get_mi_matrix(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}/correlation/{config["data_format"]}_ensembl_{wildcards.cond}.adj'

def get_mi_matrix_boot(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}/correlation/bootstrap_spearman_{wildcards.s}/spearman_bootstrap_{wildcards.n}.tsv'

def getGeneUniverse(wildcards):
    if config["data_format"] == "tpm" or config["data_format"] == "arsyn_tpm":
        return f'{config["datadir"]}/{wildcards.tissue}/results/tpm_ensembl_genes.tsv'
    elif config["data_format"] == "deseq2" or config["data_format"] == "arsyn_deseq2":
        return f'{config["datadir"]}/{wildcards.tissue}/results/deseq2_ensembl_genes.tsv'

def getDEGFile(wildcards):
        return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] +
        "/deg/*_*_*_si-arsyn_deg_results.tsv")]
