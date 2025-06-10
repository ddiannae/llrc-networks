# Loss of long-range co-expression is a common trait in cancer

## Co-expression Networks and Distance Analysis Workflow

This repository contains a Snakemake-based workflow for constructing and analysing gene co-expression networks from mutual information (MI) matrices. The workflow is designed to compare cancer and normal tissue networks, identify network communities, perform functional enrichment, and generate a variety of summary statistics and plots.

The input for this workflow comes from the results from [RNA-Seq pipeline for TCGA and USCS Xena datasets](https://github.com/ddiannae/llrc-pipeline)

## Overview

The workflow takes as input MI matrices (gene x gene), gene annotation files, and expression data, and produces:

- Annotated interaction and vertex tables for co-expression networks
- Metrics for gene-pair distance vs co-expression for intra-chromosomal interactions
- Network statistics and community detection
- Functional enrichment analyses (GO, KEGG, oncogenic signatures, NCG) in networks
- Assortativity and other network metrics
- Plots and summary tables

## Main Steps

1. **Network Table Generation**
   - Annotate MI matrices and extract top interactions (`networkTables.R`)
   - Filter intra-chromosomal interactions and calculate distances (`intraInteractions.R`)

2. **Network Statistics**
   - Calculate network, node, and edge attributes (`networkStats.R`)
   - Compare cancer and normal networks (`networksComparison.R`)

3. **Community Detection and Analysis**
   - Detect network communities using various algorithms (`communities.R`)
   - Build community-level networks and calculate assortativity (`communitiesNetwork.R`, `communitiesAssortativity.R`)
   - Summarize and plot community statistics (`communitiesStatsPlots.R`, `assortativityEnrichmentPlot.R`)

4. **Functional Enrichment**
   - Perform GO, KEGG, oncogenic signature, and NCG enrichment for communities (`go_enrichments.R`, `other_enrichments.R`)

5. **Distribution and Distance Analysis**
   - Plot MI and degree distributions (`MIDistributionPlots.R`, `degreeDistributionPlots.R`)
   - Bin intra-chromosomal interactions and calculate statistics (`binStats.R`, `binFittingByChr.R`, `binDistanceByChrPlot.R`)
   - Perform statistical tests between bins (`binTest.R`, `heatmapTtest.R`)

6. **Intra/Inter Chromosomal Analysis**
   - Calculate and plot intra/inter chromosomal interaction fractions (`intraInterCount.R`, `intraInterPlot.R`, `intraInterKS.R`)


## Repository Structure

- `workflow` — Snakemake rules and R scripts used in the workflow.
- `config.yaml` — Example configuration file for Snakemake pipeline.
- `Snakefile` — Main workflow definition.

