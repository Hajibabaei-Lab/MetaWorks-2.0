---
title: 'MetaWorks 2.0: A Modular and Containerized Workflow Platform for Multi-Marker Metabarcoding'
tags:
  - Python
  - Snakemake
  - bioinformatics
  - metabarcoding
  - amplicon sequencing
  - environmental DNA
  - exact sequence variants
  - operational taxonomic units
authors:
  - name: Yaye Song
    orcid: 0009-0009-3579-8976
    affiliation: 1
  - name: Mehrdad Hajibabaei
    orcid: 0000-0002-8859-7977
    affiliation: 1
    corresponding: true
affiliations:
  - name: Centre for Biodiversity Genomics, Biodiversity Institute of Ontario, University of Guelph, Guelph, Ontario, Canada
    index: 1
date: 31 August 2026
bibliography: paper.bib
---

# Summary

MetaWorks is an open-source workflow platform for processing short-read amplicon sequencing data, with the current release primarily validated using Illumina datasets. It converts raw reads into Exact Sequence Variants (ESVs, also called Amplicon Sequence Variants or ASVs) and Operational Taxonomic Units (OTUs), and supports metabarcoding markers commonly used in biodiversity and ecological studies, including COI, 12S, 16S, 18S, ITS, 28S, and rbcL. A run can include read pairing, adapter and primer trimming, denoising, chimera removal, optional ITS-region extraction and pseudogene filtering, taxonomic classification, and results assembly.

MetaWorks 2.0 reorganizes the six workflow-specific Snakefiles used in MetaWorks v1.13 into a modular architecture built from reusable functional rule components. It also adds layered configuration with Pydantic validation, Docker and Apptainer/Singularity execution backends, a REST API, a Vue 3 web interface, and automated testing. The command-line workflow, API, and web interface use the same underlying configuration and pipeline implementation, allowing analyses to move between personal workstations and shared Linux servers without maintaining separate workflow code.

# Statement of need

Amplicon sequencing and metabarcoding target short, taxonomically informative genomic regions and underpin biodiversity monitoring, environmental DNA (eDNA) surveys, microbiome research, and biosurveillance [@taberlet2018environmental; @porter2018scalingup]. These applications increasingly combine markers that represent different taxonomic groups and different sequence types. Ribosomal genes, internal transcribed spacers, and protein-coding mitochondrial or chloroplast markers may require different preprocessing, filtering, and taxonomic-classification steps.

MetaWorks is intended for researchers in ecology, conservation biology, agriculture, environmental assessment, and public health who need to process multi-marker amplicon data in a consistent and reproducible way. MetaWorks v1 [@porter2022metaworks] introduced a unified Snakemake-based environment for these analyses. By v1.13.0 [@metaworks_v113], six workflow-specific Snakefiles were used across analysis modes and read layouts, with repeated processing logic across the paired-end, single-read, and dual-indexed ESV workflows. The v1 series also relied on command-line execution in a Conda environment, hand-written configuration files without schema validation, and limited automated testing. Version 2.0 retains the established analytical approach while redesigning this workflow organization for maintainability, reproducible deployment, and easier use, and corrects script-logic issues affecting HMM-based pseudogene filtering identified during the v2 migration.

# State of the field

Several established tools serve amplicon and metabarcoding research. QIIME 2 [@bolyen2019qiime2] and mothur [@schloss2009mothur] provide extensive ecosystems for microbial community analysis, DADA2 [@callahan2016dada2] performs high-resolution ASV inference, OBITools [@boyer2016obitools] supports eDNA and metabarcoding sequence processing, and nf-core/ampliseq [@straub2020nfcore] provides a portable Nextflow-based workflow. The broader landscape includes many additional pipelines with different assumptions and target applications [@hakimzadeh2024pipelines].

MetaWorks addresses a complementary need: applying a shared processing framework across heterogeneous marker classes used in biodiversity studies. Mixed-marker projects can require marker-specific operations such as ITS extraction or translation-aware pseudogene filtering for protein-coding loci, together with different reference classifiers and genetic codes. Implementing these requirements as extensions to several otherwise capable analysis ecosystems would still leave the research program maintaining separate marker-specific configurations and integration logic. MetaWorks therefore provides one configurable workflow in which these differences are represented as reusable modules and marker profiles. This design also preserves continuity with the analytical approach, classifiers, and research workflows developed around MetaWorks v1.

# Software design

## Modular workflow architecture

MetaWorks 2.0 is built on Snakemake [@moelderSnakemake2021] with a modular rule system. In MetaWorks v1.13.0 [@metaworks_v113], separate Snakefiles defined different analysis modes and read layouts. In particular, the paired-end, single-read, and dual-indexed ESV workflows independently implemented their read QC and trimming front ends and repeated dereplication, denoising, chimera-removal, and compression stages. OTU and global workflows instead consumed outputs generated by the ESV workflows.

Version 2.0 replaces this workflow-specific organization with a common workflow entry point and reusable function-specific `.smk` components. A central registry resolves module dependencies and determines the modules required by the selected workflow and configuration. These modules cover input preparation, read pairing with SeqPrep [@stjohn2016seqprep], adapter and primer trimming with Cutadapt [@martin2011cutadapt], denoising and chimera removal with VSEARCH [@rognes2016vsearch], OTU clustering, ITS extraction with ITSx [@bengtsson-palme2013itsx], pseudogene filtering [@porter2021pseudogene], taxonomic classification, statistics, global processing, and results assembly.

The architectural transition from workflow-specific Snakefiles in MetaWorks v1.13 to reusable, configuration-selected modules in version 2.0 is summarized in Figure 1.

![Architectural comparison of MetaWorks v1.13 and MetaWorks 2.0. MetaWorks v1.13 used six workflow-specific Snakefiles. Paired, single-read, and dual-indexed ESV workflows independently implemented read processing and repeated denoising stages, whereas OTU and global workflows consumed denoised ESV output. MetaWorks 2.0 uses a common workflow entry point, centralized module and dependency resolution, and reusable rule modules selected according to workflow and marker configuration.[]{label="fig:architecture"}](architecture_v1_v2.png)

Detailed diagrams of individual ESV, OTU, ITS, input, global, and classification workflows are available in the repository's [workflow diagram documentation](https://github.com/Hajibabaei-Lab/MetaWorks-2.0/tree/main/docs/diagrams).

## Configuration and classification

Configuration is assembled from three layers: pipeline defaults, a marker-specific profile, and sparse user overrides. Fourteen profiles cover COI, 12S, 16S, 18S, ITS, 28S, rbcL, and related variants. Profiles that have been exercised with representative data provide ready-to-use settings, while additional starter profiles establish the expected marker structure for further expert refinement. Pydantic [@pydantic] validates the merged configuration before workflow execution, allowing configuration errors to be reported before computational jobs are launched.

MetaWorks supports the RDP Naive Bayesian Classifier [@wang2007rdp] and SINTAX [@edgar2016sintax]. Marker-specific reference data can be supplied by users, and curated RDP-trained reference sets for COI [@porter2018co1; @porter_co1_rdp] and other markers used by our research group are maintained separately. Keeping reference assets outside the core software allows classifier databases to evolve independently of the pipeline release.

## Interfaces and deployment

The Snakemake workflow remains directly usable from the command line. A FastAPI [@fastapi] backend and Vue 3 interface provide optional configuration, submission, monitoring, and results-management capabilities using the same underlying pipeline. Schema-driven forms reduce duplication between frontend configuration options and backend validation. Neither the API nor the web interface provides built-in authentication; shared or network-accessible deployments should therefore be placed behind an authenticating reverse proxy or equivalent access-control layer.

MetaWorks provides Conda, Docker, and Apptainer/Singularity execution backends behind a common runtime interface. The backends handle environment selection, bind mounts, and path translation so that the workflow configuration remains portable, while Docker Compose configurations support packaged web deployment. These execution modes have been exercised on Linux workstations and on our laboratory's shared Linux server. Scheduler-based HPC integration (e.g., SLURM) has not yet been validated.

# Testing and reproducible example data

MetaWorks 2.0 includes an automated backend test suite covering the REST API surface, configuration merging and validation, module enabling and dependency resolution, target selection, JSON schema generation, and the command-line interfaces of pipeline helper scripts. Frontend checks include TypeScript type checking, unit tests, and a production build. These checks run through GitHub Actions on pushes and pull requests. Automated coverage is concentrated on unit, configuration, command-line, and API behavior; broader automated end-to-end testing with full sequencing datasets and expanded web-interface coverage remain areas for continued development.

In addition to automated software tests, the repository includes a reduced COI sequencing dataset for functional workflow testing. The files are derived from sequencing data associated with the Biomonitoring 2.0 Refined study [@riley2025biomonitoring], whose complete sequence data are deposited in the NCBI Sequence Read Archive under BioProject PRJNA1201794 [@riley2025bioproject]. The reduced dataset provides a practical end-to-end example that can be distributed with the repository, unlike the substantially larger datasets used during development and testing of MetaWorks v1.

# Research impact statement

MetaWorks 2.0 is a substantial software redesign of an established research workflow rather than a new one-off analysis pipeline. Published studies explicitly report using MetaWorks across several releases and research contexts: v1.0.0 for multi-marker urban-harbour biodiversity assessment [@robinson2022urbanharbour], v1.4.0 for multi-marker forest-soil monitoring across Canada [@smenderovac2022woodash; @smenderovac2023harvesting], and v1.13.0 for freshwater-fish eDNA metabarcoding using COI and 12S markers [@morey2024blindspots] and community-based biodiversity monitoring and metaphylogeography [@riley2025biomonitoring]. These applications demonstrate use across different ecosystems, marker combinations, and research questions.

Version 2.0 responds to experience gained through this research use by replacing six workflow-specific v1 Snakemake implementations with a composable rule-module architecture, validating configuration before execution, adding containerized runtimes and optional programmatic and graphical interfaces, and expanding automated testing to 166 tests. The published use of MetaWorks demonstrates an established research need and motivates a more reproducible and maintainable successor.

# Availability and use

MetaWorks 2.0 is released under the GNU General Public License version 3 (GPLv3) and is available at <https://github.com/Hajibabaei-Lab/MetaWorks-2.0>.

The exact MetaWorks 2.0 v2.0.0 software release corresponding to this paper is archived on Zenodo [@metaworks_v2].

# Acknowledgements

This work was supported in part by a Natural Sciences and Engineering Research Council of Canada (NSERC) Undergraduate Student Research Award (USRA).

# AI usage disclosure

Generative AI tools were used during the development of MetaWorks 2.0 and in the preparation of this paper. Claude Opus 4.x (Anthropic) and ChatGPT (OpenAI) assisted with software development, including generating and refactoring Python helper scripts, Snakemake rules, and API code, and with drafting and editing this manuscript.

All AI-assisted outputs were reviewed, edited, and validated by the human authors, who made the core scientific and architectural decisions. Generated code was evaluated through the project's automated test suite and manual review. The authors take responsibility for the accuracy, originality, licensing, and legal and ethical compliance of the submitted software and manuscript.

# References
