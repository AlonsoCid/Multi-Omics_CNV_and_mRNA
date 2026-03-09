rule all:
    input:
        "results/pca.png",
        "results/pca_exp.png",
        "results/mfa_layer.png",
        "results/mfa_circos.pdf",
        "results/DIABLO_model_perf.png",
        "results/DIABLO_plots.pdf",
        "results/DIABLO_heatmap.png"

rule eda_pca:
    input:
        data = "data/pollack.RData",
        script = "scripts/01_pca.R"
    output:
        pca_cnv = "results/pca.png",
        pca_exp = "results/pca_exp.png"
    conda:
        "environment.yml"
    shell:
        "Rscript {input.script}"

rule unsupervised_mfa:
    input:
        data = "data/pollack.RData",
        script = "scripts/02_unsupervised_mfa.R"
    output:
        mfa_layer = "results/mfa_layer.png",
        mfa_circos = "results/mfa_circos.pdf"
    conda:
        "environment.yml"
    shell:
        "Rscript {input.script}"

rule supervised_diablo:
    input:
        data = "data/pollack.RData",
        script = "scripts/03_supervised_diablo.R"
    output:
        diablo_perf = "results/DIABLO_model_perf.png",
        diablo_plots = "results/DIABLO_plots.pdf",
        diablo_heatmap = "results/DIABLO_heatmap.png"
    conda:
        "environment.yml"
    shell:
        "Rscript {input.script}"
