###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in qc.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# SAMPLES          - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

rule samtools_flagstat:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        flagstat = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam.flagstat',
    log:
        samtools_flagstat = f'{output_directory}/logs/samtools_flagstat/samtools_flagstat.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/samtools_flagstat/{{sample}}.txt'
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['samtools'],
    shell:
        """
        # Get some initial stats
        samtools flagstat {input.bam} 1> {output.flagstat} 2> {log.samtools_flagstat}
        """

rule biscuit_qc:
    input:
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
    output:
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_cv_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_mapq_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_strand_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_dup_report.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_totalBaseConversionRate.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_totalReadConversionRate.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_CpHRetentionByReadPos.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_CpGRetentionByReadPos.txt',
    params:
        assets = config['ref']['assets'],
        output_dir = f'{output_directory}/analysis/BISCUITqc',
    log:
        f'{output_directory}/logs/biscuit_qc/{{sample}}_QC.log'
    benchmark:
        f'{output_directory}/benchmarks/biscuit_qc/{{sample}}.txt'
    threads: 8
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['bedtools'],
        config['envmodules']['htslib'],
        config['envmodules']['samtools'],
        config['envmodules']['parallel'],
    shell:
        """
        set +o pipefail;
        QC.sh \
            -o {params.output_dir} \
            --vcf {input.vcf} \
            {params.assets} \
            {input.ref} \
            {wildcards.sample} \
            {input.bam} \
            2> {log}
        """

rule preseq:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        tsv = f'{output_directory}/analysis/preseq/{{sample}}.ccurve.txt',
    params:
        dir = f'{output_directory}/analysis/preseq',
        out = f'{output_directory}/analysis/preseq/{{sample}}.ccurve.txt',
        opt = config['preseq']['args_list'],
    log:
        f'{output_directory}/logs/preseq/{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/preseq/{{sample}}.txt',
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/preseq.yaml'
    envmodules:
        config['envmodules']['preseq'],
    shell:
        """
        mkdir -p {params.dir}
        preseq c_curve {params.opt} -o {params.out} -P -B {input.bam} 2> {log}
        """

def get_multiQC_params(wildcards):
    raw = config['fastqs']
    out = output_directory
    indirs = f'{raw} {out}/analysis/BISCUITqc {out}/analysis/raw_fastqc {out}/analysis/align'
    if config['fastq_screen']['run']:
        indirs += f' {out}/analysis/fastq_screen' # space needed at beginning to separate directories
    if config['trim_galore']['trim_before_biscuit']:
        indirs += f' {out}/analysis/trim_reads' # space needed at beginning to separate directories
    if config['preseq']['run']:
        indirs += f' {out}/analysis/preseq' # space needed at beginning to separate directories
    return indirs

rule multiQC:
    input:
        # raw fastqc
        get_qc_file(),
        #expand(f'{output_directory}/analysis/raw_fastqc/{{samples.sample}}-1-R{{read}}_fastqc.zip', read=[1,2], samples=SAMPLES.itertuples()),
        # flagstat
        expand(f'{output_directory}/analysis/align/{{samples.sample}}.sorted.markdup.bam.flagstat', samples=SAMPLES.itertuples()),
        # biscuit_qc
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_cv_table.txt', samples=SAMPLES.itertuples()),
        # fastq_screen
        expand(f'{output_directory}/analysis/fastq_screen/{{samples.sample}}-1-R{{read}}_screen.txt', read=[1,2], samples=SAMPLES.itertuples()) if config['fastq_screen']['run'] else [],
        # trim_galore
        get_trim_file() if config['trim_galore']['trim_before_biscuit'] else [],
        #expand(f'{output_directory}/analysis/trim_reads/{{samples.sample}}_val_1_merged.fq.gz', samples=SAMPLES.itertuples()) if config['trim_galore']['trim_before_biscuit'] else [],
        #expand(f'{output_directory}/analysis/trim_reads/{{samples.sample}}-R{{read}}_val_{{read}}_merged.fq.gz', read=[1,2], samples=SAMPLES.itertuples()) if config['trim_galore']['trim_before_biscuit'] else [],

        # preseq
        expand(f'{output_directory}/analysis/preseq/{{samples.sample}}.ccurve.txt', samples=SAMPLES.itertuples()) if config['preseq']['run'] else [],
    output:
        directory(f'{output_directory}/analysis/multiqc/multiqc_report_data',),
        f'{output_directory}/analysis/multiqc/multiqc_report.html',
    params:
        mqc_dirs = get_multiQC_params,
        output_dir = f'{output_directory}/analysis/multiqc',
    log:
        f'{output_directory}/logs/multiqc.log'
    benchmark:
        f'{output_directory}/benchmarks/multiQC.txt'
    threads: 1
    resources:
        mem_gb=config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['medium']
    conda:
        '../envs/python_packages.yaml'
    envmodules:
        config['envmodules']['multiqc'],
    shell:
        """
        multiqc -f -o {params.output_dir} -n multiqc_report.html {params.mqc_dirs} 2> {log}
        """

rule percent_covered:
    input:
        all = expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_table.txt', samples=SAMPLES.itertuples()),
        q40 = expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_table.txt', samples=SAMPLES.itertuples()),
    output:
        out = f'{output_directory}/analysis/percent_genome_covered.pdf',
    log:
        f'{output_directory}/logs/percent_covered.log',
    benchmark:
        f'{output_directory}/benchmarks/percent_covered.txt',
    threads: 1
    resources:
        mem_gb=config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['short'],
    conda:
        '../envs/python_packages.yaml'
    script:
        '../scripts/plot_percent_covered.py'

if config['control_vectors']:
    rule methylation_controls_qc:
        input:
            bed = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
        output:
            lambda_bed = f'{output_directory}/analysis/qc_vectors/lambda/{{sample}}.bed',
            puc19_bed = f'{output_directory}/analysis/qc_vectors/puc19/{{sample}}.bed',
        params:
            lambda_dir = f'{output_directory}/analysis/qc_vectors/lambda',
            puc19_dir = f'{output_directory}/analysis/qc_vectors/puc19',
        log:
            lambda_log = f'{output_directory}/logs/qc_vectors/lambda.{{sample}}_QC.log',
            puc19_log = f'{output_directory}/logs/qc_vectors/puc19.{{sample}}_QC.log',
        benchmark:
            f'{output_directory}/benchmarks/qc_vectors/{{sample}}.txt',
        threads: 1
        resources:
            mem_gb=config['hpc_parameters']['small_memory_gb'],
            time = config['runtime']['medium'],
        conda:
            '../envs/biscuit.yaml'
        envmodules:
            config['envmodules']['samtools'],
            config['envmodules']['htslib'],
        shell:
            """
            mkdir -p {params.lambda_dir}
            mkdir -p {params.puc19_dir}
           
            # >J02459.1 Escherichia phage Lambda, complete genome - UNMETHYLATED CONTROL
            zcat {input.bed} | {{ grep '^J02459.1' || true; }} > {output.lambda_bed} 2> {log.lambda_log}

            # >M77789.2 Cloning vector pUC19, complete sequence - METHYLATED CONTROL
            zcat {input.bed} | {{ grep '^M77789.2' || true; }}  > {output.puc19_bed} 2> {log.puc19_log}
            """

    rule methylation_controls_figure:
        input:
            lambda_files = expand(f'{output_directory}/analysis/qc_vectors/lambda/{{samples.sample}}.bed', samples=SAMPLES.itertuples()),
            puc19_files = expand(f'{output_directory}/analysis/qc_vectors/puc19/{{samples.sample}}.bed', samples=SAMPLES.itertuples()),
        output:
            pdf = f'{output_directory}/analysis/qc_vectors/control_vector_boxplot.pdf',
        log:
            fn = f'{output_directory}/logs/qc_vectors/control_vector_boxplot.log',
        benchmark:
            f'{output_directory}/benchmarks/qc_vectors/control_vector_boxplot.txt',
        threads: 1
        resources:
            mem_gb=config['hpc_parameters']['small_memory_gb'],
            time = config['runtime']['short'],
        conda:
            '../envs/r.yaml',
        envmodules:
            config['envmodules']['R'],
        script:
            '../scripts/control_vector.R'
