###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in preprocessing.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

checkpoint rename_fastq_files:
    output:
        symlink_dir = directory(f'{output_directory}/analysis/renamed_fastqs'),
    params:
        samplesheet = config['samples'],
    log:
        f'{output_directory}/logs/rename/rename.log',
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['short'],
    conda:
        '../envs/python_packages.yaml'
    script:
        '../scripts/rename.py'
  
def get_renamed_fastq_files(wildcards):
    sample = wildcards.sample
    print(f'Getting renamed fastq files for {sample}')
    checkpoint_sample = checkpoints.rename_fastq_files.get(sample=wildcards.sample)
    if os.path.exists(f'{checkpoint_sample.output.symlink_dir}/{sample}_R1.fastq.gz') and os.path.exists(f'{checkpoint_sample.output.symlink_dir}/{sample}_R2.fastq.gz'):
        return [f'{checkpoint_sample.output.symlink_dir}/{sample}_R1.fastq.gz', f'{checkpoint_sample.output.symlink_dir}/{sample}_R2.fastq.gz']
    else:
        return [f'{checkpoint_sample.output.symlink_dir}/{sample}.fastq.gz']

rule raw_fastqc_pe:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_R1_fastqc.html',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_R1_fastqc.zip',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_R2_fastqc.html',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_R2_fastqc.zip',
    params:
        dir = f'{output_directory}/analysis/raw_fastqc',
    log:
        stdout = f'{output_directory}/logs/raw_fastqc/{{sample}}.o',
        stderr = f'{output_directory}/logs/raw_fastqc/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/raw_fastqc/{{sample}}.txt',
    conda:
        '../envs/babraham.yaml'
    envmodules:
        config['envmodules']['fastqc'],
    threads: config['hpc_parameters']['trim_threads'],
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['medium'],
    shell:
        """
        if [[ -d {params.dir} ]]; then echo "raw_fastqc already exists"; else mkdir -p {params.dir}; fi
        fastqc --outdir {params.dir} --threads {threads} {input} 2> {log.stderr} 1> {log.stdout}
        """

rule raw_fastqc_se:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_fastqc.html',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}_fastqc.zip',
    params:
        dir = f'{output_directory}/analysis/raw_fastqc',
    log:
        stdout = f'{output_directory}/logs/raw_fastqc/{{sample}}.o',
        stderr = f'{output_directory}/logs/raw_fastqc/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/raw_fastqc/{{sample}}.txt',
    conda:
        '../envs/babraham.yaml'
    envmodules:
        config['envmodules']['fastqc'],
    threads: config['hpc_parameters']['trim_threads'],
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['medium'],
    shell:
        """
        if [[ -d {params.dir} ]]; then echo "raw_fastqc already exists"; else mkdir -p {params.dir}; fi
        fastqc --outdir {params.dir} --threads {threads} {input} 2> {log.stderr} 1> {log.stdout}
        """

rule trim_reads_pe:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/trim_reads/{{sample}}_val_1.fq.gz',
        f'{output_directory}/analysis/trim_reads/{{sample}}_val_2.fq.gz',
    params:
        outdir = f'{output_directory}/analysis/trim_reads',
        args_list = config['trim_galore']['args_list'],
    log:
        stdout = f'{output_directory}/logs/trim_reads/{{sample}}.o',
        stderr = f'{output_directory}/logs/trim_reads/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/trim_reads/{{sample}}.txt',
    conda:
        '../envs/babraham.yaml'
    envmodules:
        config['envmodules']['trim_galore'],
        config['envmodules']['fastqc'],
    threads: config['hpc_parameters']['trim_threads']
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['medium'],
    shell:
        """
        trim_galore \
            --paired \
            {input} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample} \
            --cores {threads} \
            --fastqc \
            --gzip \
            {params.args_list} \
            2> {log.stderr} 1> {log.stdout}

       """

rule trim_reads_se:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/trim_reads/{{sample}}_trimmed.fq.gz',
    params:
        outdir = f'{output_directory}/analysis/trim_reads',
        args_list = config['trim_galore']['args_list_se'],
    log:
        stdout = f'{output_directory}/logs/trim_reads/{{sample}}.o',
        stderr = f'{output_directory}/logs/trim_reads/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/trim_reads/{{sample}}.txt',
    conda:
        '../envs/babraham.yaml'
    envmodules:
        config['envmodules']['trim_galore'],
        config['envmodules']['fastqc'],
    threads: config['hpc_parameters']['trim_threads']
    resources:
        mem_gb = config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['medium'],
    shell:
        """
        trim_galore \
            {input} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample} \
            --cores {threads} \
            --fastqc \
            --gzip \
            {params.args_list} \
            2> {log.stderr} 1> {log.stdout}

        """


rule fastq_screen:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/fastq_screen/{{sample}}_R1_screen.html',
        f'{output_directory}/analysis/fastq_screen/{{sample}}_R2_screen.html',
        f'{output_directory}/analysis/fastq_screen/{{sample}}_R1_screen.txt',
        f'{output_directory}/analysis/fastq_screen/{{sample}}_R2_screen.txt',
    params:
        conf = config['fastq_screen']['conf'],
        output_dir = f'{output_directory}/analysis/fastq_screen/',
    log:
        fastq_screen = f'{output_directory}/logs/fastq_screen/{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/fastq_screen/{{sample}}.txt'
    threads: 8
    resources:
        nodes = 1,
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/babraham.yaml'
    envmodules: 
        config['envmodules']['fastq_screen'],
        config['envmodules']['bismark'],
    shell:
        """
        fastq_screen --bisulfite --conf {params.conf} --outdir {params.output_dir} {input} 2> {log.fastq_screen}
        """

rule fastq_screen_se:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/fastq_screen/{{sample}}_screen.html',
        f'{output_directory}/analysis/fastq_screen/{{sample}}_screen.txt',
    params:
        conf = config['fastq_screen']['conf'],
        output_dir = f'{output_directory}/analysis/fastq_screen/',
    log:
        fastq_screen = f'{output_directory}/logs/fastq_screen/{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/fastq_screen/{{sample}}.txt'
    threads: 8
    resources:
        nodes = 1,
        mem_gb = config['hpc_parameters']['intermediate_memory_gb'],
        time = config['runtime']['medium'],
    conda:
        '../envs/babraham.yaml'
    envmodules: 
        config['envmodules']['fastq_screen'],
        config['envmodules']['bismark'],
    shell:
        """
        fastq_screen --bisulfite --conf {params.conf} --outdir {params.output_dir} {input} 2> {log.fastq_screen}
        """


rule prep_chromsizes_file:
    input:
        fai=f'{config["ref"]["fasta"]}.fai'
    output:
        f'{output_directory}/analysis/prep_chromsizes_file/chrom_sizes.tsv'
    log:
        stderr=f'{output_directory}/logs/prep_chromsizes_file/err.e',
    benchmark:
        f'{output_directory}/benchmarks/prep_chromsizes_file/bench.txt'
    threads: 1
    resources:
        nodes = 1,
        mem_gb=config['hpc_parameters']['small_memory_gb'],
        time = config['runtime']['short'],
    shell:
        """
        cut -f 1,2 {input.fai} > {output} 2> {log.stderr}
        """
