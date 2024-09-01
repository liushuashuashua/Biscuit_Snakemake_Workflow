#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import sys
import pandas as pd

# samples

def rename_fastqs(samplesheet, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    df = pd.read_csv(samplesheet, sep="\t", comment='#')
    df = df[df['Analysis'] == 1]
    for i, row in df.iterrows():
        f_id = row['FamilyID']
        s_id = row['SampleID']
        fq1 = row['FqR1']
        fq2 = row['FqR2']
        if pd.isna(fq1) and pd.isna(fq2):
            continue
        elif pd.isna(fq1):
            print(f"Error: {fq1} not found for sample {s_id}")
            sys.exit(1)
        elif pd.isna(fq2):
            assert fq1 is not None and os.path.exists(fq1)
            new_fq_path = os.path.join(outdir, f"{s_id}.fastq.gz")
            print(f"Symlinking {fq1} to {new_fq_path}")
            os.symlink(os.path.abspath(fq1), new_fq_path)
            df.loc[i, 'fq1'] = new_fq_path
        else:
            assert fq1 is not None and fq2 is not None and os.path.exists(fq1) and os.path.exists(fq2)
            new_fq1_path = os.path.join(outdir, f"{s_id}_R1.fastq.gz")
            new_fq2_path = os.path.join(outdir, f"{s_id}_R2.fastq.gz")
            print(f"Symlinking {fq1} to {new_fq1_path} and {fq2} to {new_fq2_path}")
            os.symlink(os.path.abspath(fq1), new_fq1_path)
            os.symlink(os.path.abspath(fq2), new_fq2_path)
            df.loc[i, 'fq1'] = new_fq1_path
            df.loc[i, 'fq2'] = new_fq2_path
    new_samplesheet = os.path.join(outdir, "samplesheet.csv")
    df.to_csv(new_samplesheet, sep="\t", index=False)


rename_fastqs(snakemake.params.samplesheet, snakemake.output.symlink_dir)
