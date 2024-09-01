import os
import sys
import pandas as pd

# samples
df = pd.read_csv("config/samples.tsv", sep="\t", comment='#')
df = df[df['Analysis'] == 1]

batches = {}
batch_samples = {}
for i, row in df.iterrows():
    batch = row['FamilyID']
    sample = row['SampleID']
    if batch not in batches:
        batches[batch] = {}
    assert sample not in batches[batch]
    batches[batch][sample]= {
        'fq1': row['FqR1'],
        'fq2': row['FqR2'],
        'bam': row['Bam'],
        'batch': row['FamilyID']
    }
    if batch not in batch_samples:
        batch_samples[batch] = [sample]
    else:
        batch_samples[batch].append(sample)


def get_fastq(wildcards):
    s_id = wildcards.sample
    f_id = wildcards.batch
    s_info = batches[f_id][s_id]
    fq = []
    for i in ('fq1', 'fq2'):
        if pd.isna(s_info[i]) or s_info[i] == '':
            pass
        else:
            fq.append(s_info[i])
    return fq

def get_bam(wildcards):
    s_id = wildcards.sample
    f_id = wildcards.batch
    s_info = batches[f_id][s_id]
    bam = ''
    if pd.isna(s_info['bam']) or s_info['bam'] == '':
        pass
    else:
        bam = s_info['bam']
    return bam

def get_rename_fq(wildcards, outdir):
    s_id = wildcards.sample
    f_id = wildcards.batch
    s_info = batches[f_id][s_id]
    fq = []
    for i in ('fq1', 'fq2'):
        if pd.isna(s_info[i]) or s_info[i] == '':
            pass
        else:
            fq.append(os.path.join(outdir, f"{f_id}_{s_id}_{i}.fastq.gz"))
    return fq
