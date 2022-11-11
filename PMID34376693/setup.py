import multiprocessing
import shlex
import subprocess
from pathlib import Path

import pandas as pd

base_dir = Path(__file__).parent

targets_dir = base_dir / 'targets'
targets_dir.mkdir(exist_ok=True)

data_dir = base_dir / 'data'
data_dir.mkdir(exist_ok=True)

run_table = pd.read_csv(base_dir / 'documents' / 'PRJNA641538.txt')
SRR_to_sample_name = run_table[['Run', 'Sample Name']].set_index('Run').squeeze()
excel_table = pd.read_excel(base_dir / 'documents' / '41467_2021_25154_MOESM7_ESM.xlsx')

def download_SRR(SRR_accession, dest_dir):
    subprocess.run(shlex.split(f'fastq-dump --split-3 --gzip --outdir {dest_dir} {SRR_accession}'), check=True)

def download_rows(rows, dest_dir):
    with multiprocessing.Pool(processes=4) as pool:
        pool.starmap(download_SRR, [(SRR, dest_dir) for SRR in rows.index])

def setup_Fig1(download=True):
    fig_dir = data_dir / 'Fig1'
    fig_dir.mkdir(exist_ok=True)

    rows = SRR_to_sample_name[SRR_to_sample_name.str.startswith('HPRT')]

    if download:
        download_rows(rows, fig_dir)

    excel_rows = excel_table[excel_table['name'].str.startswith('HPRT')].set_index('name')

    df = pd.DataFrame(rows)
    df['genetic_background'] = [excel_rows.loc[name, 'genetic background'] for name in df['Sample Name']]

    groups = {}
    samples = {}

    groups['eGFP'] = {
        'supplemental_indices': '',
        'experiment_type': 'prime_editing',
        'target_info': 'eGFP',
        'sgRNAs': 'eGFP_A;eGFP_B',
        'min_relevant_length': 100,
        'condition_keys': 'genetic_background',
        'baseline_condition': 'wild-type',
    }

    for genetic_background, gb_rows in df.groupby('genetic_background'):
        for rep_i, (SRR, row) in enumerate(gb_rows.iterrows(), 1):
            group_name = 'eGFP'

            samples[row['Sample Name']] = {
                'R1': f'{SRR}_1.fastq.gz',
                'R2': f'{SRR}_2.fastq.gz',
                'group': group_name,
                'genetic_background': genetic_background,
                'replicate': rep_i,
            }

    groups_df = pd.DataFrame.from_dict(groups, orient='index')
    groups_df.index.name = 'group'

    groups_csv_fn = fig_dir / 'group_descriptions.csv'
    groups_df.to_csv(groups_csv_fn)

    samples_df = pd.DataFrame.from_dict(samples, orient='index')
    samples_df.index.name = 'sample_name'

    samples_csv_fn = fig_dir / 'sample_sheet.csv'
    samples_df.to_csv(samples_csv_fn)

    return samples_df, groups_df
