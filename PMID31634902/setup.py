import multiprocessing
import logging
import shlex
import subprocess

from pathlib import Path
from collections import defaultdict

import pandas as pd

import hits.fasta

import repair_seq as rs
import knock_knock.build_targets

paper_name = 'PMID31634902'

base_dir = Path('/home/jah/projects/prime_editing_public') / paper_name

targets_dir = base_dir / 'targets'
targets_dir.mkdir(exist_ok=True)

data_dir = base_dir / 'data'
data_dir.mkdir(exist_ok=True)

run_table = pd.read_csv(base_dir / 'PRJNA565979.txt')
SRR_to_fn = run_table[['Run', 'Sample Name']].set_index('Run').squeeze()

def get_fig_rows(fig_prefix):
    return SRR_to_fn[SRR_to_fn.str.startswith(fig_prefix, na=False)].sort_values()

def download_SRR(SRR_accession, dest_dir):
    subprocess.run(shlex.split(f'fastq-dump --gzip --outdir {dest_dir} {SRR_accession}'), check=True)

def download_rows(rows, dest_dir):
    with multiprocessing.Pool(processes=12) as pool:
        pool.starmap(download_SRR, [(SRR, dest_dir) for SRR in rows.index])

def make_primers():
    full_primers = defaultdict(dict)

    for line in open(base_dir / 'primers.txt'):
        name, side, sequence = line.strip().rsplit(' ', 2)
        if 'off-target' in name:
            continue
        full_primers[name][side] = sequence

    prefixes = {
        'fwd': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNN',
        'rev': 'TGGAGTTCAGACGTGTGCTCTTCCGATCT',
    }

    primers = defaultdict(dict)

    for name in full_primers:
        for direction in ['fwd', 'rev']:
            full_primer = full_primers[name][direction]
            prefix = prefixes[direction]
            assert full_primer.startswith(prefix)
            primers[name][direction] = full_primer[len(prefixes[direction]):]

    with open(targets_dir / 'amplicon_primers.csv', 'w') as fh:
        fh.write('name,amplicon_primer_sequences\n')
        for name in sorted(primers):
            fh.write(f'{name},{primers[name]["fwd"]};{primers[name]["rev"]}\n')

def make_sgRNAs():
    nicking_sgRNAs = pd.read_csv(base_dir / 'sgRNAs.txt',
                                 sep=' ',
                                 comment='#',
                                 header=None,
                                 index_col=0,
                                ).squeeze().sort_index()


    with open(targets_dir / 'sgRNAs.csv', 'w') as fh:
        fh.write('name,sgRNA_sequence\n')
                
        for name, seq in sorted(set(nicking_sgRNAs.items())):
            fh.write(f'{name},{seq}\n')

def make_pegRNAs():
    scaffolds = hits.fasta.to_dict(base_dir / 'scaffolds.fasta')

    pegRNAs_in = pd.read_csv(base_dir / 'pegRNAs.txt',
                             sep=' ',
                             header=None,
                             index_col='name',
                             names=['name', 'spacer', 'extension', 'PBS_length', 'RT_template_length'],
                             comment='#',
                            )

    rows = {}

    for full_name, row in pegRNAs_in.iterrows():
        
        name = full_name.lstrip('*^')
        
        if full_name.startswith('*'):
            scaffold = scaffolds['sgRNA_scaffold_2']
        elif full_name.startswith('^'):
            # 22.07.17 - this appears to have been lost,
            # no line start with this.
            scaffold = scaffolds['sgRNA_scaffold_3']
        else:
            scaffold = scaffolds['sgRNA_scaffold_1']
            
        rows[name] = {
            'effector': 'SpCas9H840A',
            'protospacer': row['spacer'],
            'scaffold': scaffold,
            'extension': row['extension'],
        }

    pegRNAs = pd.DataFrame.from_dict(rows, orient='index')
    pegRNAs.index.name = 'name'

    pegRNAs.to_csv(targets_dir / f'pegRNAs.csv')

def setup_Fig4G():
    rows = get_fig_rows('Fig4_g')

    fig_dir = data_dir / 'Fig4G'
    fig_dir.mkdir(exist_ok=True)

    #download_rows(rows, fig_dir)

    def parse_fn(fn):
        _, _, del_type, rep = fn.split('.')[0].split('_')
        rep = int(rep[-1])
        return del_type, rep

    def group_name_to_pegRNA_name(group_name):
        _, del_type = group_name.split('_')
        pegRNA_name = f'HEK3_4g_{del_type}'
        return pegRNA_name

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        del_type, rep = parse_fn(original_fastq_fn)
        
        group_name = f'HEK3_{del_type}'
        exp_name = f'{group_name}_{rep}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'prime_editing',
            'target_info': 'HEK3',
            'pegRNAs': group_name_to_pegRNA_name(group_name),
            'sgRNA': 'sgRNA_HEK3_4a_+90',
        } for group_name in groups
    }    

    group_descriptions_df = pd.DataFrame.from_dict(group_descriptions, orient='index')
    group_descriptions_df.index.name = 'group'
    group_descriptions_df.sort_index(inplace=True)

    group_descriptions_df.to_csv(fig_dir / 'group_descriptions.csv')

    sample_sheet = {}

    for group_name, group in groups.items():
        for exp_name, info in group.items():
            sample_sheet[exp_name] = {
                'group': group_name,
                **info,
            }

    sample_sheet_df = pd.DataFrame.from_dict(sample_sheet, orient='index')
    sample_sheet_df.index.name = 'sample_name'
    sample_sheet_df.sort_index(inplace=True)

    sample_sheet_df.to_csv(fig_dir / 'sample_sheet.csv')

def make_targets():
    batches = rs.arrayed_experiment_group.get_all_batches(base_dir)
    all_groups = pd.concat({name: batch.group_descriptions for name, batch in batches.items()})

    targets = {}

    sgRNAs = pd.read_csv(targets_dir / 'sgRNAs.csv', index_col='name').squeeze()

    for ti_name, rows in all_groups.groupby('target_info'):
        pegRNAs = set()
        
        for pegRNA in rows['pegRNAs']:
            pegRNAs.add(pegRNA)

        targets[ti_name] = {
            'genome': 'hg38',
            'amplicon_primers': ti_name,
            'pegRNAs': ';'.join(sorted(pegRNAs)),
            'sgRNA_sequence': ';'.join(sorted(n for n in sgRNAs.index if n.startswith(ti_name))),
        }

    targets_df = pd.DataFrame(targets).T
    targets_df.index.name = 'name'
    targets_df = targets_df.sort_index()[['genome', 'amplicon_primers', 'pegRNAs', 'sgRNA_sequence']]

    targets_df.to_csv(targets_dir / 'targets.csv')

    knock_knock.build_targets.build_target_infos_from_csv(base_dir)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    #make_primers()
    #make_sgRNAs()
    #make_pegRNAs()

    #setup_Fig4G()

    make_targets()