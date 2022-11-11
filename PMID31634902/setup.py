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

base_dir = Path(__file__).parent

targets_dir = base_dir / 'targets'
targets_dir.mkdir(exist_ok=True)

data_dir = base_dir / 'data'
data_dir.mkdir(exist_ok=True)

run_table = pd.read_csv(base_dir / 'documents' / 'PRJNA565979.txt')
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

    for line in open(base_dir / 'documents' / 'primers.txt'):
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
    nicking_sgRNAs = pd.read_csv(base_dir / 'documents' / 'sgRNAs.txt',
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
    scaffolds = hits.fasta.to_dict(base_dir / 'documents' / 'scaffolds.fasta')

    pegRNAs_in = pd.read_csv(base_dir / 'documents' / 'pegRNAs.txt',
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

def setup_Fig3B(download=True):
    rows = get_fig_rows('Fig3b')

    fig_dir = data_dir / 'Fig3B'
    fig_dir.mkdir(exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    def parse_fn(fn):
        _, target, nick_offset, rep = fn.split('.')[0].split('_')
        rep = int(rep[-1])
        return target, nick_offset, rep

    groups = defaultdict(dict)

    group_descriptions = {}

    for SRR_accession, original_fastq_fn in rows.items():
        target, nick_offset, rep = parse_fn(original_fastq_fn)

        int_nick_offset = nick_offset
        if int_nick_offset == 'none':
            int_nick_offset = 0
        int_nick_offset = int(int_nick_offset)

        group_name = f'{target}_{int_nick_offset:+04d}'
        exp_name = f'{group_name}_{rep}'
        
        groups[group_name][exp_name] = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
            'target_info': target,
            'nick_offset': int_nick_offset,
        }
        
        group_descriptions[group_name] = {
            'target_info': target,
            'supplemental_indices': '',
            'experiment_type': 'prime_editing',
            'pegRNAs': f'{target}_3b',
            'sgRNA': f'{target}_3b_{nick_offset}' if nick_offset != 'none' else '',
        }


    group_descriptions_df = pd.DataFrame.from_dict(group_descriptions, orient='index')
    group_descriptions_df.index.name = 'group'
    group_descriptions_df.sort_index(inplace=True)

    group_descriptions_df.to_csv(fig_dir / 'group_descriptions.csv')

    sample_sheet = {}

    for g_i, (group_name, group) in enumerate(groups.items()):
        for exp_name, info in group.items():
            sample_sheet[exp_name] = {
                'group': group_name,
                'color': g_i,
                **info,
            }

    sample_sheet_df = pd.DataFrame.from_dict(sample_sheet, orient='index')
    sample_sheet_df.index.name = 'sample_name'
    sample_sheet_df.sort_values(by=['target_info', 'nick_offset', 'replicate'], inplace=True)
    sample_sheet_df.drop(columns=['target_info', 'nick_offset'], inplace=True)

    sample_sheet_df['color'] = 0

    for g_i, (group, rows) in enumerate(sample_sheet_df.groupby('group', sort=False)):
        sample_sheet_df.loc[sample_sheet_df['group'] == group, 'color'] = g_i + 1

    sample_sheet_df.to_csv(fig_dir / 'sample_sheet.csv')

def setup_Fig4F(download=True):
    rows = get_fig_rows('Fig4_f')

    fig_dir = data_dir / 'Fig4F'
    fig_dir.mkdir(exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    sgRNAs = pd.read_csv(targets_dir / 'sgRNAs.csv', index_col='name').squeeze()

    def parse_fn(fn):
        _, _, locus, edit, rep = fn.split('.')[0].split('_')
        rep = int(rep[-1])
        return locus, edit, rep

    def group_name_to_locus(group_name):
        locus, edit = group_name.split('_')
        return locus

    def group_name_to_pegRNA_name(group_name):
        locus, edit = group_name.split('_')
        pegRNA_name = f'{locus}_4f_{edit}'
        return pegRNA_name

    def group_name_to_sgRNA_name(group_name):
        locus, edit = group_name.split('_')
        matches = [n for n in sgRNAs.index if n.startswith(f'{locus}_4')]

        if len(matches) != 1:
            seqs = set(sgRNAs.loc[matches])
            if len(seqs) > 1:
                print('Warning: more than one possible nicking sgRNA')

        sgRNA_name = matches[0]

        return sgRNA_name

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        locus, edit, rep = parse_fn(original_fastq_fn)
        
        group_name = f'{locus}_{edit}'
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
            'target_info': group_name_to_locus(group_name),
            'pegRNAs': group_name_to_pegRNA_name(group_name),
            'sgRNA': group_name_to_sgRNA_name(group_name),
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

    #make_targets()

def setup_Fig4G(download=True):
    rows = get_fig_rows('Fig4_g')

    fig_dir = data_dir / 'Fig4G'
    fig_dir.mkdir(exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

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
            'sgRNA': 'HEK3_4a_+90',
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

def setup_Fig4H(download=True):
    rows = get_fig_rows('Fig4_h')

    fig_dir = data_dir / 'Fig4H'
    fig_dir.mkdir(exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    sgRNAs = pd.read_csv(targets_dir / 'sgRNAs.csv', index_col='name').squeeze()
    pegRNAs = pd.read_csv(targets_dir / 'pegRNAs.csv', index_col='name').squeeze()

    def parse_fn(fn):
        fields = fn.split('.')[0].split('_')
        locus = fields[2]
        rep = fields[-1] 
        edit = '_'.join(fields[3:-1])
        rep = int(rep[-1])
        return locus, edit, rep

    def group_name_to_locus(group_name):
        locus, edit = group_name.split('_', 1)
        return locus

    def group_name_to_pegRNA_name(group_name):
        locus, edit = group_name.split('_', 1)
        substrings_to_expand = [
            '1AC',
            '6GT',
            '6GA',
            '5GT',
            '2GC',
            '5GC',
            '6GT',
            '1CA',
            '5GT',
        ]

        for substring in substrings_to_expand:
            edit = edit.replace(substring, f'{substring[:2]}to{substring[2:]}')

        edit = edit.replace('2AAin', '2AAins')

        pegRNA_name = f'{locus}_4h_{edit}'

        if pegRNA_name not in pegRNAs.index:
            raise ValueError(pegRNA_name)

        return pegRNA_name

    def group_name_to_sgRNA_name(group_name):
        locus, edit = group_name.split('_', 1)
        matches = [n for n in sgRNAs.index if n.startswith(f'{locus}_4')]

        if len(matches) != 1:
            seqs = set(sgRNAs.loc[matches])
            if len(seqs) > 1:
                print('Warning: more than one possible nicking sgRNA')

        sgRNA_name = matches[0]

        return sgRNA_name

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        locus, edit, rep = parse_fn(original_fastq_fn)
        
        group_name = f'{locus}_{edit}'
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
            'target_info': group_name_to_locus(group_name),
            'pegRNAs': group_name_to_pegRNA_name(group_name),
            'sgRNA': group_name_to_sgRNA_name(group_name),
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

    #make_targets()

def make_targets():
    batches = rs.arrayed_experiment_group.get_all_batches(base_dir)
    all_groups = pd.concat({name: batch.group_descriptions for name, batch in batches.items()})

    targets = {}

    for ti_name, rows in all_groups.groupby('target_info'):
        pegRNAs = set()
        
        for pegRNA in rows['pegRNAs']:
            pegRNAs.add(pegRNA)

        sgRNAs = set()
        for sgRNA in rows['sgRNA']:
            if sgRNA is not None:
                sgRNAs.add(sgRNA)

        targets[ti_name] = {
            'genome': '',
            'amplicon_primers': ti_name,
            'pegRNAs': ';'.join(sorted(pegRNAs)),
            'sgRNA_sequence': ';'.join(sorted(sgRNAs)),
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
    #setup_Fig3B()

    make_targets()