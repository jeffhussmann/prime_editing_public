import itertools
import logging
import multiprocessing
import shlex
import subprocess

from pathlib import Path
from collections import defaultdict, Counter

import pandas as pd
import tqdm

import hits.fasta
import hits.fastq
from hits import utilities

import knock_knock.build_targets

import repair_seq as rs

base_dir = Path(__file__).parent

targets_dir = base_dir / 'targets'
targets_dir.mkdir(exist_ok=True)

data_dir = base_dir / 'data'
data_dir.mkdir(exist_ok=True)

run_table = pd.read_csv(base_dir / 'PRJNA770428.txt')
SRR_to_fn = run_table[['Run', 'Sample Name']].set_index('Run').squeeze()

def get_fig_rows(fig_prefix):
    return SRR_to_fn[SRR_to_fn.str.startswith(fig_prefix, na=False)].sort_values()

def download_SRR(SRR_accession, dest_dir):
    subprocess.run(shlex.split(f'fastq-dump --gzip --outdir {dest_dir} {SRR_accession}'), check=True)

def download_rows(rows, dest_dir):
    with multiprocessing.Pool(processes=12) as pool:
        pool.starmap(download_SRR, [(SRR, dest_dir) for SRR in rows.index])

def load_primers():
    amplicon_primers = pd.read_csv(base_dir / 'amplicon_primers.txt', sep=' ', comment='#', header=None)

    prefixes = {
        'fwd': 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNN',
        'rev': 'TGGAGTTCAGACGTGTGCTCTTCCGATCT',
    }

    relevant_primers = {}

    for name, rows in amplicon_primers.groupby(0):
        seqs = set(rows[1])
        if len(seqs) > 1:
            raise ValueError(rows)
        else:
            seq = sorted(seqs)[0]
        
        if seq.startswith(prefixes['fwd']):
            relevant_seq = seq[len(prefixes['fwd']):]
        elif seq.startswith(prefixes['rev']):
            relevant_seq = seq[len(prefixes['rev']):]
        else:
            raise ValueError(seq)
            
        relevant_primers[name] = relevant_seq

    relevant_primers = pd.Series(relevant_primers)

    if len(set(relevant_primers.index)) != len(relevant_primers.index):
        raise ValueError

    return relevant_primers

def infer_primers(fastq_fn, num_reads=None):
    name_to_primer = load_primers()
    name_to_primer_rc = {name: utilities.reverse_complement(seq) for name, seq in name_to_primer.items()}

    total = 0
    
    starts = Counter()
    ends = Counter()

    reads = hits.fastq.reads(fastq_fn)

    if num_reads is not None:
        reads = itertools.islice(reads, num_reads)

    for read in reads:
        total += 1
        read = read[4:]
        for name, seq in name_to_primer.items():
            seq_rc = name_to_primer_rc[name]
            
            if read.seq.startswith(seq):
                starts[name] += 1
                
            if read.seq.endswith(seq_rc):
                ends[name] += 1
    
    if len(starts) > 0:
        (left_primer, left_count), = starts.most_common(1)
    else:
        left_primer, left_count = None, 0

    if len(ends) > 0:
        (right_primer, right_count), = ends.most_common(1)
    else:
        right_primer, right_count = None, 0
    
    left_fraction = left_count / total
    right_fraction = right_count / total
    
    return left_primer, right_primer, left_fraction, right_fraction, total

def infer_all_primers(fig_dir, rows, num_reads=None):
    SRR_to_primers = {}

    left_fractions = []
    right_fractions = []

    for SRR_accession in tqdm.tqdm(rows.index, desc='Inferring amplicon primers'):
        fastq_fn = fig_dir / f'{SRR_accession}.fastq.gz'
        
        left_primer, right_primer, left_fraction, right_fraction, total = infer_primers(fastq_fn, num_reads)

        left_fractions.append(left_fraction)
        right_fractions.append(right_fraction)
        
        SRR_to_primers[SRR_accession] = (left_primer, right_primer)

    return SRR_to_primers

def setup_Fig1C(download=True):
    rows = get_fig_rows('Fig1c')

    fig_dir = data_dir / 'Fig1C'
    fig_dir.mkdir(parents=True, exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    # Can't infer primers because there is an early cycle that is
    # entirely N in relevant fastqs.

    name_to_primer = load_primers()

    primer_pairs_seen = {'HEK3': ';'.join([name_to_primer['HEK3_fwd'], name_to_primer['HEK3_rev']])}

    primer_pairs_seen = pd.Series(primer_pairs_seen)
    primer_pairs_seen.index.name = 'name'
    primer_pairs_seen.name = 'amplicon_primer_sequences'

    primer_pair_fn = fig_dir / 'amplicon_primers.csv'
    primer_pairs_seen.to_csv(primer_pair_fn)

    def parse_fn(fn):
        _, att_type, A_length, B_length, rep = fn.split('.')[0].split('_')
        A_length = int(A_length[1:])
        B_length = int(B_length[1:])
        rep = int(rep[-1])
        return att_type, A_length, B_length, rep

    def group_name_to_pegRNA_names(group_name):
        _, att_type, A, B = group_name.split('_')
        pegRNA_names = [
            f'HEK3_{att_type}_A_{A[1:]}',
            f'HEK3_{att_type}_B_{B[1:]}',
        ]
        return pegRNA_names

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        att_type, A_length, B_length, rep = parse_fn(original_fastq_fn)
        group_name = f'HEK3_{att_type}_A{A_length}_B{B_length}'
        exp_name = f'{group_name}_{rep}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'twin_prime',
            'target_info': 'HEK3',
            'pegRNAs': ';'.join(group_name_to_pegRNA_names(group_name)),
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

def setup_Fig2C(download=True):
    rows = get_fig_rows('Fig2c')

    fig_dir = data_dir / 'Fig2C'
    fig_dir.mkdir(parents=True, exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    def primer_pair_to_target_name(first, second):
        return f'PAH-{first[-4:]}+{second[-4:]}'

    # Note: PAH exon 4 samples appear to use an unannotated primer.
    SRR_to_primers = infer_all_primers(fig_dir, rows, num_reads=1000)
    name_to_primer = load_primers()

    primer_pairs_seen = {}
    for first, second in SRR_to_primers.values():
        name = primer_pair_to_target_name(first, second)
        seqs = ';'.join([name_to_primer[first], name_to_primer[second]])
        primer_pairs_seen[name] = seqs

    primer_pairs_seen = pd.Series(primer_pairs_seen)
    primer_pairs_seen.index.name = 'name'
    primer_pairs_seen.name = 'amplicon_primer_sequences'

    primer_pair_fn = fig_dir / 'amplicon_primers.csv'
    primer_pairs_seen.to_csv(primer_pair_fn)

    exon_and_overlap_to_pegRNA_prefixes = {
        (4, 24): ['PAH_E4.2_45', 'PAH_E4.4_43'],
        (4, 36): ['PAH_E4.2_50', 'PAH_E4.4_50'],
        (4, 59): ['PAH_E4.2_62', 'PAH_E4.4_61'],
        (7, 22): ['PAH_E7.2_34', 'PAH_E7.5_34'],
        (7, 24): ['PAH_E7.2_44', 'PAH_E7.6_44'],
        (7, 42): ['PAH_E7.2_44', 'PAH_E7.5_44'],
        (7, 47): ['PAH_E7.2_55', 'PAH_E7.6_56'],
    }

    def parse_Fig2C_fn(fn):
        fn = fn.split('.')[0]
        _, _, exon, overlap, _, pegRNA_type, _ = fn.split('_', 6)
        
        if exon == 'Ex4':
            recoded = 64
        else:
            if exon[-1] == 'a':
                recoded = 46
            elif exon[-1] == 'b':
                recoded = 64
            else:
                raise ValueError(exon)
                
        exon_num = int(exon[2])
        
        rep = int(fn[-1])
        
        overlap = int(overlap[:-2])
        
        return exon_num, recoded, overlap, pegRNA_type, rep

    def Fig2C_group_name_to_pegRNA_names(group_name):
        target_name, recoded, overlap, pegRNA_type = group_name.split('_')
        
        exon_num = 4 if '1687' in target_name else 7
        overlap = int(overlap)
        
        extra = '_EvoPreQ1' if pegRNA_type == 'epegRNA' else ''
        
        pegRNA_names = [f'{name}{extra}' for name in exon_and_overlap_to_pegRNA_prefixes[exon_num, overlap]]
        
        return pegRNA_names

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        exon, recoded, overlap, pegRNA_type, rep = parse_Fig2C_fn(original_fastq_fn)
        
        target_name = primer_pair_to_target_name(*SRR_to_primers[SRR_accession])
        
        group_name = f'{target_name}_{recoded}_{overlap}_{pegRNA_type}'
        
        exp_name = f'{group_name}_{rep}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'twin_prime',
            'target_info': group_name.split('_')[0],
            'pegRNAs': ';'.join(Fig2C_group_name_to_pegRNA_names(group_name)),
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

def setup_Fig2E(download=True):
    rows = get_fig_rows('Fig2e')

    fig_dir = data_dir / 'Fig2E'
    fig_dir.mkdir(parents=True, exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    name_to_primer = load_primers()

    primer_pairs_seen = {'HEK3': ';'.join([name_to_primer['HEK3_fwd'], name_to_primer['HEK3_rev']])}

    primer_pairs_seen = pd.Series(primer_pairs_seen)
    primer_pairs_seen.index.name = 'name'
    primer_pairs_seen.name = 'amplicon_primer_sequences'

    primer_pair_fn = fig_dir / 'amplicon_primers.csv'
    primer_pairs_seen.to_csv(primer_pair_fn)

    def parse_Fig2E_fn(fn):
        _, _, pegRNA_type, _, strategy, length, rep, _ = fn.split('.')[0].split('-')
        rep = int(rep[-1])
        length = int(length[3:-2])
        return pegRNA_type, strategy, length, rep

    def Fig2E_group_name_to_pegRNA_names(group_name):
        _, strategy, length, pegRNA_type = group_name.split('_')
        extra = '_EvoPreQ1' if pegRNA_type == 'EvoPreQ1' else ''
        pegRNA_names = [f'HEK3_DF_{letter}_{strategy}_del{length}nt{extra}' for letter in ['A', 'B']]
        return pegRNA_names

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        pegRNA_type, strategy, length, rep = parse_Fig2E_fn(original_fastq_fn)
        group_name = f'HEK3_{strategy}_{length}_{pegRNA_type}'
        exp_name = f'{group_name}_{rep}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'twin_prime',
            'target_info': 'HEK3',
            'pegRNAs': ';'.join(Fig2E_group_name_to_pegRNA_names(group_name)),
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

def setup_Fig3C(download=True):
    rows = get_fig_rows('Fig3c')

    fig_dir = data_dir / 'Fig3C'
    fig_dir.mkdir(parents=True, exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    def parse_Fig3C_fn(fn):
        fn = fn.split('.')[0]
        _, _, _, pegRNAs, rep, _ = fn.split('-')
        
        first_pegRNA, second_pegRNA = pegRNAs.split('_')
        
        # 22.07.05: empirically, there appears to be a cyclic mis-annotation of the second pegRNA
        # for samples with first_pegRNA A1615a.
        
        if first_pegRNA == 'A1615a':
            fix_annotations = {
                'B1701a': 'B1705b',
                'B1676a': 'B1701a',
                'B1705b': 'B1676a',
            }
            second_pegRNA = fix_annotations[second_pegRNA]

        rep = int(rep[-1])
        
        return first_pegRNA, second_pegRNA, rep

    def Fig3C_group_name_to_pegRNA_names(group_name):
        target_name, first_pegRNA, second_pegRNA = group_name.split('_')
            
        pegRNA_names = [f'AAVS1_{name}' for name in [first_pegRNA, second_pegRNA]]
        
        return pegRNA_names

    def primer_pair_to_target_name(first, second):
        return f'AAVS1-{first[-4:]}+{second[-4:]}'

    SRR_to_primers = infer_all_primers(fig_dir, rows)
    name_to_primer = load_primers()

    primer_pairs_seen = {}
    for first, second in SRR_to_primers.values():
        name = primer_pair_to_target_name(first, second)
        seqs = ';'.join([name_to_primer[first], name_to_primer[second]])
        primer_pairs_seen[name] = seqs

    primer_pairs_seen = pd.Series(primer_pairs_seen)
    primer_pairs_seen.index.name = 'name'
    primer_pairs_seen.name = 'amplicon_primer_sequences'

    primer_pair_fn = fig_dir / 'amplicon_primers.csv'
    primer_pairs_seen.to_csv(primer_pair_fn)

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        first_pegRNA, second_pegRNA, rep = parse_Fig3C_fn(original_fastq_fn)
        
        target_name = primer_pair_to_target_name(*SRR_to_primers[SRR_accession])
        
        group_name = f'{target_name}_{first_pegRNA}_{second_pegRNA}'
        
        exp_name = f'{group_name}_{rep}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': rep,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'twin_prime',
            'target_info': group_name.split('_')[0],
            'pegRNAs': ';'.join(Fig3C_group_name_to_pegRNA_names(group_name)),
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

def setup_FigED4(download=True):
    rows = get_fig_rows('ED_Fig4')

    fig_dir = data_dir / 'FigED4'
    fig_dir.mkdir(parents=True, exist_ok=True)

    if download:
        download_rows(rows, fig_dir)

    def primer_pair_to_target_name(first, second):
        return f'AAVS1-{first[-4:]}+{second[-4:]}'

    SRR_to_primers = infer_all_primers(fig_dir, rows)
    name_to_primer = load_primers()

    primer_pairs_seen = {}
    for first, second in SRR_to_primers.values():
        name = primer_pair_to_target_name(first, second)
        seqs = ';'.join([name_to_primer[first], name_to_primer[second]])
        primer_pairs_seen[name] = seqs

    primer_pairs_seen = pd.Series(primer_pairs_seen)
    primer_pairs_seen.index.name = 'name'
    primer_pairs_seen.name = 'amplicon_primer_sequences'

    primer_pair_fn = fig_dir / 'amplicon_primers.csv'
    primer_pairs_seen.to_csv(primer_pair_fn)

    def parse_FigED4_fn(fn):
        fn = fn.split('.')[0]
        _, _, _, first_pegRNA, second_pegRNA, _ = fn.split('_')
        
        return first_pegRNA, second_pegRNA

    def FigED4_group_name_to_pegRNA_names(group_name):
        target_name, first_pegRNA, second_pegRNA = group_name.split('_')
            
        pegRNA_names = [f'AAVS1_{name}' for name in [first_pegRNA, second_pegRNA]]
        
        return pegRNA_names

    groups = defaultdict(dict)

    for SRR_accession, original_fastq_fn in rows.items():
        first_pegRNA, second_pegRNA = parse_FigED4_fn(original_fastq_fn)
        
        target_name = primer_pair_to_target_name(*SRR_to_primers[SRR_accession])
        
        group_name = f'{target_name}_{first_pegRNA}_{second_pegRNA}'
        
        exp_name = f'{group_name}'
        
        info = {
            'R1': f'{SRR_accession}.fastq.gz',
            'replicate': 1,
        }
        
        groups[group_name][exp_name] = info

    group_descriptions = {
        group_name: {
            'supplemental_indices': '',
            'experiment_type': 'twin_prime',
            'target_info': group_name.split('_')[0],
            'pegRNAs': ';'.join(FigED4_group_name_to_pegRNA_names(group_name)),
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

    # Merge amplicon primers from each batch.
    fns = [batch.data_dir / 'amplicon_primers.csv' for batch in batches.values()]
    all_dfs = [pd.read_csv(fn, index_col=0).squeeze('columns') for fn in fns]
    amplicon_primers = pd.concat(all_dfs).drop_duplicates()
    amplicon_primers.to_csv(targets_dir / 'amplicon_primers.csv')

    all_groups = pd.concat({name: batch.group_descriptions for name, batch in batches.items()})

    targets = {}

    for ti_name, rows in all_groups.groupby('target_info'):
        pegRNAs = set()
        
        for pegRNA_pair in rows['pegRNAs']:
            for pegRNA in pegRNA_pair.split(';'):
                pegRNAs.add(pegRNA)

        targets[ti_name] = {
            'genome': 'hg38',
            'amplicon_primers': ti_name,
            'pegRNAs': ';'.join(sorted(pegRNAs)),
        }

    targets_df = pd.DataFrame(targets).T

    targets_df.index.name = 'name'
    targets_df = targets_df.sort_index()[['genome', 'amplicon_primers', 'pegRNAs']]

    targets_df.to_csv(targets_dir / 'targets.csv')

    knock_knock.build_targets.build_target_infos_from_csv(base_dir)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    make_targets()