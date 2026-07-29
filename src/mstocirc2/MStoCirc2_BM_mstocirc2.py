import time
import sys
import os
import subprocess
import csv
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
import re
import math

BANNER = r"""
 __  __  _____   _____  _____   _____  _____   ___
|  \/  || ____| | _ \  | ____| |_   _|| __ \  |__ \
| \  / || |___  | | | | | |___    | | | | | |__) |  |
| |\/| ||____ | | _  \ _ \ |    | | | | | _ /   / /
| |  | ||____) || | | | | |____ _| |_ | | \ \  / /_
|_|  |_|_____/ \_| |_| |_____|_____||_|  \_\|____|
"""


def clean_peptide_sequence(pep):
    """Unified peptide cleaning engine."""
    pep = pep.strip().upper()
    if '.' in pep:
        match = re.match(r'^[A-Z\-_]\.(.+)\.[A-Z\-_]$', pep)
        if match:
            pep = match.group(1)
    pep = re.sub(r'[^A-Z]', '', pep)
    return pep


def extract_circ_info_id(orf_base):
    """Extract circRNA ID from ORF base ID, stripping known prefixes."""
    PREFIXES = [
        'F_exon_', 'F_base_', 'F_self_', 'F_noexonic_', 'F_unknown_',
        'R_exon_', 'R_base_', 'R_self_', 'R_noexonic_', 'R_unknown_',
        'exon_',
    ]
    for pfx in PREFIXES:
        if orf_base.startswith(pfx):
            return orf_base[len(pfx):]
    return orf_base


def main(argv=None):
    start_time_total = time.time()

    def open_safe_read(file_path):
        try:
            return open(file_path, 'r', encoding='utf-8', errors='ignore')
        except Exception:
            return open(file_path, 'r', encoding='latin-1', errors='ignore')

    class CustomHelpFormatter(argparse.RawTextHelpFormatter):
        def _format_action(self, action):
            parts = super(argparse.RawTextHelpFormatter, self)._format_action(action)
            if action.nargs == argparse.PARSER:
                parts = "\n".join(parts.split("\n")[1:])
            return parts
        def _get_help_string(self, action):
            return action.help

    parser = argparse.ArgumentParser(
        description=BANNER + "\nMScirctl2: Circular RNA Mass Spectrometry Analysis Tool\n"
        "-------------------------------------------------------",
        formatter_class=CustomHelpFormatter, add_help=False
    )
    required_group = parser.add_argument_group('Required Arguments')
    optional_group = parser.add_argument_group('Optional Arguments')
    required_group.add_argument('-fi', '--circ_file_input', metavar='<FILE>', required=False,
        help='Result folder from MScirctl1.\nOptional if -cs and -corf are specified.')
    required_group.add_argument('-cs', '--circ_seq', metavar='<FILE>', required=False,
        help='circRNA full-length sequence file.\nRequired if -fi is not specified.')
    required_group.add_argument('-ci', '--circ_info', metavar='<FILE>', required=False,
        help='circRNA annotation file (.gff or .txt format).\nOptional if -fi is not specified.')
    required_group.add_argument('-corf', '--circ_orf', metavar='<FILE>', required=False,
        help='circRNA ORF file.\nRequired if -fi is not specified.')
    required_group.add_argument('-ms', '--path_ms_input', metavar='<PATH>', required=True,
        help='Folder containing Mass Spec result files.\nSupports MS-GF+, FragPipe, or summary tables.')
    optional_group.add_argument('-cm', '--circ_mapping', metavar='<FILE>', required=False,
        help='circRNA ORF mapping file.\nRequired if -fi is not specified.')
    optional_group.add_argument('-fip', '--file_input_protein', metavar='<FILE>', default='none',
        help='Linear protein sequence file (FASTA format).\n(default: none)')
    optional_group.add_argument('-fo', '--file_out', metavar='<PATH>', default='./',
        help='Output result folder path.\n(default: ./)')
    optional_group.add_argument('--deepcip_path', metavar='<PATH>', default='',
        help='Path to DeepCIP directory. Auto-detected if not specified.')
    optional_group.add_argument('--deepcircm6a_path', metavar='<PATH>', default='',
        help='Path to Deepcircm6a directory. Auto-detected if not specified.')
    optional_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
        help='Show this help message and exit.')

    try:
        args = parser.parse_args(argv)
    except SystemExit:
        raise

    if not args.circ_file_input and not (args.circ_seq and args.circ_orf):
        sys.exit("Error: circRNA sequence and ORF files are not provided via -fi or (-cs and -corf).")
    if not args.circ_file_input and not args.circ_mapping:
        print("[INFO] No circ_mapping file provided (-cm). ORF mapping/merging will be skipped; search ORF IDs used directly.")

    circ_file_input = args.circ_file_input
    path_ms_input = args.path_ms_input
    file_input_protein = args.file_input_protein
    file_out = args.file_out
    deepcip_path_arg = args.deepcip_path
    deepcircm6a_path_arg = args.deepcircm6a_path
    run_ires = 1

    if circ_file_input:
        if circ_file_input[-1] != '/':
            circ_file_input = circ_file_input + '/'
        circ_seq = circ_file_input + 'circRNA_full_length.fasta'
        circ_info_omit_file = circ_file_input + 'circ_info_omit.gff'
        circ_orf = circ_file_input + 'circRNA_bsj_corf.fasta'
        circ_mapping = circ_file_input + 'circRNA_ORF_mapping.tsv'
    else:
        circ_seq = args.circ_seq
        circ_info_omit_file = args.circ_info
        circ_orf = args.circ_orf
        circ_mapping = args.circ_mapping

    def ensure_output_dir(path):
        if path == './':
            return path
        if not path.endswith('/'):
            path += '/'
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        return path

    if file_out == './':
        now = time.strftime('%y-%m-%d', time.localtime())
        base_name = f'MScirctl2_{now}'
        ii = 1
        while os.path.exists(f"{base_name}.{ii}"):
            ii += 1
        file_out = f"{base_name}.{ii}/"
        os.makedirs(file_out, exist_ok=True)
    else:
        file_out = ensure_output_dir(file_out)

    # ---- Resolve external tool paths ----
    deepcip_path = deepcip_path_arg
    if not deepcip_path:
        _script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else '.'
        candidates = [os.path.join(_script_dir, 'DeepCIP'), os.path.expanduser('~/software/DeepCIP')]
        for c in candidates:
            if os.path.exists(c) and os.path.exists(os.path.join(c, 'DeepCIP.py')):
                deepcip_path = c
                break
    if not deepcip_path or not os.path.exists(os.path.join(deepcip_path, 'DeepCIP.py')):
        print("[WARNING] DeepCIP not found. IRES prediction will be skipped.")
        run_ires = 0

    deepcircm6a_path = deepcircm6a_path_arg
    if not deepcircm6a_path:
        _script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else '.'
        candidates = [os.path.join(_script_dir, 'Deepcircm6a'), os.path.expanduser('~/software/Deepcircm6a')]
        for c in candidates:
            if os.path.exists(c) and os.path.exists(os.path.join(c, 'predict.py')):
                deepcircm6a_path = c
                break

    # ---- Read mapping file ----
    dic_output_to_source = {}
    dic_source_aa = {}
    dic_source_orf_coords = {}
    if circ_mapping and os.path.exists(circ_mapping):
        with open(circ_mapping, 'r') as f:
            header = f.readline()
            for line in f:
                cols = line.strip('\n').split('\t')
                if len(cols) >= 2:
                    dic_output_to_source[cols[0]] = cols[1]
                    if len(cols) >= 3:
                        dic_source_aa[cols[1]] = cols[2]
                    m = re.search(r'\((-?\d+),(-?\d+)\)', cols[1])
                    if m:
                        dic_source_orf_coords[cols[0]] = (int(m.group(1)), int(m.group(2)))
                        dic_source_orf_coords[cols[1]] = (int(m.group(1)), int(m.group(2)))

    def circ_id_seq():
        dic_circ_id_seq = {}
        current_name = None
        with open(circ_seq, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    current_name = line[1:].split(':')[0]
                    dic_circ_id_seq[current_name] = []
                else:
                    if current_name is not None:
                        dic_circ_id_seq[current_name].append(line)
        for k in dic_circ_id_seq:
            dic_circ_id_seq[k] = "".join(dic_circ_id_seq[k])
        return dic_circ_id_seq

    def build_ires_input_seq(circ_orf_id):
        coords = dic_source_orf_coords.get(circ_orf_id)
        if coords:
            start_0based, end_0based = coords
        else:
            m = re.search(r'\((-?\d+),(-?\d+)\)', circ_orf_id)
            if not m:
                return []
            start_0based = int(m.group(1))
            end_0based = int(m.group(2))
        is_rct = 'RCT' in circ_orf_id
        circ_id = circ_orf_id.strip().split('-', 1)[0]
        full_seq = dic_circ_id_seq.get(circ_id, '')
        if not full_seq:
            return []
        s = full_seq.upper().replace('U', 'T')
        circ_len = len(s)
        if circ_len < 10:
            return []
        orf_length = end_0based - start_0based
        if orf_length < 0:
            orf_length = 0
        required_len = orf_length + 1000
        multiplier = max(10, int(math.ceil(required_len / circ_len)) + 4)
        s_multi = s * multiplier
        mid_copy_idx = multiplier // 2
        normalized_start = start_0based % circ_len
        anchor_start = normalized_start + mid_copy_idx * circ_len
        anchor_end = anchor_start + orf_length
        seqs = []
        s1_start = anchor_start - 137
        s1_end = anchor_start + 37
        if s1_start < 0 or s1_end > len(s_multi):
            seqs.append('N' * 174)
        else:
            seqs.append(s_multi[s1_start: s1_end])
        s2_start = anchor_start - 237
        s2_end = anchor_start - 63
        if s2_start < 0 or s2_end > len(s_multi):
            seqs.append('N' * 174)
        else:
            seqs.append(s_multi[s2_start: s2_end])
        if is_rct:
            seqs.append('')
        else:
            s3_start = anchor_end
            s3_end = anchor_end + 174
            if s3_start < 0 or s3_end > len(s_multi):
                seqs.append('')
            else:
                seqs.append(s_multi[s3_start: s3_end])
        return seqs

    def batch_ires_predict(orf_ids):
        dic_orf_seqs = {}
        print(f" -> Preparing sequences for IRES prediction...")
        skipped_count = 0
        for orf_id in orf_ids:
            seqs = build_ires_input_seq(orf_id)
            if not seqs:
                skipped_count += 1
                continue
            oid = orf_id.split()[0]
            dic_orf_seqs[oid] = seqs
        if skipped_count > 0:
            print(f" -> Note: Skipped {skipped_count} ORFs during IRES sequence extraction.")
        if not dic_orf_seqs:
            return {}
        temp_fasta = file_out + 'all_circRNA_ires.fasta'
        with open(temp_fasta, 'w+') as fp:
            for oid, seqs in dic_orf_seqs.items():
                if seqs[0] and len(seqs[0]) == 174:
                    fp.write(f'>{oid}_seq1\n{seqs[0]}\n')
                if seqs[1] and len(seqs[1]) == 174:
                    fp.write(f'>{oid}_seq2\n{seqs[1]}\n')
                if len(seqs) > 2 and seqs[2] and len(seqs[2]) == 174:
                    fp.write(f'>{oid}_seq3\n{seqs[2]}\n')
        if not deepcip_path or not os.path.exists(os.path.join(deepcip_path, 'DeepCIP.py')):
            print("[WARNING] DeepCIP path invalid, skipping IRES prediction.")
            return {}
        import shutil
        deepcip_data_dir = os.path.join(deepcip_path, 'data')
        if not os.path.exists(deepcip_data_dir):
            os.makedirs(deepcip_data_dir, exist_ok=True)
        deepcip_results_dir = os.path.join(deepcip_path, 'results')
        if not os.path.exists(deepcip_results_dir):
            os.makedirs(deepcip_results_dir, exist_ok=True)
        deepcip_input_name = f'temp_ires_{int(time.time())}.fasta'
        deepcip_input_path = os.path.join(deepcip_data_dir, deepcip_input_name)
        shutil.copy(temp_fasta, deepcip_input_path)
        deepcip_python = sys.executable
        predict_script_abs = os.path.join(deepcip_path, 'DeepCIP.py')
        dataset_name = f'run_{int(time.time())}'
        cmd_args = [deepcip_python, predict_script_abs, '-n', dataset_name, '-i', deepcip_input_name, '-m', '0']
        print(f" -> Running DeepCIP prediction using: {deepcip_python}")
        try:
            subprocess.run(cmd_args, cwd=deepcip_path, check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f" -> Error running DeepCIP (Exit code {e.returncode}): {e}")
        except Exception as e:
            print(f" -> Unexpected error running DeepCIP: {e}")
        input_base_name = os.path.splitext(deepcip_input_name)[0]
        deepcip_res_csv = os.path.join(deepcip_path, 'results', f'{input_base_name}_mode_0.csv')
        dic_ires_result = {}
        if os.path.exists(deepcip_res_csv):
            try:
                import pandas as pd
                df = pd.read_csv(deepcip_res_csv)
                temp_probs = {}
                rows_processed = 0
                for _, row in df.iterrows():
                    seq_id = str(row['Sequence_name']).strip()
                    try:
                        prob = float(row['Predict_probs'])
                    except (ValueError, KeyError):
                        continue
                    if '_seq' in seq_id:
                        parts_split = seq_id.rsplit('_seq', 1)
                        if len(parts_split) == 2:
                            oid = parts_split[0]
                            seq_idx = parts_split[1]
                            if oid not in temp_probs:
                                temp_probs[oid] = {}
                            temp_probs[oid][seq_idx] = prob
                            rows_processed += 1
                print(f" -> Parsed {rows_processed} rows from DeepCIP results")
                for oid, probs in temp_probs.items():
                    p1 = probs.get('1', -1.0)
                    p2 = probs.get('2', -1.0)
                    p3 = probs.get('3', -1.0)
                    dic_ires_result[oid] = {'p1': p1, 'p2': p2, 'p3': p3}
                print(f" -> Successfully mapped results for {len(dic_ires_result)} ORFs")
            except Exception as e:
                print(f" -> Error parsing DeepCIP results: {e}")
                import traceback
                traceback.print_exc()
        try:
            os.remove(deepcip_input_path)
            if os.path.exists(deepcip_res_csv):
                os.remove(deepcip_res_csv)
        except:
            pass
        return dic_ires_result

    def batch_m6a_predict(orf_entries):
        print(f" -> Preparing batch m6A prediction input sequences...")
        temp_fasta = file_out + 'all_circRNA_m6a_windows.fasta'
        valid_windows_count = 0
        dic_m6a_counts = {}
        with open(temp_fasta, 'w+') as fp:
            for orf_id, orf_seq in orf_entries:
                coords = dic_source_orf_coords.get(orf_id)
                if coords:
                    start_0based = coords[0]
                else:
                    m = re.search(r'\((-?\d+),-?\d+\)', orf_id)
                    if not m:
                        continue
                    start_0based = int(m.group(1))
                circ_id = orf_id.strip().split('-')[0]
                seq = dic_circ_id_seq.get(circ_id, '')
                if not seq:
                    continue
                s = seq.upper().replace('U', 'T')
                circ_len = len(s)
                if circ_len < 10:
                    continue
                s5 = s * 5
                anchor_pos = start_0based + 2 * circ_len
                region_start = anchor_pos - 126
                region_end = anchor_pos + 25
                if region_start < 0 or region_end > len(s5):
                    continue
                region = s5[region_start: region_end]
                n = len(region)
                for idx in range(25, n - 25):
                    if region[idx] == 'A':
                        window = region[idx - 25: idx + 26]
                        if len(window) != 51:
                            continue
                        rel_pos = idx - 126
                        header = f"{orf_id}###{rel_pos}"
                        fp.write(f">{header}\n{window}\n")
                        valid_windows_count += 1
        if valid_windows_count == 0:
            return {}
        print(f" -> Generated {valid_windows_count} prediction windows, calling Deepcircm6a...")
        if not deepcircm6a_path or not os.path.exists(os.path.join(deepcircm6a_path, 'predict.py')):
            print("[WARNING] Deepcircm6a path invalid, skipping m6A prediction.")
            return {}
        predict_py = os.path.join(deepcircm6a_path, 'predict.py')
        out_path = file_out + 'all_circRNA_m6a_pred.txt'
        cmd_args = [sys.executable, predict_py, '-predict_fa', temp_fasta, '-model_path', deepcircm6a_path, '-outfile', out_path]
        print(f" -> Running Deepcircm6a prediction...")
        try:
            subprocess.run(cmd_args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f" -> Error running Deepcircm6a (Exit code {e.returncode})")
        except Exception as e:
            print(f" -> Unexpected error running Deepcircm6a: {e}")
        if os.path.exists(out_path):
            with open(out_path, 'r') as fr:
                for line in fr:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    try:
                        header = parts[0]
                        prob = float(parts[-1])
                        if prob >= 0.5:
                            real_orf_id = header.split('###')[0]
                            dic_m6a_counts[real_orf_id] = dic_m6a_counts.get(real_orf_id, 0) + 1
                    except Exception:
                        continue
        return dic_m6a_counts

    def dic_circ_mgf_pep_add(circ_id, mgf_id, peptide, dic_id_peptide, count=1):
        val = f"{mgf_id}###{peptide}###{int(count)}"
        if circ_id in dic_id_peptide:
            dic_id_peptide[circ_id] += ';' + val
        else:
            dic_id_peptide[circ_id] = val
        return dic_id_peptide

    def detect_file_type(file_path):
        try:
            with open_safe_read(file_path) as f:
                header = f.readline()
                if not header:
                    return 'unknown'
                if 'MS_ID' in header and 'peptide' in header and 'cORF_ID' in header:
                    return 'summary'
                parts = [h.strip() for h in header.split('\t')]
                if 'Stripped.Sequence' in parts and ('Protein.Group' in parts or 'Protein.Ids' in parts):
                    return 'fragpipe_combined'
                
                if ('Peptide Sequence' in parts or 'Peptide' in parts or 'Sequence' in parts) and \
                   ('Protein' in parts or 'Protein ID' in parts):
                    if any('Spectral Count' in h for h in parts):
                        return 'fragpipe'
                    if 'SpecEValue' in parts or 'EValue' in parts:
                        return 'msgf'
                    return 'fragpipe'
                
                if '#SpecFile' in header or 'SpecEValue' in parts or 'EValue' in parts:
                    return 'msgf'
                return 'unknown'
        except:
            return 'unknown'

    def parse_msgf_file(file_path, file_line_circ, dic_id_peptide, global_stats):
        content = open_safe_read(file_path)
        header_line = content.readline()
        if not header_line:
            return dic_id_peptide
        headers = header_line.strip().split('\t')
        try:
            idx_evalue = headers.index('EValue')
        except ValueError:
            try:
                idx_evalue = headers.index('SpecEValue')
            except ValueError:
                idx_evalue = 13
        try:
            idx_file = headers.index('#SpecFile')
        except ValueError:
            try:
                idx_file = headers.index('SpecFile')
            except ValueError:
                idx_file = 0
        try:
            idx_peptide = headers.index('Peptide')
        except ValueError:
            idx_peptide = 9
        try:
            idx_protein = headers.index('Protein')
        except ValueError:
            idx_protein = 10
        for line in content:
            columns = line.strip().split('\t')
            if len(columns) <= max(idx_evalue, idx_file, idx_peptide, idx_protein):
                continue
            try:
                EValue = float(columns[idx_evalue].strip())
            except ValueError:
                continue
            if EValue <= 0.01:
                file_id = columns[idx_file].replace('.mzML', '.mgf')
                global_stats['psms'] += 1
                global_stats['files'].add(file_id)
                raw_pep = columns[idx_peptide].strip()
                Peptide = clean_peptide_sequence(raw_pep)
                if len(Peptide) < 6:
                    continue
                protein_col = columns[idx_protein].strip()
                circ_id_multiple = protein_col.split(', ')
                for circ_id in circ_id_multiple:
                    clean_circ_id = circ_id
                    line_out = f"{file_id}\t{Peptide}\t{clean_circ_id}\n"
                    file_line_circ.write(line_out)
                    if file_input_protein == 'none':
                        dic_id_peptide = dic_circ_mgf_pep_add(clean_circ_id, file_id, Peptide, dic_id_peptide, 1)
        return dic_id_peptide

    def parse_fragpipe_file(file_path, file_line_circ, dic_id_peptide, global_stats):
        content = open_safe_read(file_path)
        header_line = content.readline()
        headers = header_line.strip().split('\t')
        
        try:
            if 'Peptide Sequence' in headers:
                idx_peptide = headers.index('Peptide Sequence')
            elif 'Peptide' in headers:
                idx_peptide = headers.index('Peptide')
            elif 'Sequence' in headers:
                idx_peptide = headers.index('Sequence')
            else:
                return dic_id_peptide
            
            if 'Protein' in headers:
                idx_protein = headers.index('Protein')
            elif 'Protein ID' in headers:
                idx_protein = headers.index('Protein ID')
            else:
                return dic_id_peptide
            
            sample_cols = []
            for i, h in enumerate(headers):
                if h.endswith(' Spectral Count'):
                    sample_name = h.replace(' Spectral Count', '').strip()
                    sample_cols.append({'idx': i, 'name': sample_name, 'type': 'count'})
            
            if not sample_cols and 'Spectral Count' in headers:
                idx_spec = headers.index('Spectral Count')
                base_name = os.path.splitext(os.path.basename(file_path))[0]
                sample_name = base_name.split('peptide_', 1)[1] if base_name.startswith('peptide_') else base_name
                sample_cols = [{'idx': idx_spec, 'name': sample_name, 'type': 'count'}]
            
            if not sample_cols:
                for i, h in enumerate(headers):
                    if h.endswith(' MaxLFQ Intensity') or h.endswith(' Intensity'):
                        sample_name = h.replace(' MaxLFQ Intensity', '').replace(' Intensity', '').strip()
                        if not any(s['name'] == sample_name for s in sample_cols):
                            sample_cols.append({'idx': i, 'name': sample_name, 'type': 'intensity'})
            
            if not sample_cols:
                return dic_id_peptide
        except ValueError:
            return dic_id_peptide

        for line in content:
            line = line.strip()
            if not line:
                continue
            columns = line.split('\t')
            max_idx = max([s['idx'] for s in sample_cols])
            if len(columns) <= max_idx:
                continue
            
            raw_pep = columns[idx_peptide].strip().strip('"')
            peptide_seq = clean_peptide_sequence(raw_pep)
            if len(peptide_seq) < 6:
                continue
            
            protein_val = columns[idx_protein].strip().strip('"')
            clean_circ_id = protein_val
            
            for sample in sample_cols:
                try:
                    val_str = columns[sample['idx']].strip()
                    if sample['type'] == 'intensity':
                        count_val = 1 if float(val_str) > 0 else 0
                    else:
                        count_val = int(float(val_str)) if val_str else 0
                except ValueError:
                    count_val = 0
                
                if count_val > 0:
                    file_id = sample['name']
                    global_stats['files'].add(file_id)
                    for _ in range(count_val):
                        global_stats['psms'] += 1
                        line_out = f"{file_id}\t{peptide_seq}\t{clean_circ_id}\n"
                        file_line_circ.write(line_out)
                    if file_input_protein == 'none':
                        dic_id_peptide = dic_circ_mgf_pep_add(clean_circ_id, file_id, peptide_seq, dic_id_peptide, count_val)
        return dic_id_peptide

    def parse_fragpipe_combined_file(file_path, file_line_circ, dic_id_peptide, global_stats):
        content = open_safe_read(file_path)
        header_line = content.readline()
        headers = header_line.strip().split('\t')
        try:
            idx_peptide = headers.index('Stripped.Sequence')
        except ValueError:
            return dic_id_peptide
        idx_protein = -1
        for col_name in ['Protein.Group', 'Protein.Ids', 'Protein', 'Protein ID']:
            if col_name in headers:
                idx_protein = headers.index(col_name)
                break
        if idx_protein == -1:
            return dic_id_peptide
        metadata_cols = {
            'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes',
            'First.Protein.Description', 'Proteotypic', 'Stripped.Sequence',
            'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id',
            'All Mapped Proteins', 'All Mapped Genes', 'Gene',
            'Protein Description', 'Protein ID', 'Entry Name',
            'Mapped Genes', 'Mapped Proteins'
        }
        sample_cols = []
        for i, h in enumerate(headers):
            h_stripped = h.strip()
            if h_stripped in metadata_cols or h_stripped.startswith('#'):
                continue
            if i > idx_peptide:
                sample_cols.append({'idx': i, 'name': h_stripped})
        for line in content:
            line = line.strip()
            if not line:
                continue
            columns = line.split('\t')
            if len(columns) <= max(idx_peptide, idx_protein):
                continue
            peptide_seq = columns[idx_peptide].strip().strip('"')
            peptide_seq = clean_peptide_sequence(peptide_seq)
            if len(peptide_seq) < 6:
                continue
            protein_val = columns[idx_protein].strip().strip('"')
            clean_circ_id = protein_val
            for sample in sample_cols:
                if sample['idx'] >= len(columns):
                    continue
                try:
                    val = float(columns[sample['idx']].strip())
                except (ValueError, IndexError):
                    continue
                if val > 0:
                    file_id = sample['name']
                    if '/' in file_id:
                        file_id = os.path.splitext(os.path.basename(file_id))[0]
                    global_stats['psms'] += 1
                    global_stats['files'].add(file_id)
                    line_out = f"{file_id}\t{peptide_seq}\t{clean_circ_id}\n"
                    file_line_circ.write(line_out)
                    if file_input_protein == 'none':
                        dic_id_peptide = dic_circ_mgf_pep_add(clean_circ_id, file_id, peptide_seq, dic_id_peptide, 1)
                    break
        return dic_id_peptide

    def parse_summary_file(file_path, file_line_circ, dic_id_peptide, global_stats):
        content = open_safe_read(file_path)
        header = content.readline()
        for line in content:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            file_id = parts[0]
            peptide = clean_peptide_sequence(parts[1])
            if len(peptide) < 6:
                continue
            circ_id = parts[2]
            global_stats['psms'] += 1
            global_stats['files'].add(file_id)
            line_out = f"{file_id}\t{peptide}\t{circ_id}\n"
            file_line_circ.write(line_out)
            if file_input_protein == 'none':
                dic_id_peptide = dic_circ_mgf_pep_add(circ_id, file_id, peptide, dic_id_peptide, 1)
        return dic_id_peptide

    def process_ms_directory():
        dic_id_peptide = {}
        outfile_path = file_out + "tsv_result.txt"
        file_line_circ = open(outfile_path, "w+")
        file_line_circ.write('MS_ID\tpeptide\tcORF_ID\n')
        global_stats = {'psms': 0, 'files': set()}
        if os.path.isfile(path_ms_input):
            files_to_process = [path_ms_input]
        else:
            files_to_process = []
            for root, dirs, files in os.walk(path_ms_input):
                for f in files:
                    if f.lower().endswith(('.tsv', '.txt')):
                        files_to_process.append(os.path.join(root, f))
        print(f"Scanning {len(files_to_process)} files in {path_ms_input}...")
        success_count = 0
        for file_path in files_to_process:
            ftype = detect_file_type(file_path)
            if ftype == 'msgf':
                dic_id_peptide = parse_msgf_file(file_path, file_line_circ, dic_id_peptide, global_stats)
                success_count += 1
            elif ftype == 'fragpipe':
                dic_id_peptide = parse_fragpipe_file(file_path, file_line_circ, dic_id_peptide, global_stats)
                success_count += 1
            elif ftype == 'fragpipe_combined':
                dic_id_peptide = parse_fragpipe_combined_file(file_path, file_line_circ, dic_id_peptide, global_stats)
                success_count += 1
            elif ftype == 'summary':
                dic_id_peptide = parse_summary_file(file_path, file_line_circ, dic_id_peptide, global_stats)
                success_count += 1
        file_line_circ.close()
        TOTAL_INPUT_PSMS = global_stats['psms']
        TOTAL_MS_FILES = len(global_stats['files'])
        print(f"Successfully processed {success_count} files.")
        print(f"Total input PSMs: {TOTAL_INPUT_PSMS}")
        print(f"Total MS files/samples: {TOTAL_MS_FILES}")
        return dic_id_peptide

    def normalize_prot_id(pid):
        p = pid.strip()
        if '|' in p:
            p = p.split('|')[-1]
        p = p.split()[0]
        if '.' in p:
            p = p.split('.')[0]
        return p.upper()

    def normalize_circ_key(cid):
        return cid.strip()

    KMER_SIZE_LINEAR = 6

    def tsv_remove_line_protaid():
        global dic_removed_line_peptides
        dic_removed_line_peptides = {}
        file_line_protide = open(file_input_protein, 'r')
        current_id = None
        current_seq = ''
        list_line_entries = []
        for line_protide in file_line_protide:
            if line_protide[0] == '>':
                if current_id is not None:
                    list_line_entries.append((current_id, current_seq, current_seq.upper().replace('I', 'L')))
                current_id = line_protide[1:].strip().split()[0]
                current_seq = ''
            else:
                current_seq += line_protide.strip()
        if current_id is not None:
            list_line_entries.append((current_id, current_seq, current_seq.upper().replace('I', 'L')))
        file_line_protide.close()
        kmer_index = {}
        for i, (prot_id, prot_seq, prot_seq_IL) in enumerate(list_line_entries):
            if len(prot_seq_IL) < KMER_SIZE_LINEAR:
                continue
            seen = set()
            for j in range(len(prot_seq_IL) - KMER_SIZE_LINEAR + 1):
                kmer = prot_seq_IL[j:j + KMER_SIZE_LINEAR]
                if kmer not in seen:
                    seen.add(kmer)
                    kmer_index.setdefault(kmer, set()).add(i)
        file_tsv_result = None
        try:
            file_tsv_result = open(file_out + "tsv_result.txt", 'r')
        except FileNotFoundError:
            print(f"[WARNING] {file_out}tsv_result.txt does not exist")
        if file_tsv_result is None:
            return {}
        file_tsv_result_header = file_tsv_result.readline()
        file_remove = open(file_out + 'tsv_result_remove_line.txt', 'w+')
        file_remove.write(file_tsv_result_header)
        file_peptide_line = open(file_out + 'peptide_line.txt', 'w+')
        file_peptide_line.write('peptide\tlinear_protein\tMS_ID\n')
        set_line_p = set()
        set_circ_p = set()
        dic_id_peptide = {}
        for line_tsv_result in file_tsv_result:
            parts = line_tsv_result.strip().split('\t')
            if len(parts) < 3:
                continue
            mgf_id = parts[0]
            peptide = parts[1]
            circ_id = parts[2]
            peptide_IL = peptide.upper().replace('I', 'L')
            if peptide_IL in set_line_p:
                continue
            if peptide_IL in set_circ_p:
                file_remove.write(line_tsv_result)
                dic_id_peptide = dic_circ_mgf_pep_add(circ_id, mgf_id, peptide, dic_id_peptide, 1)
                continue
            matched_prot_id = None
            matched_prot_norm = None
            if len(peptide_IL) >= KMER_SIZE_LINEAR:
                candidate_indices = None
                for j in range(len(peptide_IL) - KMER_SIZE_LINEAR + 1):
                    kmer = peptide_IL[j:j + KMER_SIZE_LINEAR]
                    indices = kmer_index.get(kmer)
                    if indices is None:
                        candidate_indices = set()
                        break
                    if candidate_indices is None or len(indices) < len(candidate_indices):
                        candidate_indices = indices
                    if len(candidate_indices) == 1:
                        break
                if candidate_indices:
                    for idx in candidate_indices:
                        prot_id, prot_seq, prot_seq_IL = list_line_entries[idx]
                        if peptide_IL in prot_seq_IL:
                            matched_prot_id = prot_id
                            matched_prot_norm = normalize_prot_id(prot_id)
                            break
                else:
                    for prot_id, prot_seq, prot_seq_IL in list_line_entries:
                        if peptide_IL in prot_seq_IL:
                            matched_prot_id = prot_id
                            matched_prot_norm = normalize_prot_id(prot_id)
                            break
            else:
                for prot_id, prot_seq, prot_seq_IL in list_line_entries:
                    if peptide_IL in prot_seq_IL:
                        matched_prot_id = prot_id
                        matched_prot_norm = normalize_prot_id(prot_id)
                        break
            if matched_prot_id:
                set_line_p.add(peptide_IL)
                file_peptide_line.write(f"{peptide}\t{matched_prot_norm}\t{mgf_id}\n")
                circ_key_full = normalize_circ_key(circ_id)
                lst = dic_removed_line_peptides.get(circ_key_full, [])
                lst.append((peptide, mgf_id, matched_prot_norm))
                dic_removed_line_peptides[circ_key_full] = lst
            else:
                set_circ_p.add(peptide_IL)
                file_remove.write(line_tsv_result)
                dic_id_peptide = dic_circ_mgf_pep_add(circ_id, mgf_id, peptide, dic_id_peptide, 1)
        TOTAL_LINEAR_PEPS = len(set_line_p)
        TOTAL_CIRC_PEPS = len(set_circ_p)
        print(f"Total linear RNA peptides: {TOTAL_LINEAR_PEPS}")
        print(f"Total circular peptides: {TOTAL_CIRC_PEPS}")
        file_remove.close()
        file_peptide_line.close()
        file_tsv_result.close()
        return dic_id_peptide

    def batch_diamond_alignment_filter(peptide_list):
        if not os.path.exists(file_input_protein):
            return set()
        if not peptide_list:
            return set()
        db_name = f"{file_out}_diamond_protein_db"
        if not os.path.exists(db_name + '.dmnd'):
            print("Building Diamond database...")
            subprocess.run(['diamond', 'makedb', '--in', file_input_protein, '-d', db_name, '--quiet'], check=False)
        temp_query = file_out + 'diamond_peptides_query.fasta'
        with open(temp_query, 'w+') as f:
            for i, pep in enumerate(peptide_list):
                f.write(f'>pep_{i}\n{pep}\n')
        temp_result = file_out + 'diamond_filter_result.txt'
        cmd_args = ['diamond', 'blastp', '-q', temp_query, '-d', db_name, '-o', temp_result,
                    '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen',
                    '-p', '8', '-e', '0.001', '--id', '95', '--quiet']
        print(f"Running peptide homology scan (Diamond, {len(peptide_list)} peptides)...")
        try:
            subprocess.run(cmd_args, check=False)
        except Exception as e:
            print(f"[WARNING] Diamond execution error: {e}")
            return set()
        blacklisted_indices = set()
        if os.path.exists(temp_result):
            with open(temp_result, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 5:
                        continue
                    qseqid = parts[0]
                    try:
                        pident = float(parts[2])
                        aln_len = int(parts[3])
                        q_len = int(parts[4])
                    except (ValueError, IndexError):
                        continue
                    if q_len == 0:
                        continue
                    effective_coverage = (aln_len * (pident / 100.0)) / q_len
                    if effective_coverage < 0.95:
                        continue
                    try:
                        idx = int(qseqid.split('_')[1])
                        blacklisted_indices.add(idx)
                    except:
                        pass
        blacklisted_peptides = set()
        for idx in blacklisted_indices:
            if idx < len(peptide_list):
                blacklisted_peptides.add(peptide_list[idx])
        print(f"Homology scan complete, marked {len(blacklisted_peptides)} peptides highly similar (>=95%) to linear proteins.")
        return blacklisted_peptides

    def parse_bsj_location(bsj_location_str):
        if not bsj_location_str or str(bsj_location_str).strip() == '0':
            return []
        res = []
        for x in str(bsj_location_str).split('|'):
            x = x.strip().split()[0] 
            if x:
                try:
                    res.append(float(x))
                except ValueError:
                    pass
        return res

    def is_peptide_crossing_bsj(pep_start_0based, pep_end_0based, bsj_list):
        if not bsj_list:
            return False
        for b_val in bsj_list:
            b_int = int(b_val)
            b_dec = round(b_val - b_int, 1)
            b_0based = b_int - 1
            if b_dec == 0.3:
                if pep_start_0based <= b_0based and pep_end_0based >= b_0based + 2:
                    return True
            else:
                if pep_start_0based <= b_0based < pep_end_0based:
                    return True
        return False

    def circ_forward_reverse():
        dic_circ_strand = {}
        if not circ_info_omit_file or not os.path.exists(circ_info_omit_file):
            return dic_circ_strand
        with open(circ_info_omit_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    dic_circ_strand[parts[0]] = parts[5]
        return dic_circ_strand


    def build_source_bsj_location():
        dic_bsj = {}

        # ---- 从已加载的circRNA序列字典构建长度字典 ----
        dic_circ_len = {}
        for circ_id, seq in dic_circ_id_seq.items():
            dic_circ_len[circ_id] = len(seq)
        print(f" [BSJ] Loaded {len(dic_circ_len)} circRNA lengths for BSJ computation")

        # ---- 根据ORF坐标和circRNA长度计算BSJ的辅助函数 ----
        def compute_bsj_from_coords(header_str):
            try:
                m_coord = re.search(r'\((-?\d+),(-?\d+)\)', header_str)
                if not m_coord:
                    return None

                orf_start = int(m_coord.group(1))
                orf_end = int(m_coord.group(2))
                orf_nt_len = orf_end - orf_start

                if orf_nt_len <= 0:
                    return '0'

                # 提取circRNA ID
                base_id = header_str.split('(')[0]
                circ_id_local = base_id.split('-ORF')[0] if '-ORF' in base_id else base_id
                circ_id_local = extract_circ_info_id(circ_id_local)

                circ_len = dic_circ_len.get(circ_id_local, 0)
                if circ_len <= 0:
                    return None  # 找不到对应circRNA长度

                # ORF未跨越BSJ
                if orf_end <= circ_len:
                    return '0'

                # ORF跨越BSJ - 按A脚本RNA_protein()的算法计算断点位置
                breakpoints = []
                j1_dist = circ_len - orf_start
                current_j_dist = j1_dist
                while current_j_dist < orf_nt_len:
                    aa_idx = (current_j_dist // 3) + 1
                    rem = current_j_dist % 3
                    if aa_idx <= 0:
                        aa_idx = 1
                    bp_str = f"{aa_idx - 1}.3" if rem == 0 else f"{aa_idx}.{rem}"
                    breakpoints.append(bp_str)
                    current_j_dist += circ_len

                return '|'.join(breakpoints) if breakpoints else '0'
            except Exception as e:
                print(f" [BSJ] Error computing BSJ for '{header_str}': {e}")
                return '0'

        # ---- 解析circ_orf文件 ----
        orf_count = 0
        computed_count = 0
        if os.path.exists(circ_orf):
            with open(circ_orf, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        line = line.strip()
                        header = line[1:].split()[0]
                        base_id = header.split('(')[0]
                        parts = header.split(':')
                        orf_count += 1

                        # 优先从头中提取合法的BSJ信息(含小数点的N.X格式)
                        bsj_found = False
                        if len(parts) >= 3:
                            for c in [parts[2], parts[-1]]:
                                if c and re.match(r'^[\d.|]+$', c) and '.' in c:
                                    dic_bsj[base_id] = c
                                    bsj_found = True
                                    break
                        elif len(parts) == 2:
                            c = parts[1]
                            if c and re.match(r'^[\d.|]+$', c) and '.' in c:
                                dic_bsj[base_id] = c
                                bsj_found = True

                        # 头中无合法BSJ -> 从坐标计算
                        if not bsj_found:
                            computed = compute_bsj_from_coords(header)
                            if computed is not None:
                                dic_bsj[base_id] = computed
                                computed_count += 1
                            else:
                                dic_bsj[base_id] = '0'

        print(f" [BSJ] Parsed {orf_count} ORFs, computed BSJ for {computed_count}, total BSJ entries: {len(dic_bsj)}")

        # ---- 解析mapping文件 ----
        if circ_mapping and os.path.exists(circ_mapping):
            with open(circ_mapping, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.strip('\n').split('\t')
                    if len(cols) < 2:
                        continue
                    source_header = cols[1]
                    source_base = source_header.split('(')[0]
                    parts = source_header.split(':')
                    src_bsj = None

                    if len(parts) >= 3:
                        for c in [parts[-1], parts[2]]:
                            if c and re.match(r'^[\d.|]+$', c) and '.' in c:
                                src_bsj = c
                                break
                    elif len(parts) == 2:
                        c = parts[1]
                        if c and re.match(r'^[\d.|]+$', c) and '.' in c:
                            src_bsj = c

                    if src_bsj is None:
                        computed = compute_bsj_from_coords(source_header)
                        if computed is not None:
                            src_bsj = computed

                    if src_bsj is not None:
                        dic_bsj[source_base] = src_bsj
                    elif source_base not in dic_bsj:
                        out_base = cols[0].split('(')[0]
                        dic_bsj[source_base] = dic_bsj.get(out_base, '0')

        return dic_bsj





    def map_to_source_orfs():
        """Map output ORFs to source ORFs and merge peptides."""
        dic_source_peptide_counts = {}
        dic_source_search_orfs = {}
        for output_orf_id, pep_string in dic_id_peptide.items():
            source_header = dic_output_to_source.get(output_orf_id, output_orf_id)
            dic_source_search_orfs.setdefault(source_header, []).append(output_orf_id)
            source_counts = dic_source_peptide_counts.setdefault(source_header, {})
            for raw_entry in pep_string.strip(';').split(';'):
                if not raw_entry.strip():
                    continue
                parts = raw_entry.strip().split('###')
                if len(parts) < 2:
                    continue
                ms_id = parts[0].strip()
                peptide = parts[1].strip().upper()
                try:
                    count = max(1, int(float(parts[2]))) if len(parts) >= 3 else 1
                except (TypeError, ValueError):
                    count = 1
                key = (ms_id, peptide)
                source_counts[key] = source_counts.get(key, 0) + count
        dic_source_peptides = {}
        for k, v in dic_source_peptide_counts.items():
            dic_source_peptides[k] = ';'.join(
                f"{ms_id}###{peptide}###{count}"
                for (ms_id, peptide), count in sorted(v.items(), key=lambda item: (item[0][1], item[0][0]))
            )
        return dic_source_peptides, dic_source_search_orfs

    def circ_orf_tsv_peptide_map(circ_orf_id, corf):
        if circ_orf_id not in dic_id_peptide:
            return 'none', 'none', 'none', 'none', 'none', 'none', 'none', {}, 0, ''
        mgf_pep_score_assemble = dic_id_peptide[circ_orf_id]
        dic_pro_sore = {}
        dic_pro_ms_counts = {}
        dic_mgf_score = {}
        list_mgf_id = []
        stats_cross = {'pep_count': 0, 'unique_pep': set(), 'ms_ids': set()}
        stats_non = {'pep_count': 0, 'unique_pep': set(), 'ms_ids': set()}
        corf_clean = corf.strip('*').upper()
        corf_clean_IL = corf_clean.replace('I', 'L')
        corf_s = len(corf_clean)
        coverage_mask = [False] * corf_s
        coverage_map = ['x'] * corf_s
        bsj_location_str = dic_circ_orf_bsj_location.get(circ_orf_id.split('(')[0], '0')
        bsj_list = parse_bsj_location(bsj_location_str)
        entries = [x.strip().split('###') for x in mgf_pep_score_assemble.strip(';').split(';') if x.strip()]
        for parts in entries:
            if len(parts) < 2:
                continue
            mgf_id = parts[0]
            pro = parts[1].upper()
            pro_IL = pro.replace('I', 'L')
            try:
                score = max(1, int(float(parts[2]))) if len(parts) >= 3 else 1
            except (TypeError, ValueError):
                score = 1
            list_mgf_id.append(mgf_id)
            dic_mgf_score[mgf_id] = dic_mgf_score.get(mgf_id, 0) + score
            if pro not in dic_pro_sore:
                dic_pro_sore[pro] = {}
            dic_pro_sore[pro][mgf_id] = dic_pro_sore[pro].get(mgf_id, 0) + score
            dic_pro_ms_counts.setdefault(pro, set()).add(mgf_id)
            start = 0
            is_crossing = False
            while True:
                idx = corf_clean_IL.find(pro_IL, start)
                if idx == -1:
                    break
                end = idx + len(pro)
                for k in range(idx, min(end, corf_s)):
                    coverage_mask[k] = True
                current_instance_crossing = is_peptide_crossing_bsj(idx, end, bsj_list)
                if current_instance_crossing:
                    is_crossing = True
                for i in range(idx, min(end, corf_s)):
                    char = corf[i]
                    if current_instance_crossing:
                        coverage_map[i] = char.upper()
                    else:
                        if coverage_map[i] == 'x' or coverage_map[i].islower():
                            coverage_map[i] = char.lower()
                start = idx + 1
            if is_crossing:
                stats_cross['pep_count'] += score
                stats_cross['unique_pep'].add(pro)
                stats_cross['ms_ids'].add(mgf_id)
            else:
                stats_non['pep_count'] += score
                stats_non['unique_pep'].add(pro)
                stats_non['ms_ids'].add(mgf_id)
        pro_number_score_parts = []
        number_all = stats_cross['pep_count'] + stats_non['pep_count']
        unique_peptide_seqs_list = []
        all_unique_peps = set(dic_pro_sore.keys())
        unique_peptide_count = len(all_unique_peps)
        for pro in all_unique_peps:
            total_count = sum(dic_pro_sore[pro].values())
            ms_hit_count = len(dic_pro_ms_counts.get(pro, set()))
            avg_score = (total_count / ms_hit_count) if ms_hit_count > 0 else 0.0
            pro_number_score_parts.append(f"{pro}|{total_count}|{avg_score:.3f}")
            if pro in stats_cross['unique_pep']:
                seq_str = pro.upper()
            else:
                seq_str = pro.lower()
            unique_peptide_seqs_list.append(f"{seq_str}({ms_hit_count}|{total_count})")
        pro_number_score = ";".join(pro_number_score_parts) + ";"
        peptide_seq_str = ";".join(unique_peptide_seqs_list)
        mgf_score_all = "".join([f"{k}:{v};" for k, v in dic_mgf_score.items()])
        raw_sum = len(set(list_mgf_id))
        covered_len = sum(coverage_mask)
        corf_peptide_vatio = covered_len / corf_s if corf_s > 0 else 0.0
        marge_peptide_multiple = "".join(coverage_map)
        if '*' in corf:
            marge_peptide_multiple += '*'
        marge_peptide_bsj_multiple = marge_peptide_multiple
        if bsj_list:
            seen_floats = set()
            unique_bsj = []
            for b in bsj_list:
                if b not in seen_floats:
                    seen_floats.add(b)
                    unique_bsj.append(b)
            from collections import defaultdict as _defaultdict
            pos_groups = _defaultdict(list)
            for b in unique_bsj:
                pos_groups[int(b)].append(b)
            sorted_positions = sorted(pos_groups.keys(), reverse=True)
            temp_list = list(marge_peptide_multiple)
            for pos in sorted_positions:
                if 0 < pos < len(temp_list):
                    count = len(pos_groups[pos])
                    for _ in range(count):
                        temp_list.insert(pos, '|')
            marge_peptide_bsj_multiple = "".join(temp_list)
        detailed_stats = {
            'cross': {'pep_count': stats_cross['pep_count'], 'unique_pep': len(stats_cross['unique_pep']),
                      'ms_count': len(stats_cross['ms_ids'])},
            'non': {'pep_count': stats_non['pep_count'], 'unique_pep': len(stats_non['unique_pep']),
                    'ms_count': len(stats_non['ms_ids'])}
        }
        return pro_number_score, marge_peptide_bsj_multiple, corf_peptide_vatio, bsj_location_str, raw_sum, number_all, mgf_score_all, detailed_stats, unique_peptide_count, peptide_seq_str

    def calculate_and_write_scores(results_buffer, file_path_out):
        def get_stat_list(key1, key2):
            return [d['stats']['detailed_stats'][key1][key2] for d in results_buffer]

        Cross_Pep_Vals = get_stat_list('cross', 'pep_count')
        Cross_Uni_Vals = get_stat_list('cross', 'unique_pep')
        Cross_MS_Vals = get_stat_list('cross', 'ms_count')
        Non_Pep_Vals = get_stat_list('non', 'pep_count')
        Non_Uni_Vals = get_stat_list('non', 'unique_pep')
        Non_MS_Vals = get_stat_list('non', 'ms_count')
        C_P_max = max(Cross_Pep_Vals) if Cross_Pep_Vals else 0
        C_U_max = max(Cross_Uni_Vals) if Cross_Uni_Vals else 0
        C_M_max = max(Cross_MS_Vals) if Cross_MS_Vals else 0
        N_P_max = max(Non_Pep_Vals) if Non_Pep_Vals else 0
        N_U_max = max(Non_Uni_Vals) if Non_Uni_Vals else 0
        N_M_max = max(Non_MS_Vals) if Non_MS_Vals else 0
        all_orf_lengths = [d['stats']['L'] for d in results_buffer]
        median_orf_len = float(np.median(all_orf_lengths)) if all_orf_lengths else 1.0

        def get_nonzero_percentile(values, q=10, default=1.0):
            nonzero = [v for v in values if v > 0]
            if not nonzero:
                return default
            return float(np.percentile(nonzero, q))

        K_C_P = get_nonzero_percentile(Cross_Pep_Vals, q=10)
        K_C_U = get_nonzero_percentile(Cross_Uni_Vals, q=10)
        K_C_M = get_nonzero_percentile(Cross_MS_Vals, q=10)
        K_N_P = get_nonzero_percentile(Non_Pep_Vals, q=10)
        K_N_U = get_nonzero_percentile(Non_Uni_Vals, q=10)
        K_N_M = get_nonzero_percentile(Non_MS_Vals, q=10)

        def mixed_score(val, k, max_val):
            if val <= 0:
                return 0.0
            if max_val <= 0:
                return 0.0
            hill = val / (val + k) if k > 0 else 0.0
            log_ratio = math.log(1 + val) / math.log(1 + max_val)
            return 0.7 * hill + 0.3 * log_ratio

        has_info = bool(circ_info_omit_file and os.path.exists(circ_info_omit_file))
        header_cols = ['circ_ORF_ID', 'search_ORF_IDs', 'circRNA_ID', 'gene']
        if has_info:
            header_cols.append('strand')
        header_cols.extend(['IRES_Up50_Center', 'IRES_Up150_Center', 'IRES_Down_Stop',
                            'm6A_Count', 'cORF_seq', 'cORF_len', 'cORF_map',
                            'peptide_coverage', 'cORF_bsj_position', 'MS_count',
                            'total_peptide_count', 'unique_peptide_count', 'peptide_seq',
                            'ms_sample_evidence', 'coding_function_score'])
        scored_rows = []
        for item in results_buffer:
            stats = item['stats']
            ds = stats['detailed_stats']
            IRES_SCORE = stats['IRES_SCORE']
            m = stats['m']
            cov = stats['cov']
            orf_len = stats['L']
            try:
                score_ires = 0.15 * float(IRES_SCORE)
            except (ValueError, TypeError):
                score_ires = 0.0
            score_m6a = 0.15 * (1 - math.exp(-0.15 * m))
            len_adjustment = (orf_len / median_orf_len) ** 0.1 if median_orf_len > 0 else 1.0
            score_cov = 0.10 * cov * len_adjustment
            s_c_p = mixed_score(ds['cross']['pep_count'], K_C_P, C_P_max)
            s_c_u = mixed_score(ds['cross']['unique_pep'], K_C_U, C_U_max)
            s_c_m = mixed_score(ds['cross']['ms_count'], K_C_M, C_M_max)
            score_cross = 0.2 * s_c_p + 0.4 * s_c_u + 0.4 * s_c_m
            s_n_p = mixed_score(ds['non']['pep_count'], K_N_P, N_P_max)
            s_n_u = mixed_score(ds['non']['unique_pep'], K_N_U, N_U_max)
            s_n_m = mixed_score(ds['non']['ms_count'], K_N_M, N_M_max)
            score_non = 0.2 * s_n_p + 0.4 * s_n_u + 0.4 * s_n_m
            P = score_ires + score_m6a + score_cov + 0.4 * score_cross + 0.2 * score_non
            parts = item['data_row'].strip().split('\t')
            # Remove strand column if no info file (strand is at index 4)
            if not has_info:
                parts.pop(4)
            n_base = len(header_cols) - 1
            base_parts = parts[:n_base]
            base_row = '\t'.join(base_parts)
            scored_rows.append((P, f"{base_row}\t{P:.4f}\n"))
        scored_rows.sort(key=lambda x: x[0], reverse=True)
        with open(file_path_out, 'w+') as file_circ_predict:
            header = '\t'.join(header_cols) + '\n'
            file_circ_predict.write(header)
            for _, row_text in scored_rows:
                file_circ_predict.write(row_text)

    def circ_peptide_result():
        list_circ_orf_id_function = []
        results_buffer = []
        if not dic_id_peptide:
            return
        valid_source_ids = set(dic_id_peptide.keys())
        if not valid_source_ids:
            return
        print(f" -> Found {len(valid_source_ids)} unique source ORFs with peptide support.", flush=True)

        dic_orf_seqs_fallback = {}
        if os.path.exists(circ_orf):
            cur_id, cur_seq, cur_base = '', '', ''
            with open(circ_orf, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if cur_id and cur_seq:
                            dic_orf_seqs_fallback[cur_id] = cur_seq
                            dic_orf_seqs_fallback[cur_base] = cur_seq
                        cur_id = line[1:].strip()
                        cur_base = cur_id.split('(')[0]
                        cur_seq = ''
                    else:
                        cur_seq += line
                if cur_id and cur_seq:
                    dic_orf_seqs_fallback[cur_id] = cur_seq
                    dic_orf_seqs_fallback[cur_base] = cur_seq


        precomputed = {}
        candidate_orfs = []
        for source_header in valid_source_ids:
            if source_header in list_circ_orf_id_function:
                continue
            list_circ_orf_id_function.append(source_header)
            source_seq = dic_source_aa.get(source_header, '')
            if not source_seq:
                source_base = source_header.split('(')[0]
                source_seq = dic_orf_seqs_fallback.get(source_header, '') or dic_orf_seqs_fallback.get(source_base, '')
            if not source_seq:
                for out_id in dic_source_search_orfs.get(source_header, []):
                    source_seq = dic_orf_seqs_fallback.get(out_id, '') or dic_orf_seqs_fallback.get(out_id.split('(')[0], '')
                    if source_seq:
                        break
            if not source_seq:
                continue
            (pro_number_score, marge_peptide_bsj_multiple, corf_peptide_vatio,
             _, raw_sum, number_all, mgf_score_all, detailed_stats,
             unique_peptide_count, peptide_seq_str) = circ_orf_tsv_peptide_map(source_header, source_seq)
            if marge_peptide_bsj_multiple == 'none' and corf_peptide_vatio == 'none':
                continue
            # Extract metadata from source header
            colon_parts = source_header.split(':')
            if len(colon_parts) >= 3:
                gene = colon_parts[1] if colon_parts[1] else 'unknown'
            else:
                gene = 'unknown'
            base_id = source_header.split('(')[0] if '(' in source_header else source_header
            orf_base = base_id.split('-')[0] if '-' in base_id else base_id
            circ_info_id = extract_circ_info_id(orf_base)
            circ_strand = dic_circ_strand.get(circ_info_id, 'NA')
            # BSJ location from pre-built dict
            bsj_base = source_header.split('(')[0]
            bsj_location = dic_circ_orf_bsj_location.get(bsj_base, '0')
            # Search ORF IDs
            search_orf_ids = dic_source_search_orfs.get(source_header, [source_header])
            search_orf_str = ','.join(search_orf_ids)
            precomputed[source_header] = {
                'pro_number_score': pro_number_score,
                'marge_peptide_bsj_multiple': marge_peptide_bsj_multiple,
                'corf_peptide_vatio': corf_peptide_vatio,
                'bsj_location': bsj_location,
                'raw_sum': raw_sum,
                'number_all': number_all,
                'mgf_score_all': mgf_score_all,
                'detailed_stats': detailed_stats,
                'circ_info_id': circ_info_id,
                'circ_strand': circ_strand,
                'unique_peptide_count': unique_peptide_count,
                'peptide_seq_str': peptide_seq_str,
                'gene': gene,
                'search_orf_ids': search_orf_str,
            }
            candidate_orfs.append((source_header, source_seq))
        if not candidate_orfs:
            return
        accepted_ids = set(x[0] for x in candidate_orfs)
        # MS sample report
        ms_orf_report_file = file_out + 'ms_sample_orf_identification.txt'
        print(f" -> Generating MS sample level ORF identification report: {ms_orf_report_file}...")
        all_ms_files = set()
        ms_to_orfs = {}
        for orf_id, _ in candidate_orfs:
            if orf_id not in dic_id_peptide:
                continue
            entries = dic_id_peptide[orf_id].strip(';').split(';')
            for entry in entries:
                if not entry:
                    continue
                parts = entry.split('###')
                if len(parts) < 2:
                    continue
                ms_id = parts[0]
                all_ms_files.add(ms_id)
                if ms_id not in ms_to_orfs:
                    ms_to_orfs[ms_id] = set()
                ms_to_orfs[ms_id].add(orf_id)
        with open(ms_orf_report_file, 'w+') as f_rep:
            f_rep.write("MS_Sample_ID\tIdentified_ORF_Count\tORF_IDs(semicolon_separated)\n")
            for ms_id in sorted(list(all_ms_files)):
                orfs = ms_to_orfs.get(ms_id, set())
                count = len(orfs)
                orf_str = ";".join(sorted(list(orfs))) if orfs else ""
                f_rep.write(f"{ms_id}\t{count}\t{orf_str}\n")
        # IRES prediction
        dic_ires_results = {}
        if run_ires == 1:
            orf_ids = [oid for oid, _ in candidate_orfs]
            print(f" -> Performing batch IRES prediction for {len(orf_ids)} ORFs...")
            dic_ires_results = batch_ires_predict(orf_ids)
            print(f" -> IRES prediction complete, identified {len(dic_ires_results)} IRES signals")
        else:
            print(" -> Skipping IRES prediction (run_ires=0)")
        # m6A prediction
        print(f" -> Performing batch m6A prediction for {len(candidate_orfs)} ORFs...")
        dic_m6a_results = batch_m6a_predict(candidate_orfs)
        print(f" -> m6A prediction complete, identified {len(dic_m6a_results)} ORFs with m6A sites")
        # Score calculation
        print(" -> Integrating features and calculating coding function score...")
        for orf_id, orf_seq in candidate_orfs:
            res = process_orf(orf_id, orf_seq, dic_ires_results, precomputed, dic_m6a_results)
            if res:
                results_buffer.append(res)
        print(f" -> Feature integration complete, valid results: {len(results_buffer)}")
        calculate_and_write_scores(results_buffer, file_out + 'circ_predict.txt')

    def process_orf(circ_orf_id, circ_orf_seq, dic_ires_results, precomputed, dic_m6a_results):
        data = precomputed.get(circ_orf_id)
        if not data:
            return None
        if run_ires == 0:
            ires_data = {'p1': 0.0, 'p2': 0.0, 'p3': 0.0}
        else:
            oid = circ_orf_id.split()[0]
            ires_data = dic_ires_results.get(oid, {'p1': 0.0, 'p2': 0.0, 'p3': 0.0})
        p1 = ires_data.get('p1', 0.0)
        p2 = ires_data.get('p2', 0.0)
        p3 = ires_data.get('p3', 0.0)
        val_p1 = p1 if p1 >= 0 else 0.0
        val_p2 = p2 if p2 >= 0 else 0.0
        val_p3 = p3 if p3 >= 0 else 0.0
        weighted_ires_score = 0.5 * val_p1 + 0.3 * val_p2 + 0.2 * val_p3
        circ_gene = data.get('gene', 'unknown')
        circ_m6a_count = str(dic_m6a_results.get(circ_orf_id, '0'))
        orf_len = len(circ_orf_seq.strip('*'))
        str_p1 = f"{p1:.4f}" if p1 >= 0 else "na"
        str_p2 = f"{p2:.4f}" if p2 >= 0 else "na"
        str_p3 = f"{p3:.4f}" if p3 >= 0 else "na"
        ires_str = f"{str_p1}\t{str_p2}\t{str_p3}"
        search_orf_str = data.get('search_orf_ids', circ_orf_id)
        row_circ_predict = (
            f"{circ_orf_id}\t{search_orf_str}\t{data['circ_info_id']}\t{circ_gene}\t{data['circ_strand']}\t{ires_str}\t{circ_m6a_count}\t{circ_orf_seq}\t{orf_len}\t"
            f"{data['marge_peptide_bsj_multiple']}\t{data['corf_peptide_vatio']:.4f}\t{data['bsj_location']}\t{data['raw_sum']}\t{data['number_all']}\t{data['unique_peptide_count']}\t{data['peptide_seq_str']}\t{data['mgf_score_all']}"
        )
        try:
            m_val = float(circ_m6a_count)
        except:
            m_val = 0.0
        stats = {
            'IRES_SCORE': weighted_ires_score,
            'm': m_val,
            'cov': float(data['corf_peptide_vatio']),
            'L': len(circ_orf_seq.strip('*')),
            'MS': data['raw_sum'],
            'PEP': data['number_all'],
            'detailed_stats': data['detailed_stats']
        }
        return {'data_row': row_circ_predict, 'stats': stats}

    def round_up_to_nice_number(n):
        if n <= 0:
            return 1
        s = str(int(n))
        if len(s) == 1:
            return 10
        first_digit = int(s[0])
        rest = int(s[1:])
        if rest == 0:
            return int(n)
        return (first_digit + 1) * (10 ** (len(s) - 1))

    def calculate_global_max_stats(top_entries, dic_id_peptide_ref, dic_stripped_to_full, idx_corf_id, idx_corf_seq):
        max_pep_count = 0
        max_ms_count = 0
        for _, line_input in top_entries:
            parts = line_input.strip().split('\t')
            if len(parts) <= max(idx_corf_id, idx_corf_seq):
                continue
            corf_id = parts[idx_corf_id]
            key = None
            if corf_id in dic_id_peptide_ref:
                key = corf_id
            elif corf_id in dic_stripped_to_full:
                key = dic_stripped_to_full[corf_id]
            else:
                base = corf_id.split('(')[0]
                if base in dic_id_peptide_ref:
                    key = base
                elif base in dic_stripped_to_full:
                    key = dic_stripped_to_full[base]
            if not key:
                continue
            pep_info = dic_id_peptide_ref[key]
            corf_seq = parts[idx_corf_seq]
            corf_seq_clean = corf_seq.strip('*')
            seq_len = len(corf_seq_clean)
            counts = [0] * seq_len
            ms_sets = [set() for _ in range(seq_len)]
            pep_entries = pep_info.strip(';').split(';')
            for entry in pep_entries:
                if not entry:
                    continue
                subparts = entry.split('###')
                if len(subparts) < 2:
                    continue
                mgf_id = subparts[0]
                pep_seq = subparts[1]
                positions = []
                start = 0
                while True:
                    idx = corf_seq_clean.find(pep_seq, start)
                    if idx == -1:
                        break
                    positions.append(idx)
                    start = idx + 1
                for pos in positions:
                    end = pos + len(pep_seq)
                    for i in range(pos, end):
                        if i < seq_len:
                            counts[i] += 1
                            ms_sets[i].add(mgf_id)
            if counts:
                max_pep_count = max(max_pep_count, max(counts))
            current_max_ms = max(len(s) for s in ms_sets) if ms_sets else 0
            max_ms_count = max(max_ms_count, current_max_ms)
        return max_pep_count, max_ms_count

    def get_top_entries(file_input, count):
        entries = []
        try:
            with open(file_input, 'r') as f:
                header = f.readline()
                header_parts = header.strip().split('\t')
                try:
                    idx_score = header_parts.index('coding_function_score')
                except ValueError:
                    idx_score = len(header_parts) - 1
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) <= idx_score:
                        continue
                    try:
                        score = float(parts[idx_score])
                    except:
                        score = 0.0
                    entries.append((score, line))
        except Exception as e:
            print(f"[WARNING] Failed to read result file: {e}")
            return [], None
        entries.sort(key=lambda x: x[0], reverse=True)
        return entries[:count], header_parts

    def draw_msgf_second(file_input, excellent_count):
        if not os.path.exists(file_input):
            return
        peptide_map_dir = f"{file_out}Peptide_coverage_map/"
        if not os.path.exists(peptide_map_dir):
            os.makedirs(peptide_map_dir, exist_ok=True)
        top_entries, header_parts = get_top_entries(file_input, excellent_count)
        if not top_entries:
            return
        try:
            idx_corf_id = header_parts.index('circ_ORF_ID')
            idx_circ_gene = header_parts.index('gene')
            idx_corf_seq = header_parts.index('cORF_seq')
            idx_seq_cov = header_parts.index('peptide_coverage')
            idx_bsj_pos = header_parts.index('cORF_bsj_position')
        except ValueError:
            idx_corf_id = 0
            idx_circ_gene = 3
            idx_corf_seq = 8
            idx_seq_cov = 10
            idx_bsj_pos = 11
        list_start_end = ['(', ')', '|']
        count_drawn = 0
        dic_stripped_to_full = {}
        if dic_id_peptide:
            for key in dic_id_peptide:
                stripped = key.split('(')[0]
                if stripped != key:
                    dic_stripped_to_full[stripped] = key
        global_max_pep_count, global_max_ms_count = calculate_global_max_stats(
            top_entries, dic_id_peptide, dic_stripped_to_full, idx_corf_id, idx_corf_seq)
        scale_max_pep = round_up_to_nice_number(global_max_pep_count) if global_max_pep_count > 0 else 1
        scale_max_ms = round_up_to_nice_number(global_max_ms_count) if global_max_ms_count > 0 else 1
        for _, line_input in top_entries:
            parts = line_input.strip().split('\t')
            max_idx = max(idx_corf_id, idx_circ_gene, idx_corf_seq, idx_seq_cov, idx_bsj_pos)
            if len(parts) <= max_idx:
                continue
            corf_id = parts[idx_corf_id]
            circ_gene = parts[idx_circ_gene]
            corf_seq = parts[idx_corf_seq]
            try:
                float_circ_seq_coverage = float(parts[idx_seq_cov]) * 100
            except:
                float_circ_seq_coverage = 0.0
            circ_orf_id_key = None
            if corf_id in dic_id_peptide:
                circ_orf_id_key = corf_id
            elif corf_id in dic_stripped_to_full:
                circ_orf_id_key = dic_stripped_to_full[corf_id]
            else:
                base_id = corf_id.split('(')[0]
                if base_id in dic_id_peptide:
                    circ_orf_id_key = base_id
                elif base_id in dic_stripped_to_full:
                    circ_orf_id_key = dic_stripped_to_full[base_id]
            if not circ_orf_id_key:
                continue
            pep_info_str = dic_id_peptide[circ_orf_id_key]
            corf_seq_clean = corf_seq.strip('*')
            corf_seq_len = len(corf_seq_clean)
            site_pep_counts = [0] * corf_seq_len
            site_ms_sets = [set() for _ in range(corf_seq_len)]
            pep_entries = pep_info_str.strip(';').split(';')
            for entry in pep_entries:
                if not entry:
                    continue
                subparts = entry.split('###')
                if len(subparts) < 2:
                    continue
                mgf_id = subparts[0]
                pep_seq = subparts[1]
                positions = []
                start = 0
                while True:
                    idx = corf_seq_clean.find(pep_seq, start)
                    if idx == -1:
                        break
                    positions.append(idx)
                    start = idx + 1
                for pos in positions:
                    end = pos + len(pep_seq)
                    for i in range(pos, end):
                        if i < corf_seq_len:
                            site_pep_counts[i] += 1
                            site_ms_sets[i].add(mgf_id)
            list_pep_count_multiple = site_pep_counts
            list_ms_count_multiple = [len(s) for s in site_ms_sets]
            if corf_seq_len <= 20:
                column_num = max(corf_seq_len, 10)
            else:
                column_num = max(20, min(60, corf_seq_len // 10))
            corf_list_int = corf_seq_len // column_num
            corf_list_remainder = corf_seq_len - corf_list_int * column_num
            x = list(range(column_num)) * corf_list_int + list(range(corf_list_remainder))
            yyy_list = list(range(1, corf_list_int + 1))
            yyy_list.sort(reverse=True)
            y_list = []
            for yy in yyy_list:
                y_list.extend([yy] * column_num)
            y = y_list + [0] * corf_list_remainder
            norm = colors.Normalize(vmin=0, vmax=scale_max_pep)
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list('custom_teal', ['#13bbe9', '#035166'])
            fig_width = max(12, column_num * 0.4)
            fig = plt.figure(figsize=(fig_width, 10))
            plt.subplots_adjust(right=0.8)
            ax = fig.add_subplot(111)
            font_size = max(4, min(10, 300 / column_num))
            for aa in range(len(corf_seq_clean)):
                if aa < len(x) and aa < len(y):
                    ax.text(x[aa], y[aa], corf_seq_clean[aa], fontsize=font_size, color='red',
                            ha='center', va='center', zorder=10)
            try:
                bsj_positions = parts[idx_bsj_pos].split('|')
                for pos_str in bsj_positions:
                    try:
                        pos = int(float(pos_str))
                        if pos < len(x) and pos < len(y):
                            ax.text(x[pos] - 0.5, y[pos], list_start_end[2],
                                    fontsize=600 / column_num, color='black',
                                    ha='center', va='center', zorder=11)
                    except:
                        pass
            except:
                pass
            base_s = 50 * 1.5
            max_s = (600 * 1.5) / 1.2
            valid_indices = [i for i, c in enumerate(list_pep_count_multiple) if c > 0]
            if valid_indices:
                plot_x = [x[i] for i in valid_indices]
                plot_y = [y[i] for i in valid_indices]
                plot_c = [list_pep_count_multiple[i] for i in valid_indices]
                plot_ms = [list_ms_count_multiple[i] for i in valid_indices]
                plot_s = []
                plot_c_clipped = []
                for ms_c in plot_ms:
                    effective_ms = min(ms_c, scale_max_ms)
                    ratio = effective_ms / scale_max_ms if scale_max_ms > 0 else 0
                    if ratio > 1:
                        ratio = 1
                    size = base_s + ratio * (max_s - base_s)
                    plot_s.append(size)
                for pep_c in plot_c:
                    effective_pep = min(pep_c, scale_max_pep)
                    plot_c_clipped.append(effective_pep)
                segments = []
                line_widths = []
                line_colors = []
                for k in range(len(valid_indices) - 1):
                    idx_curr = valid_indices[k]
                    idx_next = valid_indices[k + 1]
                    if idx_next == idx_curr + 1:
                        if plot_y[k] == plot_y[k + 1]:
                            segments.append([(plot_x[k], plot_y[k]), (plot_x[k + 1], plot_y[k + 1])])
                            avg_s = (plot_s[k] + plot_s[k + 1]) / 2
                            width = 2 * np.sqrt(avg_s / np.pi)
                            line_widths.append(width)
                            avg_c = (plot_c_clipped[k] + plot_c_clipped[k + 1]) / 2
                            line_colors.append(avg_c)
                if segments:
                    lc = LineCollection(segments, linewidths=line_widths, cmap=cmap, norm=norm,
                                        alpha=0.9, zorder=4)
                    lc.set_array(np.array(line_colors))
                    ax.add_collection(lc)
                sc = ax.scatter(plot_x, plot_y, s=plot_s, marker='o', c=plot_c_clipped,
                                norm=norm, cmap=cmap, alpha=0.9, zorder=5)
                cax_x = 0.82
                cax_width = 0.02
                cax = fig.add_axes([cax_x, 0.45, cax_width, 0.3])
                cbar = plt.colorbar(sc, cax=cax)
                cbar.set_label('Coverage Peptide Count', fontsize=10)
                ticks = [scale_max_pep * 0.3, scale_max_pep * 0.6, scale_max_pep * 0.9]
                cbar.set_ticks(ticks)
                cbar.set_ticklabels([f'{int(t)}' for t in ticks])
                ms_ticks = [0.2, 0.4, 0.6, 0.8]
                legend_elements = []
                for t in ms_ticks:
                    val = int(scale_max_ms * t)
                    if val == 0:
                        continue
                    size = base_s + t * (max_s - base_s)
                    legend_elements.append(plt.scatter([], [], s=size, facecolors='none',
                                                       edgecolors='#035166', alpha=0.9, marker='o',
                                                       label=f'{val}'))
                if legend_elements:
                    ax.legend(handles=legend_elements, title='Coverage Peptide-MS Count',
                              bbox_to_anchor=(0.95, 0.08), loc='lower left',
                              borderaxespad=0., labelspacing=1.5, frameon=False)
            ax.set_xlim(xmin=-1, xmax=column_num + 1)
            y_max = 2 if corf_list_int <= 20 else corf_list_int * 0.1
            ax.set_ylim(ymax=corf_list_int + y_max, ymin=-y_max * 2)
            ax.text(1, -y_max / 2, f'cORF:{corf_id}', fontsize=12)
            ax.text(1, -y_max / 1, f'Gene_name:{circ_gene}', fontsize=12)
            ax.text(1, -y_max / 0.66, f'Sequence cover:{float_circ_seq_coverage:.2f}%', fontsize=12)
            ax.axis('off')
            save_name = corf_id
            match = re.search(r'(.*-ORF\d+)_', corf_id)
            if match:
                save_name = match.group(1) + "_coverage"
            else:
                save_name = corf_id.replace('|', '_').replace(':', '_')
            plt.savefig(f"{peptide_map_dir}{save_name}.pdf", bbox_inches='tight')
            plt.close()
            count_drawn += 1

    def draw_circ_corf():
        if not os.path.exists(file_out + 'circ_predict.txt'):
            return
        top_entries, header_parts = get_top_entries(file_out + 'circ_predict.txt', 10)
        if not top_entries:
            return
        temp_request_map = {}
        try:
            idx_circ_gene = header_parts.index('gene')
            idx_full_orf = header_parts.index('circ_ORF_ID')
            idx_cORF_seq = header_parts.index('cORF_seq')
            idx_cORF_map = header_parts.index('cORF_map')
        except ValueError:
            idx_circ_gene = 3
            idx_full_orf = 0
            idx_cORF_seq = 8
            idx_cORF_map = 10
        for p_val, row_circ_predict in top_entries:
            parts = row_circ_predict.strip().split('\t')
            if len(parts) <= max(idx_circ_gene, idx_full_orf, idx_cORF_seq, idx_cORF_map):
                continue
            circ_gene = parts[idx_circ_gene]
            full_orf_id = parts[idx_full_orf]
            circ_corf_seq = parts[idx_cORF_seq]
            circ_corf_seq_x = parts[idx_cORF_map]
            value = f"{circ_gene}###{circ_corf_seq}###{circ_corf_seq_x}###{full_orf_id}"
            base_key = full_orf_id.split('-')[0]
            keys_to_try = [base_key]
            for prefix in ['P_base_', 'P_self_', 'F_base_', 'F_self_', 'F_exon_',
                           'R_base_', 'R_self_', 'R_exon_', 'F_noexonic_', 'R_noexonic_',
                           'F_unknown_', 'R_unknown_']:
                if base_key.startswith(prefix):
                    simple = base_key[len(prefix):]
                    keys_to_try.append(simple)
                    break
            for k in keys_to_try:
                if k not in temp_request_map:
                    temp_request_map[k] = []
                if value not in temp_request_map[k]:
                    temp_request_map[k].append(value)
        dic_circ_corf_txt = {}
        if os.path.exists(circ_seq):
            with open(circ_seq, 'r') as f:
                current_circ_id = ''
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        current_circ_id = line[1:].split(':')[0]
                        if current_circ_id in temp_request_map:
                            if current_circ_id not in dic_circ_corf_txt:
                                dic_circ_corf_txt[current_circ_id] = []
                            for v in temp_request_map[current_circ_id]:
                                if v not in dic_circ_corf_txt[current_circ_id]:
                                    dic_circ_corf_txt[current_circ_id].append(v)
        for k in dic_circ_corf_txt:
            dic_circ_corf_txt[k] = ";".join(dic_circ_corf_txt[k])
        circ_map_dir = f"{file_out}circRNA_ORF_map/"
        if not os.path.exists(circ_map_dir):
            os.makedirs(circ_map_dir, exist_ok=True)
        for circ_rna_key in dic_circ_corf_txt:
            if circ_rna_key not in dic_circ_id_seq:
                continue
            circ_sequence = dic_circ_id_seq[circ_rna_key]
            circ_seq_len = len(circ_sequence)
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            limit = 6.0
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)
            ax.axis('off')
            R = 5.0
            theta = np.linspace(0, 2 * np.pi, 1000)
            x_circ = R * np.cos(theta)
            y_circ = R * np.sin(theta)
            ax.plot(x_circ, y_circ, color='black', linewidth=4.5, zorder=1)
            ax.plot(0, R, marker='>', markersize=18, color='black', clip_on=False, zorder=10)
            ax.text(0, R + 0.6, 'BSJ', ha='center', va='center', fontsize=12,
                    fontweight='bold', color='black')
            tick_percents = [0.2, 0.4, 0.6, 0.8]
            for p in tick_percents:
                bp_val = int(round(circ_seq_len * p / 10.0)) * 10
                angle = np.pi / 2 - (p * 2 * np.pi)
                x_out = R * np.cos(angle)
                y_out = R * np.sin(angle)
                x_in = (R - 0.2) * np.cos(angle)
                y_in = (R - 0.2) * np.sin(angle)
                ax.plot([x_out, x_in], [y_out, y_in], color='black', linewidth=1.5)
                text_r = R - 0.6
                x_txt = text_r * np.cos(angle)
                y_txt = text_r * np.sin(angle)
                ax.text(x_txt, y_txt, f"{bp_val}bp", ha='center', va='center',
                        fontsize=6, fontweight='bold', family='sans-serif', rotation=0)
            orf_list = dic_circ_corf_txt[circ_rna_key].split(';')
            gene_name_display = "NA"
            has_gene = False
            for orf_info in orf_list:
                parts_orf = orf_info.split('###')
                if len(parts_orf) < 4:
                    continue
                c_gene = parts_orf[0]
                c_seq = parts_orf[1]
                c_map = parts_orf[2]
                c_full_id = parts_orf[3]
                if c_gene not in ['NA', 'nan', 'None', '', 'unknown']:
                    gene_name_display = c_gene
                    has_gene = True
                source_orf_id = dic_output_to_source.get(c_full_id, c_full_id)
                m = re.search(r'\((-?\d+),(-?\d+)\)', source_orf_id)
                if not m:
                    continue
                start_pos = int(m.group(1))
                if circ_seq_len == 0:
                    continue
                orf_start_abs = start_pos % circ_seq_len
                base_r_orf = 4.2
                r_drop_per_turn = 0.36
                points = []
                coverage_status = []
                clean_map = c_map.replace('|', '')
                for k in range(len(c_seq)):
                    dist_from_start_nt = k * 3
                    current_abs_nt = (orf_start_abs + dist_from_start_nt) % circ_seq_len
                    angle = np.pi / 2 - (current_abs_nt / circ_seq_len) * 2 * np.pi
                    r = base_r_orf - (dist_from_start_nt / circ_seq_len) * r_drop_per_turn
                    points.append((r * np.cos(angle), r * np.sin(angle)))
                    if k < len(clean_map):
                        char = clean_map[k]
                        if char.isupper():
                            coverage_status.append(1)
                        else:
                            coverage_status.append(0)
                    else:
                        coverage_status.append(0)
                if not points:
                    continue
                xs = [p[0] for p in points]
                ys = [p[1] for p in points]
                ax.plot(xs, ys, color='#E0E0E0', linewidth=6, zorder=2, solid_capstyle='round')
                current_segment_indices = []
                for i in range(len(coverage_status)):
                    if coverage_status[i] == 1:
                        current_segment_indices.append(i)
                    else:
                        if current_segment_indices:
                            seg_xs = [points[ix][0] for ix in current_segment_indices]
                            seg_ys = [points[ix][1] for ix in current_segment_indices]
                            if len(current_segment_indices) == 1:
                                ax.scatter(seg_xs[0], seg_ys[0], color='#13bbe9', s=40, zorder=3)
                            else:
                                ax.plot(seg_xs, seg_ys, color='#13bbe9', linewidth=6,
                                        zorder=3, solid_capstyle='round')
                        current_segment_indices = []
                if current_segment_indices:
                    seg_xs = [points[ix][0] for ix in current_segment_indices]
                    seg_ys = [points[ix][1] for ix in current_segment_indices]
                    if len(current_segment_indices) == 1:
                        ax.scatter(seg_xs[0], seg_ys[0], color='#13bbe9', s=40, zorder=3)
                    else:
                        ax.plot(seg_xs, seg_ys, color='#13bbe9', linewidth=6,
                                zorder=3, solid_capstyle='round')
            display_id = circ_rna_key
            if not has_gene:
                ax.text(0, 0, display_id, ha='center', va='center', fontsize=14,
                        fontweight='bold', color='black')
            else:
                ax.text(0, 0.5, display_id, ha='center', va='center', fontsize=14,
                        fontweight='bold', color='black')
                ax.text(0, -0.5, gene_name_display, ha='center', va='center',
                        fontsize=14, style='italic', color='black')
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='black', lw=4.5, label='circRNA'),
                Line2D([0], [0], color='#E0E0E0', lw=6, label='ORF(Predicted)'),
                Line2D([0], [0], color='#13bbe9', lw=6, label='ORF(peptide)')
            ]
            import matplotlib.font_manager as fm
            font_prop = fm.FontProperties(family='Arial', size=10)
            ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(1.0, 0.0),
                      frameon=False, prop=font_prop)
            save_name = circ_rna_key.replace('|', '_').replace(':', '_')
            plt.savefig(f"{circ_map_dir}{save_name}_circ.pdf", bbox_inches='tight')
            plt.close()

    # ==================== MAIN FLOW ====================
    print(BANNER)
    print(f"Output Directory: {file_out if file_out.startswith(('/', './')) else './' + file_out}")
    print(f"circRNA Result: {circ_file_input}")
    print(f"Mass Spec Input Directory: {path_ms_input}")
    print(f"Linear Protein File: {file_input_protein}")

    print("\n[Step 1/6] Loading circRNA full-length sequences...")
    dic_circ_id_seq = circ_id_seq()

    print("\n[Step 2/6] Parsing mass spec data and filtering linear RNA peptides...")
    dic_id_peptide = process_ms_directory()
    if dic_id_peptide is None:
        dic_id_peptide = {}
    if file_input_protein != 'none':
        print(" -> Performing linear RNA peptide filtering...")
        dic_id_peptide = tsv_remove_line_protaid()
        if dic_id_peptide:
            all_peptides_set = set()
            for circ_id, info_str in dic_id_peptide.items():
                entries = info_str.strip(';').split(';')
                for entry in entries:
                    if not entry:
                        continue
                    parts = entry.split('###')
                    if len(parts) >= 2:
                        all_peptides_set.add(parts[1])
            if all_peptides_set:
                blacklisted_peptides = batch_diamond_alignment_filter(list(all_peptides_set))
                if blacklisted_peptides:
                    cleaned_dic_id_peptide = {}
                    removed_pep_count = 0
                    for circ_id, info_str in dic_id_peptide.items():
                        entries = info_str.strip(';').split(';')
                        new_entries = []
                        for entry in entries:
                            if not entry:
                                continue
                            parts = entry.split('###')
                            if len(parts) >= 2:
                                pep = parts[1]
                                if pep in blacklisted_peptides:
                                    removed_pep_count += 1
                                    continue
                                new_entries.append(entry)
                        if new_entries:
                            cleaned_dic_id_peptide[circ_id] = ';'.join(new_entries)
                    print(f" -> Filtering complete. Removed {removed_pep_count} peptide entries.")
                    print(f" -> Remaining valid circORFs: {len(cleaned_dic_id_peptide)} (Original: {len(dic_id_peptide)})")
                    dic_id_peptide = cleaned_dic_id_peptide
                else:
                    print(" -> No homologous peptides found to filter.")
            else:
                print(" -> No valid peptide results after linear filtering.")
        elif dic_id_peptide:
            print(" -> No linear protein file specified, skipping linear filtering.")
        else:
            print(" -> No valid peptide results, skipping filtering.")

    # ---- Source ORF mapping ----
    print("\n[Step 3/6] Mapping to source ORFs and merging peptides...")
    dic_source_peptides, dic_source_search_orfs = map_to_source_orfs()
    dic_id_peptide = dic_source_peptides
    print(f" -> Mapped to {len(dic_id_peptide)} source ORFs.")
    for src, out_ids in dic_source_search_orfs.items():
        if len(out_ids) > 1:
            print(f"   Source {src.split('(')[0]}: merged from {len(out_ids)} search ORFs")

    print("\n[Step 4/6] Parsing BSJ positions and strand information...")
    dic_circ_orf_bsj_location = build_source_bsj_location()
    dic_circ_strand = circ_forward_reverse()

    print("\n[Step 5/6] Performing core analysis (ORF screening + IRES/m6A prediction)...", flush=True)
    circ_peptide_result()

    print("\n[Step 6/6] Peptide coverage visualization...")
    draw_msgf_second(file_out + 'circ_predict.txt', 10)
    draw_circ_corf()

    print("\nSuccessfully completed visualization of the top ten circRNAs!")
    print(f"\nAnalysis fully completed!")
    print(f"MScirctl2 analysis finished. Results output to: {file_out if file_out.startswith(('/', './')) else './' + file_out}")
    end_time_total = time.time()
    elapsed_time = end_time_total - start_time_total
    print(f"Total script runtime: {elapsed_time:.2f} seconds ({elapsed_time / 60:.2f} minutes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
