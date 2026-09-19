#!/usr/bin/env python3
import os
import sys
import csv
import re
import threading
import platform
from datetime import datetime
from collections import defaultdict

import pandas as pd

import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
from tkinter import ttk


def get_platform_info():
    return {
        'system': platform.system(),
        'release': platform.release(),
        'version': platform.version(),
        'python_version': platform.python_version()
    }


def safe_open_file(filename, mode='r', encoding='utf-8', errors='replace'):
    try:
        return open(filename, mode, encoding=encoding, errors=errors)
    except UnicodeDecodeError:
        try:
            return open(filename, mode, encoding=sys.getfilesystemencoding(), errors=errors)
        except Exception:
            if 'r' in mode:
                return open(filename, mode, encoding='latin-1', errors=errors)
            raise


def safe_write_text(filename, content, mode='w', encoding='utf-8'):
    try:
        with open(filename, mode, encoding=encoding) as f:
            f.write(content)
    except UnicodeEncodeError:
        with open(filename, mode, encoding=encoding, errors='replace') as f:
            f.write(content)
    except Exception:
        try:
            with open(filename, mode, encoding=sys.getfilesystemencoding(), errors='replace') as f:
                f.write(content)
        except Exception:
            with open(filename, mode + 'b') as f:
                f.write(content.encode('utf-8', errors='replace'))


def normalize_text(text):
    if isinstance(text, str):
        replacements = {
            '\u2011': '-', '\u2012': '-', '\u2013': '-', '\u2014': '-', '\u2015': '-',
            '\u2018': "'", '\u2019': "'", '\u201C': '"', '\u201D': '"',
            '\u2022': '*', '\u00A0': ' ',
        }
        for unicode_char, ascii_char in replacements.items():
            text = text.replace(unicode_char, ascii_char)
    return text


def get_file_encoding(filename):
    for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1', 'utf-16']:
        try:
            with open(filename, 'r', encoding=encoding) as f:
                f.read()
            return encoding
        except (UnicodeDecodeError, UnicodeError):
            continue
    return sys.getfilesystemencoding()


class GFFBasedFilter:
    def __init__(self, gff_file, hints_file, eggnog_file, interpro_file, output_dir,
                 min_score=1.0, margin=500, min_protein_len=80, min_exons=1,
                 min_braker_score=0.0, use_abinitio=False, run_analysis=True,
                 callback_log=None):

        self.gff_file = gff_file
        self.hints_file = hints_file
        self.eggnog_file = eggnog_file
        self.interpro_file = interpro_file
        self.output_dir = output_dir

        self.min_score = min_score
        self.margin = margin
        self.min_protein_len = min_protein_len
        self.min_exons = min_exons
        self.min_braker_score = min_braker_score if use_abinitio else 0.0
        self.use_abinitio = use_abinitio
        self.run_analysis = run_analysis
        self.callback_log = callback_log

        self.genes = {}
        self.transcripts = {}
        self.hints_index = defaultdict(list)
        self.transcript_ids = set()
        self.gene_to_transcripts = defaultdict(list)

        self.eggnog_header = []
        self.eggnog_lines = []
        self.interpro_lines = []
        self.eggnog_transcripts = set()
        self.interpro_transcripts = set()

        self.transcripts_with_hints = set()
        self.transcripts_without_hints = set()
        self.transcripts_quality_pass = set()

        self.confident_transcripts = set()
        self.eggnog_confident = set()
        self.interpro_confident = set()
        self.orphan_confident = set()
        self.confident_genes = set()

        os.makedirs(output_dir, exist_ok=True)

        platform_info = get_platform_info()
        self.log(f"\nPlatform: {platform_info['system']} {platform_info['release']}")
        self.log(f"Python: {platform_info['python_version']}")

    def log(self, msg):
        msg = normalize_text(msg)
        if self.callback_log:
            try:
                self.callback_log(msg)
            except Exception:
                print(msg)
        else:
            print(msg)

    def parse_attributes(self, attr_string):
        attrs = {}
        if pd.isna(attr_string):
            return attrs
        attr_string = normalize_text(str(attr_string))
        for item in attr_string.split(';'):
            item = item.strip()
            if not item:
                continue
            if '=' in item:
                parts = item.split('=', 1)
                if len(parts) == 2:
                    attrs[parts[0].strip()] = parts[1].strip('"')
            elif ' ' in item:
                parts = item.split(' ', 1)
                if len(parts) == 2:
                    attrs[parts[0].strip()] = parts[1].strip('"')
        return attrs

    def load_hints_index(self):
        if self.use_abinitio or not os.path.exists(self.hints_file):
            self.log("\nAb initio mode enabled or hints file not found. Skipping indexation.")
            return True

        self.log("\nIndexing hints...")
        columns = ['seqname', 'source', 'feature', 'start', 'end',
                   'score', 'strand', 'phase', 'attribute']

        try:
            encoding = get_file_encoding(self.hints_file)
            self.log(f"  Detected encoding: {encoding}")
            df = pd.read_csv(self.hints_file, sep='\t', comment='#',
                             header=None, names=columns, dtype=str,
                             encoding=encoding)
            self.log(f"  {len(df):,} hints loaded")
        except Exception as e:
            self.log(f"  Error: {e}")
            return False

        for col in ['start', 'end', 'score']:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        for idx, row in df.iterrows():
            if pd.isna(row['start']) or pd.isna(row['end']):
                continue
            chr_name = str(row['seqname'])
            start = int(row['start'])
            end = int(row['end'])
            self.hints_index[chr_name].append({
                'start': start,
                'end': end,
                'feature': str(row['feature']),
                'score': float(row['score']) if not pd.isna(row['score']) else 0,
                'source': str(row['source'])
            })

        for chr_name in self.hints_index:
            self.hints_index[chr_name].sort(key=lambda x: x['start'])

        self.log(f"  {len(self.hints_index)} chromosomes indexed")
        return True

    def parse_gff(self):
        self.log("\nParsing GFF (genes and transcripts)...")
        columns = ['seqname', 'source', 'feature', 'start', 'end',
                   'score', 'strand', 'phase', 'attribute']

        try:
            encoding = get_file_encoding(self.gff_file)
            self.log(f"  Detected encoding: {encoding}")
            df = pd.read_csv(self.gff_file, sep='\t', comment='#',
                             header=None, names=columns, dtype=str,
                             encoding=encoding)
            self.log(f"  {len(df):,} lines loaded")
        except Exception as e:
            self.log(f"  Error: {e}")
            return False

        for col in ['start', 'end']:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        transcript_exon_count = defaultdict(int)
        transcript_cds_length = defaultdict(int)

        for idx, row in df.iterrows():
            feature = str(row['feature'])
            attrs = self.parse_attributes(row['attribute'])
            if feature == 'exon':
                parent = attrs.get('Parent') or attrs.get('parent')
                if parent:
                    transcript_exon_count[parent] += 1
            elif feature in ['CDS', 'cDNA_match']:
                parent = attrs.get('Parent') or attrs.get('parent')
                if parent:
                    start = int(row['start'])
                    end = int(row['end'])
                    transcript_cds_length[parent] += abs(end - start) + 1

        for idx, row in df.iterrows():
            feature = str(row['feature'])
            attrs = self.parse_attributes(row['attribute'])

            if feature == 'gene':
                gene_id = attrs.get('ID') or attrs.get('gene_id')
                if gene_id:
                    gene_id = normalize_text(gene_id)
                    self.genes[gene_id] = {
                        'id': gene_id,
                        'chr': normalize_text(str(row['seqname'])),
                        'start': int(row['start']),
                        'end': int(row['end']),
                        'strand': normalize_text(str(row['strand'])),
                        'transcripts': set(),
                    }

            elif feature in ['mRNA', 'transcript']:
                tid = attrs.get('ID') or attrs.get('transcript_id')
                gid = attrs.get('Parent') or attrs.get('parent') or attrs.get('gene_id')
                if tid and gid:
                    tid = normalize_text(tid)
                    gid = normalize_text(gid)
                    try:
                        score_val = float(row['score'])
                    except (ValueError, TypeError):
                        score_val = 0.0
                    cds_len = transcript_cds_length.get(tid, 0)
                    self.transcripts[tid] = {
                        'id': tid,
                        'gene_id': gid,
                        'chr': normalize_text(str(row['seqname'])),
                        'start': int(row['start']),
                        'end': int(row['end']),
                        'strand': normalize_text(str(row['strand'])),
                        'score': score_val,
                        'exon_count': transcript_exon_count.get(tid, 0),
                        'cds_length': cds_len,
                        'protein_len': cds_len // 3,
                        'has_hint': False,
                        'hints_count': 0,
                        'passes_quality': False,
                        'in_eggnog': False,
                        'in_interpro': False,
                        'is_orphan': False,
                        'confident': False,
                    }
                    self.transcript_ids.add(tid)
                    if gid in self.genes:
                        self.genes[gid]['transcripts'].add(tid)
                        self.gene_to_transcripts[gid].append(tid)

        self.log(f"  {len(self.genes):,} genes identified")
        self.log(f"  {len(self.transcripts):,} transcripts identified")
        return True

    def apply_hints_to_transcripts(self):
        if self.use_abinitio:
            self.log("\nAb initio mode: hints ignored, all transcripts treated as hinted.")
            for tid, tdata in self.transcripts.items():
                tdata['has_hint'] = True
                tdata['hints_count'] = 0
                self.transcripts_with_hints.add(tid)
            return True

        self.log("\nApplying hints to transcripts...")
        for tid, tdata in self.transcripts.items():
            chr_name = tdata['chr']
            t_start = tdata['start']
            t_end = tdata['end']
            hits = []
            if chr_name in self.hints_index:
                for hint in self.hints_index[chr_name]:
                    if hint['start'] <= t_end + self.margin and hint['end'] >= t_start - self.margin:
                        if hint['score'] >= self.min_score:
                            hits.append(hint)
            tdata['hints_count'] = len(hits)
            tdata['has_hint'] = len(hits) > 0
            if len(hits) > 0:
                self.transcripts_with_hints.add(tid)
            else:
                self.transcripts_without_hints.add(tid)

        total = len(self.transcripts) if self.transcripts else 1
        self.log(f"  Transcripts with hints:    {len(self.transcripts_with_hints):,} ({len(self.transcripts_with_hints)/total*100:.1f}%)")
        self.log(f"  Transcripts without hints: {len(self.transcripts_without_hints):,} ({len(self.transcripts_without_hints)/total*100:.1f}%)")
        return True

    def apply_quality_filters(self):
        self.log("\nApplying quality filters to transcripts...")
        for tid, tdata in self.transcripts.items():
            ok = True
            if tdata['protein_len'] < self.min_protein_len:
                ok = False
            if tdata['exon_count'] < self.min_exons:
                ok = False
            if self.use_abinitio and tdata['score'] < self.min_braker_score:
                ok = False
            tdata['passes_quality'] = ok
            if ok:
                self.transcripts_quality_pass.add(tid)

        total = len(self.transcripts) if self.transcripts else 1
        self.log(f"  Transcripts passing quality: {len(self.transcripts_quality_pass):,} ({len(self.transcripts_quality_pass)/total*100:.1f}%)")
        return True

    def match_transcript_id(self, protein_id):
        if not protein_id:
            return None
        pid = normalize_text(protein_id.strip())

        if pid in self.transcript_ids:
            return pid

        if '.' in pid:
            base = pid.rsplit('.', 1)[0]
            if base in self.transcript_ids:
                return base

        for suffix in ['.p1', '.P1', '.protein', '_p1', '_P1']:
            if pid.endswith(suffix):
                base = pid[:-len(suffix)]
                if base in self.transcript_ids:
                    return base

        if pid in self.gene_to_transcripts:
            return self.gene_to_transcripts[pid][0]

        match = re.search(r'(g\d+\.t\d+|g\d+|gene\d+|maker-[^\s]+)', pid, re.IGNORECASE)
        if match:
            base = match.group(1)
            if base in self.transcript_ids:
                return base

        return None

    def load_eggnog(self):
        self.log("\nReading EggNOG...")
        if not os.path.exists(self.eggnog_file):
            self.log(f"  File not found: {self.eggnog_file}")
            return False

        try:
            encoding = get_file_encoding(self.eggnog_file)
            with safe_open_file(self.eggnog_file, 'r', encoding=encoding) as f:
                lines = f.readlines()

            self.eggnog_header = [l for l in lines if l.startswith('#')]
            self.eggnog_lines = [l for l in lines if not l.startswith('#') and l.strip()]

            for line in self.eggnog_lines:
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                pid = normalize_text(fields[0].strip())
                tid = self.match_transcript_id(pid)
                if tid:
                    self.eggnog_transcripts.add(tid)

            self.log(f"  {len(self.eggnog_lines):,} EggNOG lines")
            self.log(f"  {len(self.eggnog_transcripts):,} transcripts matched to EggNOG")
            return True
        except Exception as e:
            self.log(f"  Error: {e}")
            return False

    def load_interpro(self):
        self.log("\nReading InterProScan...")
        if not os.path.exists(self.interpro_file):
            self.log(f"  File not found: {self.interpro_file}")
            return False

        try:
            encoding = get_file_encoding(self.interpro_file)
            with safe_open_file(self.interpro_file, 'r', encoding=encoding) as f:
                reader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
                self.interpro_lines = [row for row in reader if row and not row[0].startswith('#')]

            for row in self.interpro_lines:
                pid = normalize_text(row[0].strip())
                tid = self.match_transcript_id(pid)
                if tid:
                    self.interpro_transcripts.add(tid)

            self.log(f"  {len(self.interpro_lines):,} InterProScan rows")
            self.log(f"  {len(self.interpro_transcripts):,} transcripts matched to InterProScan")
            return True
        except Exception as e:
            self.log(f"  Error: {e}")
            return False

    def classify(self):
        self.log("\nClassifying transcripts (quality AND hints)...")

        for tid, tdata in self.transcripts.items():
            tdata['in_eggnog'] = tid in self.eggnog_transcripts
            tdata['in_interpro'] = tid in self.interpro_transcripts

            has_hint = tdata['has_hint'] or self.use_abinitio
            passes_quality = tdata['passes_quality']
            is_confident = has_hint and passes_quality
            has_annotation = tdata['in_eggnog'] or tdata['in_interpro']

            tdata['is_orphan'] = is_confident and not has_annotation
            tdata['confident'] = is_confident

            if is_confident:
                self.confident_transcripts.add(tid)
                if tdata['in_eggnog']:
                    self.eggnog_confident.add(tid)
                if tdata['in_interpro']:
                    self.interpro_confident.add(tid)
                if not has_annotation:
                    self.orphan_confident.add(tid)

        self.confident_genes = set()
        for tid in self.confident_transcripts:
            gid = self.transcripts[tid]['gene_id']
            self.confident_genes.add(gid)

        total_t = len(self.transcripts) if self.transcripts else 1
        total_g = len(self.genes) if self.genes else 1

        self.log("\n  Transcript-level results:")
        self.log(f"    Total transcripts:       {len(self.transcripts):,}")
        self.log(f"    Confident transcripts:   {len(self.confident_transcripts):,} ({len(self.confident_transcripts)/total_t*100:.1f}%)")
        self.log(f"    EggNOG confident:        {len(self.eggnog_confident):,}")
        self.log(f"    InterProScan confident:  {len(self.interpro_confident):,}")
        self.log(f"    Orphan confident:        {len(self.orphan_confident):,}")

        self.log("\n  Gene-level results:")
        self.log(f"    Total genes:             {len(self.genes):,}")
        self.log(f"    Confident genes:         {len(self.confident_genes):,} ({len(self.confident_genes)/total_g*100:.1f}%)")

        return True

    def save_lists(self):
        self.log("\nSaving gene and transcript lists...")

        tfile = os.path.join(self.output_dir, 'transcripts_confident.txt')
        with safe_open_file(tfile, 'w') as f:
            for tid in sorted(self.confident_transcripts):
                f.write(normalize_text(tid) + '\n')
        self.log(f"  {tfile} ({len(self.confident_transcripts):,} transcripts)")

        gfile = os.path.join(self.output_dir, 'genes_confident.txt')
        with safe_open_file(gfile, 'w') as f:
            for gid in sorted(self.confident_genes):
                f.write(normalize_text(gid) + '\n')
        self.log(f"  {gfile} ({len(self.confident_genes):,} genes)")

        nohint = os.path.join(self.output_dir, 'transcripts_without_hints.txt')
        with safe_open_file(nohint, 'w') as f:
            for tid in sorted(self.transcripts_without_hints):
                f.write(normalize_text(tid) + '\n')
        self.log(f"  {nohint} ({len(self.transcripts_without_hints):,} transcripts)")

        return tfile, gfile

    def filter_gff(self):
        self.log("\nWriting filtered GFF3 (only confident genes and transcripts)...")
        output_file = os.path.join(self.output_dir, 'braker_confident.gff3')

        try:
            encoding = get_file_encoding(self.gff_file)

            with safe_open_file(self.gff_file, 'r', encoding=encoding) as infile, \
                 safe_open_file(output_file, 'w') as outfile:

                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                        continue

                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 9:
                        continue

                    feature = fields[2]
                    attrs = self.parse_attributes(fields[8])
                    keep = False

                    if feature == 'gene':
                        gid = attrs.get('ID') or attrs.get('gene_id')
                        if gid and normalize_text(gid) in self.confident_genes:
                            keep = True

                    elif feature in ['mRNA', 'transcript']:
                        tid = attrs.get('ID') or attrs.get('transcript_id')
                        if tid and normalize_text(tid) in self.confident_transcripts:
                            keep = True

                    else:
                        parent = attrs.get('Parent') or attrs.get('parent')
                        if parent:
                            p = normalize_text(parent)
                            if p in self.confident_transcripts:
                                keep = True
                            elif p in self.confident_genes:
                                keep = True

                    if keep:
                        outfile.write(line)

            self.log(f"  {output_file}")
            return output_file
        except Exception as e:
            self.log(f"  Error: {e}")
            return None

    def save_eggnog_confident(self):
        self.log("\nWriting EggNOG confident file...")
        output_file = os.path.join(self.output_dir, 'eggnog_confident.tsv')

        try:
            kept = 0
            with safe_open_file(output_file, 'w') as f:
                f.writelines(self.eggnog_header)
                for line in self.eggnog_lines:
                    fields = line.strip().split('\t')
                    if len(fields) < 2:
                        continue
                    pid = normalize_text(fields[0].strip())
                    tid = self.match_transcript_id(pid)
                    if tid and tid in self.eggnog_confident:
                        f.write(line)
                        kept += 1

            self.log(f"  {output_file} ({kept:,} lines)")
            return output_file
        except Exception as e:
            self.log(f"  Error: {e}")
            return None

    def save_interpro_confident(self):
        self.log("\nWriting InterProScan confident file...")
        output_file = os.path.join(self.output_dir, 'interpro_confident.tsv')

        try:
            kept = 0
            with safe_open_file(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                for row in self.interpro_lines:
                    pid = normalize_text(row[0].strip())
                    tid = self.match_transcript_id(pid)
                    if tid and tid in self.interpro_confident:
                        writer.writerow(row)
                        kept += 1

            self.log(f"  {output_file} ({kept:,} rows)")
            return output_file
        except Exception as e:
            self.log(f"  Error: {e}")
            return None

    def save_orphans_confident(self):
        self.log("\nWriting orphan confident transcripts...")
        output_file = os.path.join(self.output_dir, 'orphans_confident.tsv')

        try:
            with safe_open_file(output_file, 'w') as f:
                f.write("transcript_id\tgene_id\tchromosome\tstart\tend\tstrand\t")
                f.write("protein_length\texon_count\tbraker_score\thints_count\n")
                for tid in sorted(self.orphan_confident):
                    t = self.transcripts[tid]
                    f.write(f"{tid}\t{t['gene_id']}\t{t['chr']}\t{t['start']}\t{t['end']}\t")
                    f.write(f"{t['strand']}\t{t['protein_len']}\t{t['exon_count']}\t")
                    f.write(f"{t['score']}\t{t['hints_count']}\n")

            self.log(f"  {output_file} ({len(self.orphan_confident):,} transcripts)")
            return output_file
        except Exception as e:
            self.log(f"  Error: {e}")
            return None

    def analyze_filtered_gff(self, filtered_gff_path):
        if not self.run_analysis or not os.path.exists(filtered_gff_path):
            return

        self.log("\nRunning GFF analysis...")
        try:
            class GFFAnalyzer:
                def __init__(self):
                    self.genes = defaultdict(dict)
                    self.chromosomes = set()
                    self.total_genes = 0
                    self.total_transcripts = 0
                    self.genes_without_start = set()
                    self.genes_without_stop = set()
                    self.genes_without_both = set()
                    self.single_exon_genes = set()
                    self.multi_exon_genes = set()
                    self.strand_distribution = defaultdict(int)
                    self.gene_lengths = []
                    self.transcripts_per_gene = []
                    self.gene_transcript_details = defaultdict(list)

                def parse_attributes(self, attr_str):
                    attrs = {}
                    attr_str = normalize_text(str(attr_str))
                    for item in attr_str.split(';'):
                        if '=' in item:
                            k, v = item.split('=', 1)
                            attrs[normalize_text(k.strip())] = normalize_text(v.strip())
                    return attrs

                def parse_gff(self, filename):
                    current_gene = None
                    current_chromosome = None
                    current_transcript = None
                    encoding = get_file_encoding(filename)

                    with safe_open_file(filename, 'r', encoding=encoding) as f:
                        for line in f:
                            line = normalize_text(line.strip())
                            if not line or line.startswith('#'):
                                continue
                            fields = line.split('\t')
                            if len(fields) < 9:
                                continue
                            chrom, source, feature, start, end, score, strand, phase, attrs_str = fields
                            start, end = int(start), int(end)
                            attrs = self.parse_attributes(attrs_str)
                            self.chromosomes.add(chrom)

                            if feature == 'gene':
                                gene_id = attrs.get('ID', '')
                                if gene_id:
                                    self.genes[chrom][gene_id] = {
                                        'id': gene_id,
                                        'chrom': chrom,
                                        'start': start,
                                        'end': end,
                                        'strand': strand,
                                        'has_start': False,
                                        'has_stop': False,
                                        'exon_count': 0,
                                        'cds_count': 0,
                                        'intron_count': 0,
                                        'transcript_count': 0,
                                        'transcripts': set()
                                    }
                                    current_gene = gene_id
                                    current_chromosome = chrom
                                    self.total_genes += 1
                                    current_transcript = None

                            elif feature in ['mRNA', 'transcript']:
                                transcript_id = attrs.get('ID', '')
                                if transcript_id and current_gene:
                                    if current_gene in self.genes[current_chromosome]:
                                        self.genes[current_chromosome][current_gene]['transcripts'].add(transcript_id)
                                        self.genes[current_chromosome][current_gene]['transcript_count'] += 1
                                        self.gene_transcript_details[current_gene].append({
                                            'id': transcript_id,
                                            'start': start,
                                            'end': end,
                                            'strand': strand,
                                            'exon_count': 0,
                                            'cds_count': 0
                                        })
                                    current_transcript = transcript_id
                                    self.total_transcripts += 1

                            elif feature in ['CDS', 'exon', 'intron', 'start_codon', 'stop_codon']:
                                if current_gene and current_chromosome:
                                    gene_data = self.genes[current_chromosome].get(current_gene)
                                    if gene_data:
                                        if feature == 'CDS':
                                            gene_data['cds_count'] += 1
                                        elif feature == 'exon':
                                            gene_data['exon_count'] += 1
                                        elif feature == 'intron':
                                            gene_data['intron_count'] += 1
                                        elif feature == 'start_codon':
                                            gene_data['has_start'] = True
                                        elif feature == 'stop_codon':
                                            gene_data['has_stop'] = True

                def analyze(self):
                    for chrom, genes_dict in self.genes.items():
                        for gene_id, data in genes_dict.items():
                            tc = data.get('transcript_count', 0)
                            self.transcripts_per_gene.append(tc)

                            if not data['has_start'] and not data['has_stop']:
                                self.genes_without_both.add(gene_id)
                            elif not data['has_start']:
                                self.genes_without_start.add(gene_id)
                            elif not data['has_stop']:
                                self.genes_without_stop.add(gene_id)

                            exon_count = data.get('exon_count', 0)
                            if exon_count == 1:
                                self.single_exon_genes.add(gene_id)
                            elif exon_count > 1:
                                self.multi_exon_genes.add(gene_id)

                            strand = data.get('strand', '')
                            if strand:
                                self.strand_distribution[strand] += 1

                            if 'start' in data and 'end' in data:
                                self.gene_lengths.append(data['end'] - data['start'] + 1)

                def _pct(self, value):
                    if self.total_genes <= 0:
                        return "0.0%"
                    return f"{value/self.total_genes*100:.1f}%"

                def generate_report(self):
                    lines = []
                    lines.append("=" * 80)
                    lines.append("GENOME ANNOTATION REPORT (FILTERED SET)")
                    lines.append("=" * 80)
                    lines.append("")
                    lines.append("1. GENERAL SUMMARY")
                    lines.append("-" * 40)
                    lines.append(f"Total scaffolds/chromosomes: {len(self.chromosomes)}")
                    lines.append(f"Total predicted genes: {self.total_genes}")
                    lines.append(f"Total transcripts: {self.total_transcripts}")
                    if self.total_genes > 0:
                        lines.append(f"Average transcripts per gene: {self.total_transcripts/self.total_genes:.2f}")
                    lines.append("")
                    lines.append("2. TRANSCRIPT DISTRIBUTION")
                    lines.append("-" * 40)
                    if self.transcripts_per_gene:
                        tc = self.transcripts_per_gene
                        lines.append(f"Minimum transcripts per gene: {min(tc)}")
                        lines.append(f"Maximum transcripts per gene: {max(tc)}")
                        lines.append(f"Mean transcripts per gene: {sum(tc)/len(tc):.2f}")
                        single = sum(1 for x in tc if x == 1)
                        multi = sum(1 for x in tc if x > 1)
                        high = sum(1 for x in tc if x >= 5)
                        lines.append(f"Genes with single transcript: {single} ({single/len(tc)*100:.1f}%)")
                        lines.append(f"Genes with multiple transcripts: {multi} ({multi/len(tc)*100:.1f}%)")
                        lines.append(f"Genes with >=5 transcripts: {high} ({high/len(tc)*100:.1f}%)")
                    lines.append("")
                    lines.append("3. GENE DISTRIBUTION PER CHROMOSOME")
                    lines.append("-" * 40)
                    for chrom in sorted(self.chromosomes):
                        cnt = len(self.genes.get(chrom, {}))
                        tcnt = sum(d.get('transcript_count', 0) for d in self.genes.get(chrom, {}).values())
                        lines.append(f"  {chrom}: {cnt} genes, {tcnt} transcripts")
                    lines.append("")
                    lines.append("4. POTENTIAL PSEUDOGENE ANALYSIS")
                    lines.append("-" * 40)
                    lines.append(f"Genes missing start codon: {len(self.genes_without_start)}")
                    lines.append(f"Genes missing stop codon: {len(self.genes_without_stop)}")
                    lines.append(f"Genes missing both: {len(self.genes_without_both)}")
                    total_pseudo = len(self.genes_without_start | self.genes_without_stop)
                    lines.append(f"Total potential pseudogenes: {total_pseudo}")
                    lines.append("")
                    lines.append("5. GENE STRUCTURE")
                    lines.append("-" * 40)
                    lines.append(f"Single-exon genes: {len(self.single_exon_genes)} ({self._pct(len(self.single_exon_genes))})")
                    lines.append(f"Multi-exon genes: {len(self.multi_exon_genes)} ({self._pct(len(self.multi_exon_genes))})")
                    total_exons = sum(g.get('exon_count', 0) for gd in self.genes.values() for g in gd.values())
                    avg_exons = total_exons / self.total_genes if self.total_genes else 0
                    lines.append(f"Average exons per gene: {avg_exons:.2f}")
                    total_cds = sum(g.get('cds_count', 0) for gd in self.genes.values() for g in gd.values())
                    avg_cds = total_cds / self.total_genes if self.total_genes else 0
                    lines.append(f"Average CDS per gene: {avg_cds:.2f}")
                    lines.append("")
                    lines.append("6. STRAND DISTRIBUTION")
                    lines.append("-" * 40)
                    total_strand = sum(self.strand_distribution.values())
                    for strand, cnt in sorted(self.strand_distribution.items()):
                        pct = cnt/total_strand*100 if total_strand else 0
                        lines.append(f"Strand {strand}: {cnt} ({pct:.1f}%)")
                    lines.append("")
                    if self.gene_lengths:
                        lines.append("7. GENE LENGTH STATISTICS")
                        lines.append("-" * 40)
                        lines.append(f"Minimum gene length: {min(self.gene_lengths)} bp")
                        lines.append(f"Maximum gene length: {max(self.gene_lengths)} bp")
                        lines.append(f"Average gene length: {sum(self.gene_lengths)/len(self.gene_lengths):.1f} bp")
                        lines.append(f"Median gene length: {sorted(self.gene_lengths)[len(self.gene_lengths)//2]} bp")
                        lines.append("")
                    lines.append("=" * 80)
                    lines.append("END OF REPORT")
                    lines.append("=" * 80)
                    return "\n".join(lines)

                def save_csv(self, filename):
                    with safe_open_file(filename, 'w') as f:
                        f.write("Gene_ID\tChromosome\tStart\tEnd\tStrand\tExons\tCDS\t")
                        f.write("Transcripts\tHas_Start\tHas_Stop\tLength\n")
                        for chrom, genes_dict in self.genes.items():
                            for gene_id, data in genes_dict.items():
                                f.write(f"{gene_id}\t{chrom}\t{data['start']}\t{data['end']}\t")
                                f.write(f"{data['strand']}\t{data['exon_count']}\t{data['cds_count']}\t")
                                f.write(f"{data['transcript_count']}\t")
                                f.write(f"{data['has_start']}\t{data['has_stop']}\t")
                                f.write(f"{data['end'] - data['start'] + 1}\n")

            analyzer = GFFAnalyzer()
            analyzer.parse_gff(filtered_gff_path)
            analyzer.analyze()

            report_text = normalize_text(analyzer.generate_report())
            report_file = os.path.join(self.output_dir, 'filtered_gff_analysis_report.txt')
            safe_write_text(report_file, report_text)
            self.log(f"  Analysis report: {report_file}")

            csv_file = os.path.join(self.output_dir, 'filtered_gene_details.csv')
            analyzer.save_csv(csv_file)
            self.log(f"  Gene details CSV: {csv_file}")

            summary_lines = report_text.split('\n')[:20]
            summary = "\n".join(summary_lines) + "\n...\n[Full report saved in file]"
            self.log("\n" + summary)

        except Exception as e:
            self.log(f"  Error during analysis: {e}")
            import traceback
            self.log(traceback.format_exc())

    def _pct(self, part, total):
        if total <= 0:
            return "0.0%"
        return f"{part/total*100:.1f}%"

    def generate_report(self):
        self.log("\nGenerating final report...")
        report_file = os.path.join(self.output_dir, 'pipeline_report.txt')

        total_t = len(self.transcripts)
        total_g = len(self.genes)
        union_annot = self.eggnog_transcripts | self.interpro_transcripts
        both = self.eggnog_confident & self.interpro_confident
        only_e = self.eggnog_confident - self.interpro_confident
        only_i = self.interpro_confident - self.eggnog_confident

        all_conf_genes = 0
        all_orphan_genes = 0
        any_hint_genes = 0
        for gid, gdata in self.genes.items():
            tids = gdata['transcripts']
            if not tids:
                continue
            conf_count = sum(1 for t in tids if t in self.confident_transcripts)
            orph_count = sum(1 for t in tids if t in self.orphan_confident)
            hint_count = sum(1 for t in tids if self.transcripts[t]['has_hint'])
            if conf_count == len(tids):
                all_conf_genes += 1
            if orph_count == len(tids):
                all_orphan_genes += 1
            if hint_count > 0:
                any_hint_genes += 1

        lines = []
        lines.append("=" * 80)
        lines.append("HINTSFILTER PIPELINE REPORT")
        lines.append("=" * 80)
        lines.append("")
        lines.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Platform: {platform.system()} {platform.release()}")
        lines.append(f"Python: {platform.python_version()}")
        lines.append("")
        lines.append("INPUT FILES:")
        lines.append(f"  GFF:          {self.gff_file}")
        lines.append(f"  Hints:        {self.hints_file}")
        lines.append(f"  EggNOG:       {self.eggnog_file}")
        lines.append(f"  InterProScan: {self.interpro_file}")
        lines.append("")
        lines.append("PARAMETERS:")
        lines.append(f"  Minimum hints score:    {self.min_score}")
        lines.append(f"  Margin (bp):            {self.margin}")
        lines.append(f"  Minimum protein length: {self.min_protein_len} aa")
        lines.append(f"  Minimum exon count:     {self.min_exons}")
        lines.append(f"  Minimum BRAKER score:   {self.min_braker_score} (only in ab initio mode)")
        lines.append(f"  Ab initio mode:         {self.use_abinitio}")
        lines.append(f"  Run GFF analysis:       {self.run_analysis}")
        lines.append("")

        lines.append("-" * 80)
        lines.append("TRANSCRIPT-LEVEL SUMMARY")
        lines.append("-" * 80)
        lines.append(f"Total transcripts in GFF:                 {total_t:,}")
        lines.append(f"  With hints:                             {len(self.transcripts_with_hints):,} ({self._pct(len(self.transcripts_with_hints), total_t)})")
        lines.append(f"  Passing quality filters:                {len(self.transcripts_quality_pass):,} ({self._pct(len(self.transcripts_quality_pass), total_t)})")
        lines.append(f"  With EggNOG annotation:                 {len(self.eggnog_transcripts):,} ({self._pct(len(self.eggnog_transcripts), total_t)})")
        lines.append(f"  With InterProScan annotation:           {len(self.interpro_transcripts):,} ({self._pct(len(self.interpro_transcripts), total_t)})")
        lines.append(f"  With functional annotation (union):     {len(union_annot):,} ({self._pct(len(union_annot), total_t)})")
        lines.append(f"  Orphans (quality, no annotation):       {len(self.orphan_confident):,} ({self._pct(len(self.orphan_confident), total_t)})")
        lines.append("")
        lines.append(f"CONFIDENT transcripts (final):            {len(self.confident_transcripts):,} ({self._pct(len(self.confident_transcripts), total_t)})")
        lines.append(f"  EggNOG only:                            {len(only_e):,}")
        lines.append(f"  InterProScan only:                      {len(only_i):,}")
        lines.append(f"  EggNOG and InterProScan:                {len(both):,}")
        lines.append(f"  Orphans:                                {len(self.orphan_confident):,}")
        lines.append("")

        lines.append("-" * 80)
        lines.append("GENE-LEVEL SUMMARY")
        lines.append("-" * 80)
        lines.append(f"Total genes in GFF:                       {total_g:,}")
        lines.append(f"  Genes with >=1 confident transcript:    {len(self.confident_genes):,} ({self._pct(len(self.confident_genes), total_g)})")
        lines.append(f"  Genes with ALL transcripts confident:   {all_conf_genes:,} ({self._pct(all_conf_genes, total_g)})")
        lines.append(f"  Genes with ALL transcripts orphan:      {all_orphan_genes:,} ({self._pct(all_orphan_genes, total_g)})")
        lines.append(f"  Genes with >=1 hinted transcript:       {any_hint_genes:,} ({self._pct(any_hint_genes, total_g)})")
        lines.append("")

        lines.append("FILTERING LOGIC:")
        lines.append("  A transcript is CONFIDENT if:")
        lines.append("    passes quality filters AND (has hint OR ab initio mode)")
        lines.append("  Confident transcripts are split by evidence into:")
        lines.append("    - EggNOG confident:   confident AND has EggNOG hit")
        lines.append("    - InterPro confident: confident AND has InterProScan hit")
        lines.append("    - Orphan confident:   confident AND no functional annotation")
        lines.append("  Orphan confident transcripts are candidates for lineage-specific")
        lines.append("  genes but require additional validation (they may also reflect")
        lines.append("  fragmented models, non-coding transcripts, or database gaps).")
        lines.append("")
        lines.append("OUTPUT FILES:")
        lines.append(f"  {self.output_dir}/genes_confident.txt")
        lines.append(f"  {self.output_dir}/transcripts_confident.txt")
        lines.append(f"  {self.output_dir}/transcripts_without_hints.txt")
        lines.append(f"  {self.output_dir}/braker_confident.gff3")
        lines.append(f"  {self.output_dir}/eggnog_confident.tsv")
        lines.append(f"  {self.output_dir}/interpro_confident.tsv")
        lines.append(f"  {self.output_dir}/orphans_confident.tsv")
        lines.append(f"  {self.output_dir}/pipeline_report.txt")
        if self.run_analysis:
            lines.append(f"  {self.output_dir}/filtered_gff_analysis_report.txt")
            lines.append(f"  {self.output_dir}/filtered_gene_details.csv")
        lines.append("")
        lines.append("=" * 80)

        report_text = normalize_text("\n".join(lines))
        safe_write_text(report_file, report_text)
        self.log(f"  {report_file}")

    def run(self):
        self.log("\n" + "=" * 80)
        self.log("HINTSFILTER PIPELINE (TRANSCRIPT-LEVEL CONFIDENCE)")
        self.log("=" * 80)

        if not self.load_hints_index():
            return False
        if not self.parse_gff():
            return False

        self.apply_hints_to_transcripts()
        self.apply_quality_filters()
        self.load_eggnog()
        self.load_interpro()
        self.classify()
        self.save_lists()

        filtered_gff = self.filter_gff()
        self.save_eggnog_confident()
        self.save_interpro_confident()
        self.save_orphans_confident()

        if filtered_gff and self.run_analysis:
            self.analyze_filtered_gff(filtered_gff)

        self.generate_report()

        self.log("\n" + "=" * 80)
        self.log("PIPELINE COMPLETED SUCCESSFULLY")
        self.log(f"Output directory: {self.output_dir}")
        self.log("=" * 80)
        return True


class PipelineGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("hintsFilter GFF")
        self.root.geometry("820x860")
        self.root.resizable(True, True)

        self.bg_color = "#F0F8F0"
        self.fg_color = "#1B5E20"
        self.frame_bg = "#E8F5E9"
        self.label_color = "#2E7D32"
        self.entry_bg = "#FFFFFF"
        self.btn_bg = "#A5D6A7"
        self.btn_active = "#81C784"
        self.btn_fg = "#1B5E20"
        self.text_bg = "#FFFFFF"
        self.text_fg = "#1B5E20"

        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TFrame', background=self.frame_bg)
        style.configure('TLabel', background=self.frame_bg, foreground=self.label_color, font=('Segoe UI', 9))
        style.configure('TButton', background=self.btn_bg, foreground=self.btn_fg, font=('Segoe UI', 9, 'bold'))
        style.map('TButton', background=[('active', self.btn_active)])
        style.configure('TEntry', fieldbackground=self.entry_bg, foreground=self.fg_color)

        self.gff_file = tk.StringVar()
        self.hints_file = tk.StringVar()
        self.eggnog_file = tk.StringVar()
        self.interpro_file = tk.StringVar()
        self.output_dir = tk.StringVar(value="gff_filtered_output")

        self.min_score = tk.DoubleVar(value=1.0)
        self.margin = tk.IntVar(value=500)
        self.min_protein_len = tk.IntVar(value=80)
        self.min_exons = tk.IntVar(value=1)
        self.min_braker_score = tk.DoubleVar(value=0.0)
        self.use_abinitio = tk.BooleanVar(value=False)
        self.run_analysis = tk.BooleanVar(value=True)

        self.running = False

        self.build_widgets()

    def build_widgets(self):
        main_frame = ttk.Frame(self.root, style='TFrame')
        main_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)

        title_label = tk.Label(main_frame, text="hintsFilter GFF (for BRAKER annotation)",
                               font=('Segoe UI', 14, 'bold'), bg=self.bg_color, fg=self.fg_color)
        title_label.grid(row=0, column=0, columnspan=4, pady=(0, 15), sticky='w')

        platform_info = get_platform_info()
        platform_label = tk.Label(main_frame,
                                 text=f"{platform_info['system']} | Python {platform_info['python_version']}",
                                 font=('Segoe UI', 8), bg=self.bg_color, fg=self.label_color)
        platform_label.grid(row=0, column=3, pady=(0, 15), sticky='e')

        tk.Label(main_frame, text="BRAKER prediction (GFF3/GTF):", font=('Segoe UI', 9),
                 bg=self.bg_color, fg=self.label_color).grid(row=1, column=0, sticky='e', padx=5, pady=5)
        tk.Entry(main_frame, textvariable=self.gff_file, width=50, bg=self.entry_bg, fg=self.fg_color).grid(row=1, column=1, padx=5, pady=5)
        tk.Button(main_frame, text="Browse...", bg=self.btn_bg, fg=self.btn_fg, activebackground=self.btn_active, command=self.browse_gff).grid(row=1, column=2, padx=5, pady=5)

        tk.Label(main_frame, text="Hints file (GFF):", font=('Segoe UI', 9),
                 bg=self.bg_color, fg=self.label_color).grid(row=2, column=0, sticky='e', padx=5, pady=5)
        tk.Entry(main_frame, textvariable=self.hints_file, width=50, bg=self.entry_bg, fg=self.fg_color).grid(row=2, column=1, padx=5, pady=5)
        tk.Button(main_frame, text="Browse...", bg=self.btn_bg, fg=self.btn_fg, activebackground=self.btn_active, command=self.browse_hints).grid(row=2, column=2, padx=5, pady=5)

        tk.Label(main_frame, text="EggNOG file:", font=('Segoe UI', 9),
                 bg=self.bg_color, fg=self.label_color).grid(row=3, column=0, sticky='e', padx=5, pady=5)
        tk.Entry(main_frame, textvariable=self.eggnog_file, width=50, bg=self.entry_bg, fg=self.fg_color).grid(row=3, column=1, padx=5, pady=5)
        tk.Button(main_frame, text="Browse...", bg=self.btn_bg, fg=self.btn_fg, activebackground=self.btn_active, command=self.browse_eggnog).grid(row=3, column=2, padx=5, pady=5)

        tk.Label(main_frame, text="InterProScan file:", font=('Segoe UI', 9),
                 bg=self.bg_color, fg=self.label_color).grid(row=4, column=0, sticky='e', padx=5, pady=5)
        tk.Entry(main_frame, textvariable=self.interpro_file, width=50, bg=self.entry_bg, fg=self.fg_color).grid(row=4, column=1, padx=5, pady=5)
        tk.Button(main_frame, text="Browse...", bg=self.btn_bg, fg=self.btn_fg, activebackground=self.btn_active, command=self.browse_interpro).grid(row=4, column=2, padx=5, pady=5)

        tk.Label(main_frame, text="Output directory:", font=('Segoe UI', 9),
                 bg=self.bg_color, fg=self.label_color).grid(row=5, column=0, sticky='e', padx=5, pady=5)
        tk.Entry(main_frame, textvariable=self.output_dir, width=50, bg=self.entry_bg, fg=self.fg_color).grid(row=5, column=1, padx=5, pady=5)
        tk.Button(main_frame, text="Browse...", bg=self.btn_bg, fg=self.btn_fg, activebackground=self.btn_active, command=self.browse_output).grid(row=5, column=2, padx=5, pady=5)

        param_frame = tk.LabelFrame(main_frame, text="Filtering Parameters",
                                    bg=self.frame_bg, fg=self.fg_color, font=('Segoe UI', 10, 'bold'))
        param_frame.grid(row=6, column=0, columnspan=4, pady=10, padx=5, sticky='ew')

        tk.Label(param_frame, text="Min. hints score:", bg=self.frame_bg, fg=self.label_color).grid(row=0, column=0, padx=5, pady=5, sticky='e')
        tk.Entry(param_frame, textvariable=self.min_score, width=8, bg=self.entry_bg, fg=self.fg_color).grid(row=0, column=1, padx=5, pady=5, sticky='w')

        tk.Label(param_frame, text="Margin (bp):", bg=self.frame_bg, fg=self.label_color).grid(row=0, column=2, padx=5, pady=5, sticky='e')
        tk.Entry(param_frame, textvariable=self.margin, width=8, bg=self.entry_bg, fg=self.fg_color).grid(row=0, column=3, padx=5, pady=5, sticky='w')

        tk.Label(param_frame, text="Min. protein (aa):", bg=self.frame_bg, fg=self.label_color).grid(row=1, column=0, padx=5, pady=5, sticky='e')
        tk.Entry(param_frame, textvariable=self.min_protein_len, width=8, bg=self.entry_bg, fg=self.fg_color).grid(row=1, column=1, padx=5, pady=5, sticky='w')

        tk.Label(param_frame, text="Min. exons:", bg=self.frame_bg, fg=self.label_color).grid(row=1, column=2, padx=5, pady=5, sticky='e')
        tk.Entry(param_frame, textvariable=self.min_exons, width=8, bg=self.entry_bg, fg=self.fg_color).grid(row=1, column=3, padx=5, pady=5, sticky='w')

        tk.Label(param_frame, text="BRAKER score >=:", bg=self.frame_bg, fg=self.label_color).grid(row=1, column=4, padx=5, pady=5, sticky='e')
        self.braker_score_entry = tk.Entry(param_frame, textvariable=self.min_braker_score, width=8,
                                           bg=self.entry_bg, fg=self.fg_color, state='disabled')
        self.braker_score_entry.grid(row=1, column=5, padx=5, pady=5, sticky='w')

        self.abinitio_check = tk.Checkbutton(param_frame, text="Ab initio mode (ignore hints)",
                                             variable=self.use_abinitio,
                                             bg=self.frame_bg, fg=self.fg_color, selectcolor=self.bg_color,
                                             font=('Segoe UI', 9), command=self.toggle_braker_score)
        self.abinitio_check.grid(row=2, column=0, columnspan=6, pady=5, sticky='w', padx=10)

        self.analysis_check = tk.Checkbutton(param_frame, text="Run GFF analysis (gene statistics)",
                                             variable=self.run_analysis,
                                             bg=self.frame_bg, fg=self.fg_color, selectcolor=self.bg_color,
                                             font=('Segoe UI', 9))
        self.analysis_check.grid(row=3, column=0, columnspan=6, pady=5, sticky='w', padx=10)

        self.run_btn = tk.Button(main_frame, text="RUN PIPELINE", command=self.run_pipeline,
                                 bg="#66BB6A", fg="white", font=('Segoe UI', 11, 'bold'),
                                 activebackground="#43A047", activeforeground="white",
                                 relief=tk.RAISED, bd=3, padx=20, pady=8)
        self.run_btn.grid(row=7, column=0, columnspan=4, pady=15)

        log_frame = tk.LabelFrame(main_frame, text="Execution Log",
                                  bg=self.frame_bg, fg=self.fg_color, font=('Segoe UI', 10, 'bold'))
        log_frame.grid(row=8, column=0, columnspan=4, sticky='nsew', padx=5, pady=5)

        self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, height=18,
                                                  bg=self.text_bg, fg=self.text_fg,
                                                  font=('Consolas', 9))
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        main_frame.grid_rowconfigure(8, weight=1)
        main_frame.grid_columnconfigure(1, weight=1)

    def toggle_braker_score(self):
        if self.use_abinitio.get():
            self.braker_score_entry.config(state='normal')
        else:
            self.braker_score_entry.config(state='disabled')

    def browse_gff(self):
        f = filedialog.askopenfilename(title="Select GFF/GTF file",
                                       filetypes=[("GFF/GTF", "*.gff *.gff3 *.gtf"), ("All files", "*.*")])
        if f:
            self.gff_file.set(f)

    def browse_hints(self):
        f = filedialog.askopenfilename(title="Select hints file",
                                       filetypes=[("GFF", "*.gff"), ("All files", "*.*")])
        if f:
            self.hints_file.set(f)

    def browse_eggnog(self):
        f = filedialog.askopenfilename(title="Select EggNOG file",
                                       filetypes=[("Annotations", "*.annotations"), ("All files", "*.*")])
        if f:
            self.eggnog_file.set(f)

    def browse_interpro(self):
        f = filedialog.askopenfilename(title="Select InterProScan file",
                                       filetypes=[("TSV", "*.tsv"), ("All files", "*.*")])
        if f:
            self.interpro_file.set(f)

    def browse_output(self):
        d = filedialog.askdirectory(title="Select output directory")
        if d:
            self.output_dir.set(d)

    def log_message(self, msg):
        try:
            msg = normalize_text(msg)
            self.log_text.insert(tk.END, msg + "\n")
            self.log_text.see(tk.END)
            self.root.update_idletasks()
        except Exception as e:
            print(f"Log error: {e}")
            print(msg)

    def run_pipeline(self):
        if self.running:
            return

        if not self.gff_file.get():
            messagebox.showerror("Error", "Please select a GFF file.")
            return
        if not self.use_abinitio.get() and not self.hints_file.get():
            messagebox.showerror("Error", "Please select a hints file or enable Ab initio mode.")
            return
        if not self.eggnog_file.get() or not self.interpro_file.get():
            messagebox.showerror("Error", "Please select both EggNOG and InterProScan files.")
            return

        self.log_text.delete(1.0, tk.END)
        self.run_btn.config(state=tk.DISABLED, text="PROCESSING...")
        self.running = True

        thread = threading.Thread(target=self._run_pipeline_thread, daemon=True)
        thread.start()

    def _run_pipeline_thread(self):
        try:
            min_braker = self.min_braker_score.get() if self.use_abinitio.get() else 0.0

            pipeline = GFFBasedFilter(
                gff_file=self.gff_file.get(),
                hints_file=self.hints_file.get(),
                eggnog_file=self.eggnog_file.get(),
                interpro_file=self.interpro_file.get(),
                output_dir=self.output_dir.get(),
                min_score=self.min_score.get(),
                margin=self.margin.get(),
                min_protein_len=self.min_protein_len.get(),
                min_exons=self.min_exons.get(),
                min_braker_score=min_braker,
                use_abinitio=self.use_abinitio.get(),
                run_analysis=self.run_analysis.get(),
                callback_log=self.log_message
            )
            success = pipeline.run()
            if success:
                self.root.after(0, lambda: messagebox.showinfo("Success", "Pipeline completed successfully!"))
            else:
                self.root.after(0, lambda: messagebox.showerror("Error", "Pipeline failed. Check the log."))
        except Exception as e:
            self.log_message(f"\nCRITICAL ERROR: {e}")
            import traceback
            self.log_message(traceback.format_exc())
            self.root.after(0, lambda: messagebox.showerror("Error", f"Critical error: {e}"))
        finally:
            self.root.after(0, self._finish_run)

    def _finish_run(self):
        self.run_btn.config(state=tk.NORMAL, text="RUN PIPELINE")
        self.running = False


def main():
    try:
        root = tk.Tk()
        app = PipelineGUI(root)
        root.mainloop()
    except Exception as e:
        print(f"Error starting GUI: {e}")
        import traceback
        traceback.print_exc()
        input("Press Enter to exit...")


if __name__ == "__main__":
    main()
