#!/usr/bin/env python3
"""
多序列比对与保守性分析模块
使用MAFFT进行多序列比对，分析序列保守性和结构保守性
"""

import os
import subprocess
import tempfile
from pathlib import Path
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment

# Biopython 新版本已移除 Bio.Align.Applications；这里不再依赖其封装
try:
    from Bio.Align.Applications import MafftCommandline  # type: ignore
except Exception:
    MafftCommandline = None
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import json

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ConservationAnalyzer:
    """保守性分析器"""
    
    def __init__(self):
        self.results_dir = Path("results/tables")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir = Path("results/figures")
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.alignments_dir = Path("data/alignments")
        self.alignments_dir.mkdir(parents=True, exist_ok=True)
        
        # 遗传密码表
        self.codon_table = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def perform_multiple_sequence_alignment(self, sequences: Dict[str, str], gene_name: str) -> Optional[MultipleSeqAlignment]:
        """
        使用MAFFT进行多序列比对
        Args:
            sequences: 序列字典 {序列名: 序列}
            gene_name: 基因名称
        Returns:
            比对结果
        """
        logger.info(f"开始{gene_name}基因的多序列比对...")
        
        try:
            # 创建临时FASTA文件
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_file:
                for seq_name, seq in sequences.items():
                    tmp_file.write(f">{seq_name}\n{seq}\n")
                tmp_file_path = tmp_file.name
            
            # 输出文件路径
            output_file = self.alignments_dir / f"{gene_name}_alignment.fasta"
            
            # 运行MAFFT（优先直接调用系统命令，避免依赖 Biopython 的旧封装）
            try:
                result = subprocess.run(
                    ["mafft", "--auto", "--thread", "4", tmp_file_path],
                    capture_output=True,
                    text=True,
                    timeout=600,
                )
            except FileNotFoundError:
                logger.error("未找到mafft命令，请先安装MAFFT")
                os.unlink(tmp_file_path)
                return None

            if result.returncode != 0:
                logger.error(f"{gene_name}多序列比对失败: {result.stderr}")
                os.unlink(tmp_file_path)
                return None

            stdout = result.stdout

            # 清理临时文件
            os.unlink(tmp_file_path)
            
            # 保存比对结果
            with open(output_file, 'w') as f:
                f.write(stdout)
            
            # 读取比对结果
            alignment = AlignIO.read(output_file, "fasta")
            
            logger.info(f"{gene_name}多序列比对完成，比对长度: {alignment.get_alignment_length()}")
            return alignment
            
        except Exception as e:
            logger.error(f"{gene_name}多序列比对失败: {e}")
            return None
    
    def calculate_sequence_conservation(self, alignment: MultipleSeqAlignment) -> Dict:
        """
        计算序列保守性
        Args:
            alignment: 多序列比对结果
        Returns:
            保守性分析结果
        """
        logger.info("计算序列保守性...")
        
        alignment_length = alignment.get_alignment_length()
        n_sequences = len(alignment)
        
        # 初始化统计变量
        conservation_scores = []
        position_stats = []
        gap_positions = []
        
        for pos in range(alignment_length):
            # 获取当前位置的所有碱基
            column = alignment[:, pos]
            
            # 统计碱基频率（忽略gap）
            bases = [base for base in column if base != '-']
            gaps = column.count('-')
            
            if bases:
                # 计算保守性分数（Shannon熵的倒数）
                base_counts = {}
                for base in bases:
                    base_counts[base] = base_counts.get(base, 0) + 1
                
                # 计算频率
                frequencies = [count / len(bases) for count in base_counts.values()]
                
                # 计算Shannon熵
                entropy = -sum(f * np.log2(f) for f in frequencies if f > 0)
                
                # 保守性分数（熵越小，保守性越高）
                max_entropy = np.log2(min(len(set(bases)), 4))  # 最大可能熵
                conservation_score = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
                
                # 最常见的碱基
                most_common_base = max(base_counts, key=base_counts.get)
                most_common_count = base_counts[most_common_base]
                consensus_percentage = most_common_count / len(bases) * 100
            else:
                # 全是gap
                conservation_score = 0
                most_common_base = '-'
                consensus_percentage = 0
            
            conservation_scores.append(conservation_score)
            
            position_stats.append({
                'position': pos + 1,
                'most_common_base': most_common_base,
                'consensus_percentage': consensus_percentage,
                'gap_percentage': gaps / n_sequences * 100,
                'conservation_score': conservation_score,
                'base_counts': base_counts if bases else {}
            })
            
            if gaps > n_sequences * 0.5:  # 超过50%是gap的位置
                gap_positions.append(pos)
        
        # 计算整体统计
        avg_conservation = np.mean(conservation_scores)
        high_conservation_threshold = np.percentile(conservation_scores, 80)
        high_conservation_positions = [i for i, score in enumerate(conservation_scores) 
                                     if score >= high_conservation_threshold]
        
        return {
            'alignment_length': alignment_length,
            'n_sequences': n_sequences,
            'average_conservation': avg_conservation,
            'conservation_scores': conservation_scores,
            'position_stats': position_stats,
            'high_conservation_threshold': high_conservation_threshold,
            'high_conservation_positions': high_conservation_positions,
            'gap_positions': gap_positions,
            'consensus_sequence': self._build_consensus_sequence(alignment)
        }
    
    def _build_consensus_sequence(self, alignment: MultipleSeqAlignment) -> str:
        """构建一致序列"""
        consensus = []
        for pos in range(alignment.get_alignment_length()):
            column = alignment[:, pos]
            bases = [base for base in column if base != '-']
            if bases:
                most_common = max(set(bases), key=bases.count)
                consensus.append(most_common)
            else:
                consensus.append('-')
        return ''.join(consensus)
    
    def analyze_codon_conservation(self, alignment: MultipleSeqAlignment) -> Dict:
        """
        分析密码子保守性（针对编码序列）
        Args:
            alignment: 多序列比对结果
        Returns:
            密码子保守性分析结果
        """
        logger.info("分析密码子保守性...")
        
        # 确保比对长度是3的倍数
        alignment_length = alignment.get_alignment_length()
        if alignment_length % 3 != 0:
            logger.warning("比对长度不是3的倍数，将忽略末尾不完整的密码子")
            alignment_length = (alignment_length // 3) * 3
        
        codon_conservation = []
        amino_acid_conservation = []
        
        for i in range(0, alignment_length, 3):
            codon_positions = [i, i+1, i+2]
            
            # 获取所有序列的密码子
            codons = []
            for record in alignment:
                codon = ''.join([record.seq[pos] for pos in codon_positions])
                if '-' not in codon and len(codon) == 3:
                    codons.append(codon.replace('T', 'U'))  # 转换为RNA
            
            if codons:
                # 分析密码子保守性
                codon_counts = {}
                for codon in codons:
                    codon_counts[codon] = codon_counts.get(codon, 0) + 1
                
                most_common_codon = max(codon_counts, key=codon_counts.get)
                codon_consensus_percentage = codon_counts[most_common_codon] / len(codons) * 100
                
                # 翻译为氨基酸
                amino_acids = []
                for codon in codons:
                    aa = self.codon_table.get(codon, 'X')  # X表示未知氨基酸
                    amino_acids.append(aa)
                
                # 分析氨基酸保守性
                aa_counts = {}
                for aa in amino_acids:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
                
                most_common_aa = max(aa_counts, key=aa_counts.get)
                aa_consensus_percentage = aa_counts[most_common_aa] / len(amino_acids) * 100
                
                codon_conservation.append({
                    'position': i // 3 + 1,
                    'codon_position': i + 1,
                    'most_common_codon': most_common_codon,
                    'codon_consensus_percentage': codon_consensus_percentage,
                    'most_common_amino_acid': most_common_aa,
                    'aa_consensus_percentage': aa_consensus_percentage,
                    'synonymous_variation': len(set(codons)) > 1 and len(set(amino_acids)) == 1
                })
                
                amino_acid_conservation.append(aa_consensus_percentage / 100)
        
        return {
            'codon_conservation': codon_conservation,
            'average_aa_conservation': np.mean(amino_acid_conservation) if amino_acid_conservation else 0,
            'synonymous_sites': [c for c in codon_conservation if c['synonymous_variation']],
            'n_synonymous_sites': len([c for c in codon_conservation if c['synonymous_variation']])
        }
    
    def identify_conserved_motifs(self, conservation_data: Dict, min_length: int = 6) -> List[Dict]:
        """
        识别保守基序
        Args:
            conservation_data: 保守性数据
            min_length: 最小基序长度
        Returns:
            保守基序列表
        """
        logger.info("识别保守基序...")
        
        scores = conservation_data['conservation_scores']
        consensus = conservation_data['consensus_sequence']
        threshold = conservation_data['high_conservation_threshold']
        
        motifs = []
        current_motif_start = None
        
        for i, score in enumerate(scores):
            if score >= threshold:
                if current_motif_start is None:
                    current_motif_start = i
            else:
                if current_motif_start is not None:
                    motif_length = i - current_motif_start
                    if motif_length >= min_length:
                        motif_sequence = consensus[current_motif_start:i]
                        
                        # 计算基序的平均保守性
                        motif_scores = scores[current_motif_start:i]
                        avg_score = np.mean(motif_scores)
                        
                        motifs.append({
                            'start': current_motif_start + 1,  # 转换为1-based
                            'end': i,
                            'length': motif_length,
                            'sequence': motif_sequence,
                            'average_conservation': avg_score,
                            'type': self._classify_motif_type(motif_sequence)
                        })
                    
                    current_motif_start = None
        
        # 处理末尾的基序
        if current_motif_start is not None:
            motif_length = len(scores) - current_motif_start
            if motif_length >= min_length:
                motif_sequence = consensus[current_motif_start:]
                motif_scores = scores[current_motif_start:]
                avg_score = np.mean(motif_scores)
                
                motifs.append({
                    'start': current_motif_start + 1,
                    'end': len(scores),
                    'length': motif_length,
                    'sequence': motif_sequence,
                    'average_conservation': avg_score,
                    'type': self._classify_motif_type(motif_sequence)
                })
        
        # 按长度和保守性排序
        motifs.sort(key=lambda x: (x['length'], x['average_conservation']), reverse=True)
        
        logger.info(f"发现{len(motifs)}个保守基序")
        return motifs
    
    def _classify_motif_type(self, sequence: str) -> str:
        """分类基序类型"""
        # 简单的基序分类规则
        if 'AAAA' in sequence or 'UUUU' in sequence:
            return 'poly_A/U'
        elif 'CCGG' in sequence or 'GGCC' in sequence:
            return 'GC_rich'
        elif sequence.count('A') + sequence.count('U') > len(sequence) * 0.7:
            return 'AU_rich'
        elif sequence.count('G') + sequence.count('C') > len(sequence) * 0.7:
            return 'GC_rich'
        else:
            return 'mixed'
    
    def map_conservation_to_structure(self, conservation_data: Dict, structure: str) -> Dict:
        """
        将保守性映射到RNA结构上
        Args:
            conservation_data: 保守性数据
            structure: RNA二级结构
        Returns:
            结构-保守性映射结果
        """
        logger.info("映射保守性到RNA结构...")
        
        scores = conservation_data['conservation_scores']
        high_cons_positions = set(conservation_data['high_conservation_positions'])
        
        # 分析结构元件的保守性
        stem_positions = []
        loop_positions = []
        unpaired_positions = []
        
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
                stem_positions.append(i)
            elif char == ')':
                if stack:
                    start = stack.pop()
                    stem_positions.append(i)
            elif char == '.':
                unpaired_positions.append(i)
        
        # 计算各结构元件的平均保守性
        stem_conservation = [scores[i] for i in stem_positions if i < len(scores)]
        loop_conservation = [scores[i] for i in unpaired_positions if i < len(scores)]
        
        # 识别高保守性的结构元件
        highly_conserved_stems = []
        for i in stem_positions:
            if i < len(scores) and i in high_cons_positions:
                highly_conserved_stems.append(i)
        
        highly_conserved_loops = []
        for i in unpaired_positions:
            if i < len(scores) and i in high_cons_positions:
                highly_conserved_loops.append(i)
        
        return {
            'stem_positions': stem_positions,
            'loop_positions': unpaired_positions,
            'stem_average_conservation': np.mean(stem_conservation) if stem_conservation else 0,
            'loop_average_conservation': np.mean(loop_conservation) if loop_conservation else 0,
            'highly_conserved_stems': highly_conserved_stems,
            'highly_conserved_loops': highly_conserved_loops,
            'structure_conservation_correlation': self._calculate_structure_conservation_correlation(
                scores, structure)
        }
    
    def _calculate_structure_conservation_correlation(self, scores: List[float], structure: str) -> float:
        """计算结构与保守性的相关性"""
        if len(scores) != len(structure):
            return 0.0
        
        # 将结构转换为数值（配对=1，未配对=0）
        structure_numeric = [1 if char in '()' else 0 for char in structure]
        
        # 计算相关系数
        if len(structure_numeric) > 1:
            correlation = np.corrcoef(scores, structure_numeric)[0, 1]
            return correlation if not np.isnan(correlation) else 0.0
        
        return 0.0
    
    def analyze_gene_conservation(self, gene_name: str) -> Dict:
        """
        分析单个基因的保守性
        Args:
            gene_name: 基因名称
        Returns:
            保守性分析结果
        """
        logger.info(f"开始分析{gene_name}的保守性...")
        
        # 收集序列
        sequences = self._collect_gene_sequences(gene_name)
        if len(sequences) < 2:
            logger.warning(f"{gene_name}可用序列少于2个，跳过保守性分析")
            return None
        
        # 多序列比对
        alignment = self.perform_multiple_sequence_alignment(sequences, gene_name)
        if not alignment:
            return None
        
        # 序列保守性分析
        sequence_conservation = self.calculate_sequence_conservation(alignment)
        
        # 密码子保守性分析（仅对编码序列）
        codon_conservation = None
        if gene_name in ['RdRp', 'N', 'S', 'M', 'E']:  # 已知编码基因
            codon_conservation = self.analyze_codon_conservation(alignment)
        
        # 识别保守基序
        conserved_motifs = self.identify_conserved_motifs(sequence_conservation)
        
        # 读取结构信息（如果存在）
        structure_conservation = None
        structure_file = Path(f"data/structures/{gene_name}_structure.txt")
        if structure_file.exists():
            with open(structure_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    structure = lines[1].strip()
                    structure_conservation = self.map_conservation_to_structure(
                        sequence_conservation, structure)
        
        # 整合结果
        analysis_result = {
            'gene_name': gene_name,
            'n_sequences': len(sequences),
            'sequence_conservation': sequence_conservation,
            'codon_conservation': codon_conservation,
            'conserved_motifs': conserved_motifs,
            'structure_conservation': structure_conservation,
            'sequences_used': list(sequences.keys())
        }
        
        # 保存结果
        self._save_conservation_results(gene_name, analysis_result)
        
        return analysis_result
    
    def _collect_gene_sequences(self, gene_name: str) -> Dict[str, str]:
        """收集基因的同源序列"""
        sequences = {}
        
        # 读取SARS-CoV-2序列
        sars_cov2_file = Path(f"data/sequences/sars_cov2_{gene_name.lower()}.fasta")
        if sars_cov2_file.exists():
            record = SeqIO.read(sars_cov2_file, "fasta")
            sequences['SARS-CoV-2'] = str(record.seq)
        
        # 读取其他冠状病毒序列（需要根据基因位置提取）
        related_viruses = ['sars', 'mers', 'sars_cov_2_omicron']
        for virus in related_viruses:
            virus_file = Path(f"data/sequences/{virus}_genome.fasta")
            if virus_file.exists():
                record = SeqIO.read(virus_file, "fasta")
                genome_seq = str(record.seq)
                
                # 这里需要根据基因位置提取相应序列
                # 简化处理：假设基因位置相似（实际需要更精确的比对）
                if gene_name == 'RdRp':
                    start, end = 26244, 29533
                elif gene_name == 'N':
                    start, end = 28274, 29533
                else:
                    continue
                
                if len(genome_seq) > end:
                    gene_seq = genome_seq[start-1:end]
                    sequences[virus.upper()] = gene_seq
        
        return sequences
    
    def _save_conservation_results(self, gene_name: str, results: Dict):
        """保存保守性分析结果"""
        # 保存JSON格式结果
        json_file = self.results_dir / f"{gene_name}_conservation.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2, default=str)
        
        # 保存保守基序表格
        if results['conserved_motifs']:
            motifs_df = pd.DataFrame(results['conserved_motifs'])
            motifs_file = self.results_dir / f"{gene_name}_conserved_motifs.csv"
            motifs_df.to_csv(motifs_file, index=False, encoding='utf-8-sig')
        
        # 保存密码子保守性表格
        if results.get('codon_conservation') and results['codon_conservation']['codon_conservation']:
            codon_df = pd.DataFrame(results['codon_conservation']['codon_conservation'])
            codon_file = self.results_dir / f"{gene_name}_codon_conservation.csv"
            codon_df.to_csv(codon_file, index=False, encoding='utf-8-sig')
        
        logger.info(f"{gene_name}保守性分析结果已保存")

def main():
    """主函数"""
    logger.info("开始保守性分析...")
    
    # 创建分析器
    analyzer = ConservationAnalyzer()
    
    # 要分析的基因列表
    genes_to_analyze = ['RdRp', 'N', 'S', 'M', 'E']
    
    all_results = {}
    
    for gene in genes_to_analyze:
        result = analyzer.analyze_gene_conservation(gene)
        if result:
            all_results[gene] = result
    
    # 生成综合保守性报告
    if all_results:
        logger.info(f"完成{len(all_results)}个基因的保守性分析")
        
        # 创建保守性汇总表
        summary_data = []
        for gene, result in all_results.items():
            seq_cons = result['sequence_conservation']
            summary_data.append({
                '基因': gene,
                '序列数量': result['n_sequences'],
                '比对长度': seq_cons['alignment_length'],
                '平均保守性': round(seq_cons['average_conservation'], 3),
                '高保守性位置数': len(seq_cons['high_conservation_positions']),
                '保守基序数': len(result['conserved_motifs']),
                '密码子保守性': round(result.get('codon_conservation', {}).get('average_aa_conservation', 0), 3) if result.get('codon_conservation') else 'N/A'
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = analyzer.results_dir / "conservation_summary.csv"
        summary_df.to_csv(summary_file, index=False, encoding='utf-8-sig')
        
        logger.info(f"保守性分析汇总表已保存: {summary_file}")
    
    logger.info("保守性分析完成！")

if __name__ == "__main__":
    main()