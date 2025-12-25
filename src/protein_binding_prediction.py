#!/usr/bin/env python3
"""
蛋白结合位点预测模块
使用RNAplfold预测无规卷曲区域，结合已知的蛋白-RNA互作数据预测潜在结合位点
"""

import os
import subprocess
import tempfile
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import json
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ProteinBindingPredictor:
    """蛋白结合位点预测器"""
    
    def __init__(self):
        self.results_dir = Path("results/tables")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir = Path("results/figures")
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # 已知的SARS-CoV-2蛋白-RNA互作数据（简化版本）
        self.known_interactions = {
            'N': {
                'binding_motifs': [
                    {'sequence': 'UGAA', 'description': 'N蛋白识别基序'},
                    {'sequence': 'ACGA', 'description': 'N蛋白结合位点'},
                    {'sequence': 'UUGA', 'description': 'N蛋白高亲和力位点'}
                ],
                'binding_regions': [
                    {'start': 100, 'end': 150, 'description': 'N蛋白N端结构域结合区'},
                    {'start': 250, 'end': 300, 'description': 'N蛋白二聚化结构域结合区'}
                ]
            },
            'RdRp': {
                'binding_motifs': [
                    {'sequence': 'GGUAA', 'description': 'RdRp启动子识别序列'},
                    {'sequence': 'CUAAAC', 'description': 'RdRp模板结合位点'}
                ],
                'binding_regions': [
                    {'start': 50, 'end': 100, 'description': 'RdRp引物结合区'},
                    {'start': 200, 'end': 250, 'description': 'RdRp催化结构域结合区'}
                ]
            }
        }
        
        # RNA结合蛋白的偏好序列
        self.rbp_motifs = {
            'hnRNP': ['UA-rich', 'CA-rich'],
            'SR_protein': ['GA-rich', 'GGAA'],
            'RNA_helicase': ['G-rich', 'polyU'],
            'RIG-I_like': ['5\'-triphosphate', 'dsRNA_pattern'],
            'TIA1': ['U-rich', 'AU-rich']
        }
    
    def calculate_accessibility(self, sequence: str, gene_name: str) -> Dict:
        """
        使用RNAplfold计算RNA可及性
        Args:
            sequence: RNA序列
            gene_name: 基因名称
        Returns:
            可及性分析结果
        """
        logger.info(f"计算{gene_name}的RNA可及性...")
        
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir_path = Path(tmpdir)

                # RNAplfold 读取 stdin；输出文件名通常基于序列 header（如 RdRp_lunp / RdRp_dp.ps）
                fasta_content = f">{gene_name}\n{sequence}\n"

                result = subprocess.run(
                    ['RNAplfold', '-W', '150', '-L', '150', '-u', '30'],
                    input=fasta_content,
                    cwd=tmpdir_path,
                    capture_output=True,
                    text=True,
                    timeout=600
                )

                if result.returncode != 0:
                    logger.error(f"RNAplfold计算失败: {result.stderr}")
                    return None

                lunp_file = tmpdir_path / f"{gene_name}_lunp"
                if not lunp_file.exists():
                    # 兼容不同版本输出
                    lunp_file = tmpdir_path / "plfold_lunp"
                if not lunp_file.exists():
                    logger.warning(f"未找到{gene_name}的无规卷曲概率文件（*_lunp）")
                    return None

                # 解析无规卷曲概率
                # plfold_lunp 可能包含多列（不同片段长度），这里取最后一列作为 -u 30 的结果
                accessibility_scores: List[float] = []
                with open(lunp_file, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        if not line.strip() or line.startswith('#'):
                            continue
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                prob = float(parts[-1])
                            except ValueError:
                                continue
                            accessibility_scores.append(prob)

                if not accessibility_scores:
                    logger.warning(f"{gene_name}无规卷曲概率解析结果为空")
                    return None

                avg_accessibility = float(np.mean(accessibility_scores))
                high_accessibility_threshold = float(np.percentile(accessibility_scores, 80))
                high_accessibility_positions = [
                    i for i, score in enumerate(accessibility_scores)
                    if score >= high_accessibility_threshold
                ]

                return {
                    'accessibility_scores': accessibility_scores,
                    'average_accessibility': avg_accessibility,
                    'high_accessibility_threshold': high_accessibility_threshold,
                    'high_accessibility_positions': high_accessibility_positions,
                    'sequence_length': len(sequence)
                }

        except FileNotFoundError:
            logger.warning("未找到RNAplfold命令，可及性分析将被跳过")
            return None
        except subprocess.TimeoutExpired:
            logger.error(f"RNAplfold计算{gene_name}超时")
            return None
        except Exception as e:
            logger.error(f"计算{gene_name}可及性失败: {e}")
            return None
    
    def predict_motif_binding_sites(self, sequence: str, gene_name: str) -> List[Dict]:
        """
        基于已知基序预测蛋白结合位点
        Args:
            sequence: RNA序列
            gene_name: 基因名称
        Returns:
            预测的结合位点列表
        """
        logger.info(f"基于基序预测{gene_name}的蛋白结合位点...")
        
        binding_sites = []
        
        # 检查已知相互作用基序
        if gene_name in self.known_interactions:
            known_data = self.known_interactions[gene_name]
            
            for motif_info in known_data['binding_motifs']:
                motif_seq = motif_info['sequence']
                description = motif_info['description']
                
                # 在序列中搜索基序
                for i in range(len(sequence) - len(motif_seq) + 1):
                    if sequence[i:i+len(motif_seq)] == motif_seq:
                        binding_sites.append({
                            'start': i + 1,  # 转换为1-based
                            'end': i + len(motif_seq),
                            'sequence': motif_seq,
                            'protein': gene_name,
                            'binding_type': 'known_motif',
                            'description': description,
                            'confidence': 0.9  # 已知基序高置信度
                        })
        
        # 检查RNA结合蛋白偏好基序
        for rbp_name, motif_types in self.rbp_motifs.items():
            for motif_type in motif_types:
                if motif_type == 'UA-rich':
                    # 查找富含UA的区域
                    self._find_rich_regions(sequence, 'U', 'A', rbp_name, binding_sites)
                elif motif_type == 'CA-rich':
                    self._find_rich_regions(sequence, 'C', 'A', rbp_name, binding_sites)
                elif motif_type == 'GA-rich':
                    self._find_rich_regions(sequence, 'G', 'A', rbp_name, binding_sites)
                elif motif_type == 'G-rich':
                    self._find_g_quadruplex(sequence, rbp_name, binding_sites)
                elif motif_type == 'polyU':
                    self._find_poly_regions(sequence, 'U', rbp_name, binding_sites)
                elif motif_type == 'AU-rich':
                    self._find_rich_regions(sequence, 'A', 'U', rbp_name, binding_sites)
        
        # 按置信度排序
        binding_sites.sort(key=lambda x: x['confidence'], reverse=True)
        
        logger.info(f"为{gene_name}预测了{len(binding_sites)}个蛋白结合位点")
        return binding_sites
    
    def _find_rich_regions(self, sequence: str, base1: str, base2: str, protein: str, 
                          binding_sites: List[Dict], min_length: int = 6, threshold: float = 0.7):
        """查找富含特定碱基的区域"""
        window_size = min_length
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            rich_count = window.count(base1) + window.count(base2)
            if rich_count / window_size >= threshold:
                binding_sites.append({
                    'start': i + 1,
                    'end': i + window_size,
                    'sequence': window,
                    'protein': protein,
                    'binding_type': 'rich_region',
                    'description': f'{base1}{base2}-rich region',
                    'confidence': 0.6
                })
    
    def _find_g_quadruplex(self, sequence: str, protein: str, binding_sites: List[Dict]):
        """查找G-四链体基序"""
        # G-四链体模式：G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}
        import re
        pattern = r'G{3,}[ATCU]{1,7}G{3,}[ATCU]{1,7}G{3,}[ATCU]{1,7}G{3,}'
        matches = re.finditer(pattern, sequence)
        
        for match in matches:
            binding_sites.append({
                'start': match.start() + 1,
                'end': match.end(),
                'sequence': match.group(),
                'protein': protein,
                'binding_type': 'g_quadruplex',
                'description': 'G-quadruplex motif',
                'confidence': 0.7
            })
    
    def _find_poly_regions(self, sequence: str, base: str, protein: str, 
                          binding_sites: List[Dict], min_length: int = 4):
        """查找poly区域"""
        pattern = base + f'{{{min_length},}}'
        import re
        matches = re.finditer(pattern, sequence)
        
        for match in matches:
            binding_sites.append({
                'start': match.start() + 1,
                'end': match.end(),
                'sequence': match.group(),
                'protein': protein,
                'binding_type': 'poly_region',
                'description': f'poly-{base} region',
                'confidence': 0.5
            })
    
    def predict_structural_binding_sites(self, sequence: str, structure: str, 
                                       accessibility: Dict) -> List[Dict]:
        """
        基于结构特征预测蛋白结合位点
        Args:
            sequence: RNA序列
            structure: RNA二级结构
            accessibility: 可及性数据
        Returns:
            基于结构的结合位点预测
        """
        logger.info("基于结构特征预测蛋白结合位点...")
        
        if not accessibility:
            logger.warning("缺少可及性数据，跳过结构预测")
            return []
        
        structural_sites = []
        scores = accessibility.get('accessibility_scores', [])
        high_acc_positions = set(accessibility.get('high_accessibility_positions', []))
        
        # 识别茎环结构
        stack = []
        hairpin_loops = []
        internal_loops = []
        
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')' and stack:
                start = stack.pop()
                
                # 检查是否为发夹环
                if i - start > 3:  # 至少4个核苷酸
                    loop_start = start + 1
                    loop_end = i - 1
                    if loop_end > loop_start:
                        hairpin_loops.append({
                            'stem_start': start,
                            'stem_end': i,
                            'loop_start': loop_start,
                            'loop_end': loop_end,
                            'loop_sequence': sequence[loop_start:loop_end+1]
                        })
        
        # 分析发夹环的结合潜力
        for loop in hairpin_loops:
            loop_start = loop['loop_start']
            loop_end = loop['loop_end']
            loop_length = loop_end - loop_start + 1
            
            # 计算环区的平均可及性
            loop_accessibility = []
            for pos in range(loop_start, min(loop_end + 1, len(scores))):
                if pos in high_acc_positions:
                    loop_accessibility.append(scores[pos])
            
            if loop_accessibility:
                avg_acc = np.mean(loop_accessibility)
                
                # 评估结合潜力
                binding_potential = 0
                if 4 <= loop_length <= 8:  # 理想的发夹环大小
                    binding_potential += 0.3
                if avg_acc > 0.5:  # 高可及性
                    binding_potential += 0.4
                if 'U' in loop['loop_sequence']:  # 含尿嘧啶
                    binding_potential += 0.2
                
                if binding_potential > 0.5:
                    structural_sites.append({
                        'start': loop_start + 1,
                        'end': loop_end + 1,
                        'sequence': loop['loop_sequence'],
                        'protein': 'RNA_binding_protein',
                        'binding_type': 'hairpin_loop',
                        'description': f'Hairpin loop (length: {loop_length})',
                        'confidence': min(binding_potential, 1.0),
                        'accessibility': avg_acc
                    })
        
        # 识别无规卷曲区域
        unpaired_regions = []
        current_start = None
        
        for i, char in enumerate(structure):
            if char == '.':
                if current_start is None:
                    current_start = i
            else:
                if current_start is not None:
                    region_length = i - current_start
                    if region_length >= 5:  # 至少5个未配对碱基
                        unpaired_regions.append({
                            'start': current_start,
                            'end': i - 1,
                            'length': region_length,
                            'sequence': sequence[current_start:i]
                        })
                    current_start = None
        
        # 处理末尾的未配对区域
        if current_start is not None:
            region_length = len(structure) - current_start
            if region_length >= 5:
                unpaired_regions.append({
                    'start': current_start,
                    'end': len(structure) - 1,
                    'length': region_length,
                    'sequence': sequence[current_start:]
                })
        
        # 分析未配对区域的结合潜力
        for region in unpaired_regions:
            # 计算区域平均可及性
            region_accessibility = []
            for pos in range(region['start'], min(region['end'] + 1, len(scores))):
                if pos < len(scores):
                    region_accessibility.append(scores[pos])
            
            if region_accessibility:
                avg_acc = np.mean(region_accessibility)
                
                # 评估结合潜力
                binding_potential = 0
                if region['length'] >= 6:  # 足够长的未配对区域
                    binding_potential += 0.3
                if avg_acc > 0.6:  # 高可及性
                    binding_potential += 0.4
                if region['sequence'].count('U') > region['length'] * 0.3:  # 富含U
                    binding_potential += 0.2
                
                if binding_potential > 0.5:
                    structural_sites.append({
                        'start': region['start'] + 1,
                        'end': region['end'] + 1,
                        'sequence': region['sequence'],
                        'protein': 'RNA_binding_protein',
                        'binding_type': 'unpaired_region',
                        'description': f'Unpaired region (length: {region["length"]})',
                        'confidence': min(binding_potential, 1.0),
                        'accessibility': avg_acc
                    })
        
        return structural_sites
    
    def integrate_predictions(self, motif_sites: List[Dict], structural_sites: List[Dict],
                              sequence: str, conservation_data: Optional[Dict] = None) -> List[Dict]:
        """
        整合不同方法的预测结果
        Args:
            motif_sites: 基于基序的预测结果
            structural_sites: 基于结构的预测结果
            conservation_data: 保守性数据
        Returns:
            整合后的预测结果
        """
        logger.info("整合蛋白结合位点预测结果...")
        
        all_sites = motif_sites + structural_sites
        
        # 如果有保守性数据，用于调整置信度（兼容 main.py 传入的整包结构）
        high_positions = None
        if isinstance(conservation_data, dict):
            if 'high_conservation_positions' in conservation_data:
                high_positions = conservation_data.get('high_conservation_positions')
            elif isinstance(conservation_data.get('sequence_conservation'), dict):
                high_positions = conservation_data['sequence_conservation'].get('high_conservation_positions')

        if high_positions:
            high_cons_positions = set(high_positions)

            for site in all_sites:
                site_overlap = 0
                for pos in range(site['start'] - 1, site['end']):
                    if pos in high_cons_positions:
                        site_overlap += 1

                overlap_ratio = site_overlap / (site['end'] - site['start'] + 1)

                if overlap_ratio > 0.5:
                    site['confidence'] = min(site['confidence'] * 1.2, 1.0)
                    site['conservation_boost'] = True
                    site['conservation_overlap'] = overlap_ratio
                else:
                    site['conservation_boost'] = False
                    site['conservation_overlap'] = overlap_ratio

        # 去重：合并重叠的预测位点
        merged_sites = self._merge_overlapping_sites(all_sites, sequence)
        
        # 按置信度排序
        merged_sites.sort(key=lambda x: x['confidence'], reverse=True)
        
        logger.info(f"整合后得到{len(merged_sites)}个蛋白结合位点")
        return merged_sites
    
    def _merge_overlapping_sites(self, sites: List[Dict], sequence: str) -> List[Dict]:
        """合并重叠的结合位点"""
        if not sites:
            return []

        sorted_sites = sorted(sites, key=lambda x: x['start'])
        merged = [sorted_sites[0]]

        for current in sorted_sites[1:]:
            last = merged[-1]

            if current['start'] <= last['end']:
                start = min(last['start'], current['start'])
                end = max(last['end'], current['end'])
                merged_site = {
                    'start': start,
                    'end': end,
                    'sequence': sequence[start - 1:end],
                    'protein': f"{last['protein']}/{current['protein']}",
                    'binding_type': f"{last['binding_type']}+{current['binding_type']}",
                    'description': f"{last['description']} + {current['description']}",
                    'confidence': max(last['confidence'], current['confidence']),
                    'merged': True
                }
                merged[-1] = merged_site
            else:
                merged.append(current)

        return merged
    
    def predict_protein_binding_sites(self, gene_name: str, sequence: str, structure: str, 
                                    conservation_data: Optional[Dict] = None) -> Dict:
        """
        完整的蛋白结合位点预测流程
        Args:
            gene_name: 基因名称
            sequence: RNA序列
            structure: RNA二级结构
            conservation_data: 保守性数据
        Returns:
            蛋白结合位点预测结果
        """
        logger.info(f"开始预测{gene_name}的蛋白结合位点...")
        sequence = sequence.upper().replace('T', 'U')
        
        # 1. 计算可及性
        accessibility = self.calculate_accessibility(sequence, gene_name)
        
        # 2. 基于基序预测
        motif_sites = self.predict_motif_binding_sites(sequence, gene_name)
        
        # 3. 基于结构预测
        structural_sites = []
        if accessibility:
            structural_sites = self.predict_structural_binding_sites(sequence, structure, accessibility)
        
        # 4. 整合预测结果
        integrated_sites = self.integrate_predictions(motif_sites, structural_sites, sequence, conservation_data)
        
        # 5. 生成统计信息
        high_confidence_sites = [site for site in integrated_sites if site['confidence'] >= 0.7]
        binding_types = {}
        proteins = {}
        
        for site in integrated_sites:
            binding_type = site['binding_type']
            protein = site['protein']
            
            binding_types[binding_type] = binding_types.get(binding_type, 0) + 1
            proteins[protein] = proteins.get(protein, 0) + 1
        
        # 整合结果
        prediction_result = {
            'gene_name': gene_name,
            'accessibility': accessibility,
            'motif_based_sites': motif_sites,
            'structure_based_sites': structural_sites,
            'integrated_sites': integrated_sites,
            'summary': {
                'total_sites': len(integrated_sites),
                'high_confidence_sites': len(high_confidence_sites),
                'binding_types': binding_types,
                'target_proteins': proteins,
                'average_confidence': np.mean([site['confidence'] for site in integrated_sites]) if integrated_sites else 0
            }
        }
        
        # 保存结果
        self._save_binding_predictions(gene_name, prediction_result)
        
        return prediction_result
    
    def _save_binding_predictions(self, gene_name: str, results: Dict):
        """保存蛋白结合位点预测结果"""
        # 保存JSON格式结果
        json_file = self.results_dir / f"{gene_name}_protein_binding.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2, default=str)
        
        # 保存结合位点表格
        if results['integrated_sites']:
            sites_df = pd.DataFrame(results['integrated_sites'])
            sites_file = self.results_dir / f"{gene_name}_binding_sites.csv"
            sites_df.to_csv(sites_file, index=False, encoding='utf-8-sig')
        
        logger.info(f"{gene_name}蛋白结合位点预测结果已保存")

def main():
    """主函数"""
    logger.info("开始蛋白结合位点预测...")
    
    # 创建预测器
    predictor = ProteinBindingPredictor()
    
    # 要分析的基因列表
    genes_to_analyze = ['RdRp', 'N', 'S', 'M', 'E']
    
    all_results = {}
    
    for gene in genes_to_analyze:
        # 读取序列和结构
        sequence_file = Path(f"data/sequences/sars_cov2_{gene.lower()}.fasta")
        structure_file = Path(f"data/structures/{gene}_structure.txt")
        
        if not sequence_file.exists():
            logger.warning(f"未找到{gene}的序列文件，跳过")
            continue
        
        # 读取序列
        record = SeqIO.read(sequence_file, "fasta")
        sequence = str(record.seq)
        
        # 读取结构（如果存在）
        structure = "." * len(sequence)  # 默认为未配对
        if structure_file.exists():
            with open(structure_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    structure = lines[1].strip()
        
        # 读取保守性数据（如果存在）
        conservation_data = None
        conservation_file = Path(f"results/tables/{gene}_conservation.json")
        if conservation_file.exists():
            with open(conservation_file, 'r', encoding='utf-8') as f:
                conservation_data = json.load(f)
        
        # 预测蛋白结合位点
        result = predictor.predict_protein_binding_sites(gene, sequence, structure, conservation_data)
        if result:
            all_results[gene] = result
    
    # 生成汇总报告
    if all_results:
        logger.info(f"完成{len(all_results)}个基因的蛋白结合位点预测")
        
        # 创建汇总表
        summary_data = []
        for gene, result in all_results.items():
            summary = result['summary']
            summary_data.append({
                '基因': gene,
                '总结合位点数': summary['total_sites'],
                '高置信度位点数': summary['high_confidence_sites'],
                '平均置信度': round(summary['average_confidence'], 3),
                '主要结合类型': max(summary['binding_types'].items(), key=lambda x: x[1])[0] if summary['binding_types'] else 'N/A',
                '主要靶蛋白': max(summary['target_proteins'].items(), key=lambda x: x[1])[0] if summary['target_proteins'] else 'N/A'
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = predictor.results_dir / "protein_binding_summary.csv"
        summary_df.to_csv(summary_file, index=False, encoding='utf-8-sig')
        
        logger.info(f"蛋白结合位点预测汇总表已保存: {summary_file}")
    
    logger.info("蛋白结合位点预测完成！")

if __name__ == "__main__":
    main()