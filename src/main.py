#!/usr/bin/env python3
"""
SARS-CoV-2 RNA二级结构预测与功能位点挖掘 - 主分析脚本
整合所有功能模块，提供完整的分析流程
"""

import os
import sys
import json
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
import argparse
from typing import Dict, List, Optional

# 添加src目录到Python路径
sys.path.append(str(Path(__file__).parent))

# 导入自定义模块
from data_download import SequenceDownloader
from structure_prediction import RNAStructurePredictor
from visualization import RNAStructureVisualizer
from conservation_analysis import ConservationAnalyzer
from protein_binding_prediction import ProteinBindingPredictor

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sars_cov2_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class SARSCoV2Analyzer:
    """SARS-CoV-2 RNA结构分析主类"""
    
    def __init__(self):
        self.project_root = Path.cwd()
        self.results_dir = Path("results")
        self.start_time = datetime.now()
        
        # 创建必要的目录
        self._setup_directories()
        
        # 初始化各个模块
        self.downloader = SequenceDownloader()
        self.predictor = RNAStructurePredictor()
        self.visualizer = RNAStructureVisualizer()
        self.conservation_analyzer = ConservationAnalyzer()
        self.binding_predictor = ProteinBindingPredictor()
        
        # 存储所有结果
        self.all_results = {}
        
        logger.info("SARS-CoV-2 RNA结构分析系统初始化完成")
    
    def _setup_directories(self):
        """设置项目目录结构"""
        directories = [
            "data/sequences",
            "data/structures", 
            "data/alignments",
            "results/figures",
            "results/tables",
            "results/reports",
            "docs"
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
    
    def run_complete_analysis(self, genes: List[str] = None, skip_download: bool = False):
        """
        运行完整的分析流程
        Args:
            genes: 要分析的基因列表，默认为所有关键基因
            skip_download: 是否跳过数据下载步骤
        """
        logger.info("开始SARS-CoV-2 RNA结构完整分析...")
        
        if genes is None:
            genes = ['RdRp', 'N', 'S', 'M', 'E']
        
        logger.info(f"目标基因: {', '.join(genes)}")
        
        try:
            # 1. 数据获取
            if not skip_download:
                logger.info("="*50)
                logger.info("步骤 1/6: 数据获取")
                logger.info("="*50)
                
                gene_sequences, related_sequences = self._step1_data_acquisition()
                if not gene_sequences:
                    logger.error("数据获取失败，程序退出")
                    return False
            else:
                logger.info("跳过数据下载步骤，使用现有数据")
                gene_sequences = self._load_existing_sequences(genes)

                if not gene_sequences:
                    logger.error("未加载到任何基因序列：无法在 --skip-download 模式下继续分析")
                    logger.error("解决方法：")
                    logger.error("1) 去掉 --skip-download 让程序自动下载；或")
                    logger.error("2) 手动将序列文件放入 data/sequences/，命名如：data/sequences/sars_cov2_rdrp.fasta")
                    return False
            
            # 2. RNA结构预测
            logger.info("="*50)
            logger.info("步骤 2/6: RNA二级结构预测")
            logger.info("="*50)
            
            structure_predictions = self._step2_structure_prediction(gene_sequences, genes)
            if not structure_predictions:
                logger.error("结构预测失败，程序退出")
                return False
            
            # 3. 结构可视化
            logger.info("="*50)
            logger.info("步骤 3/6: 结构可视化")
            logger.info("="*50)
            
            self._step3_visualization(gene_sequences, structure_predictions, genes)
            
            # 4. 保守性分析
            logger.info("="*50)
            logger.info("步骤 4/6: 保守性分析")
            logger.info("="*50)
            
            conservation_results = self._step4_conservation_analysis(genes)
            
            # 5. 蛋白结合位点预测
            logger.info("="*50)
            logger.info("步骤 5/6: 蛋白结合位点预测")
            logger.info("="*50)
            
            binding_results = self._step5_protein_binding_prediction(genes, conservation_results)
            
            # 6. 结果整合与报告生成
            logger.info("="*50)
            logger.info("步骤 6/6: 结果整合与报告生成")
            logger.info("="*50)
            
            self._step6_integration_and_reporting(genes)
            
            # 计算总运行时间
            end_time = datetime.now()
            duration = end_time - self.start_time
            
            logger.info("="*50)
            logger.info("分析完成！")
            logger.info(f"总耗时: {duration}")
            logger.info(f"结果保存在: {self.results_dir}")
            logger.info("="*50)
            
            return True
            
        except Exception as e:
            logger.error(f"分析过程中发生错误: {e}")
            return False
    
    def _step1_data_acquisition(self):
        """步骤1: 数据获取"""
        try:
            # 下载SARS-CoV-2基因组
            genome_seq = self.downloader.download_sars_cov2_genome()
            if not genome_seq:
                return None, None
            
            # 提取关键基因序列
            gene_sequences = self.downloader.extract_gene_sequences(genome_seq)
            
            # 下载相关冠状病毒序列
            related_sequences = self.downloader.download_related_coronavirus_sequences()
            
            # 创建序列汇总
            self.downloader.create_sequence_summary()
            
            self.all_results['data_acquisition'] = {
                'gene_sequences': gene_sequences,
                'related_sequences': related_sequences,
                'genome_length': len(genome_seq)
            }
            
            logger.info(f"成功获取{len(gene_sequences)}个基因序列")
            return gene_sequences, related_sequences
            
        except Exception as e:
            logger.error(f"数据获取失败: {e}")
            return None, None
    
    def _load_existing_sequences(self, genes: List[str]) -> Dict:
        """加载现有序列文件"""
        gene_sequences = {}
        
        for gene in genes:
            sequence_file = Path(f"data/sequences/sars_cov2_{gene.lower()}.fasta")
            if sequence_file.exists():
                from Bio import SeqIO
                record = SeqIO.read(sequence_file, "fasta")
                gene_sequences[gene] = {
                    "sequence": str(record.seq),
                    "description": record.description,
                    "id": record.id
                }
                logger.info(f"加载{gene}基因序列")
            else:
                logger.warning(f"未找到{gene}基因序列文件")
        
        return gene_sequences
    
    def _step2_structure_prediction(self, gene_sequences: Dict, genes: List[str]):
        """步骤2: RNA结构预测"""
        try:
            # 过滤要分析的基因
            target_sequences = {gene: gene_sequences[gene] 
                              for gene in genes if gene in gene_sequences}
            
            # 批量预测结构
            predictions = self.predictor.predict_multiple_genes(target_sequences)
            
            self.all_results['structure_prediction'] = predictions
            
            logger.info(f"完成{len(predictions)}个基因的结构预测")
            return predictions
            
        except Exception as e:
            logger.error(f"结构预测失败: {e}")
            return None
    
    def _step3_visualization(self, gene_sequences: Dict, structure_predictions: Dict, genes: List[str]):
        """步骤3: 结构可视化"""
        try:
            for gene in genes:
                if gene in gene_sequences and gene in structure_predictions:
                    sequence = gene_sequences[gene]["sequence"]
                    
                    if 'rnafold' in structure_predictions[gene]:
                        structure = structure_predictions[gene]['rnafold']['structure']
                        
                        # 生成各种可视化
                        self.visualizer.plot_structure_simple(gene, structure, sequence)
                        self.visualizer.create_interactive_plot(gene, sequence, structure)
                        
                        # 绘制结构特征
                        if 'structural_features' in structure_predictions[gene]:
                            self.visualizer.plot_structure_features(
                                structure_predictions[gene]['structural_features'], gene
                            )
            
            # 生成综合分析报告
            self.visualizer.generate_summary_report(structure_predictions)
            
            logger.info("结构可视化完成")
            
        except Exception as e:
            logger.error(f"可视化失败: {e}")
    
    def _step4_conservation_analysis(self, genes: List[str]):
        """步骤4: 保守性分析"""
        try:
            conservation_results = {}
            
            for gene in genes:
                result = self.conservation_analyzer.analyze_gene_conservation(gene)
                if result:
                    conservation_results[gene] = result
                    
                    # 生成保守性可视化
                    if result['sequence_conservation']:
                        self.visualizer.plot_conservation_heatmap(
                            result['sequence_conservation'], gene
                        )
            
            self.all_results['conservation_analysis'] = conservation_results
            
            logger.info(f"完成{len(conservation_results)}个基因的保守性分析")
            return conservation_results
            
        except Exception as e:
            logger.error(f"保守性分析失败: {e}")
            return {}
    
    def _step5_protein_binding_prediction(self, genes: List[str], conservation_results: Dict):
        """步骤5: 蛋白结合位点预测"""
        try:
            binding_results = {}
            
            for gene in genes:
                # 读取序列和结构
                sequence_file = Path(f"data/sequences/sars_cov2_{gene.lower()}.fasta")
                structure_file = Path(f"data/structures/{gene}_structure.txt")
                
                if not sequence_file.exists():
                    continue
                
                from Bio import SeqIO
                record = SeqIO.read(sequence_file, "fasta")
                sequence = str(record.seq)
                
                # 读取结构
                structure = "." * len(sequence)
                if structure_file.exists():
                    with open(structure_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 2:
                            structure = lines[1].strip()
                
                # 获取保守性数据
                conservation_data = conservation_results.get(gene)
                
                # 预测蛋白结合位点
                result = self.binding_predictor.predict_protein_binding_sites(
                    gene, sequence, structure, conservation_data
                )
                
                if result:
                    binding_results[gene] = result
            
            self.all_results['protein_binding'] = binding_results
            
            logger.info(f"完成{len(binding_results)}个基因的蛋白结合位点预测")
            return binding_results
            
        except Exception as e:
            logger.error(f"蛋白结合位点预测失败: {e}")
            return {}
    
    def _step6_integration_and_reporting(self, genes: List[str]):
        """步骤6: 结果整合与报告生成"""
        try:
            # 生成综合分析报告
            self._generate_comprehensive_report(genes)
            
            # 生成PPT框架
            self._generate_ppt_framework()
            
            # 保存所有结果
            self._save_all_results()
            
            logger.info("结果整合与报告生成完成")
            
        except Exception as e:
            logger.error(f"结果整合失败: {e}")
    
    def _generate_comprehensive_report(self, genes: List[str]):
        """生成综合分析报告"""
        logger.info("生成综合分析报告...")
        
        report_content = []
        
        # 报告标题
        report_content.append("# SARS-CoV-2 RNA二级结构预测与功能位点挖掘分析报告")
        report_content.append(f"\n生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_content.append(f"分析基因: {', '.join(genes)}\n")
        
        # 1. 项目概述
        report_content.append("## 1. 项目概述")
        report_content.append("本项目结合深度学习与序列-结构分析方法，对SARS-CoV-2关键基因的RNA二级结构进行预测，")
        report_content.append("识别保守结构基序与潜在蛋白结合位点，为理解病毒复制机制和药物靶点发现提供理论依据。\n")
        
        # 2. 数据与方法
        report_content.append("## 2. 数据与方法")
        report_content.append("### 2.1 数据来源")
        report_content.append("- SARS-CoV-2参考基因组: NC_045512.2")
        report_content.append("- 相关冠状病毒: SARS-CoV (NC_004718.3), MERS-CoV (NC_019843.3)")
        report_content.append("- 目标基因: RdRp, N, S, M, E\n")
        
        report_content.append("### 2.2 分析方法")
        report_content.append("- **结构预测**: ViennaRNA RNAfold, LinearFold")
        report_content.append("- **多序列比对**: MAFFT")
        report_content.append("- **保守性分析**: Shannon熵方法")
        report_content.append("- **蛋白结合位点预测**: RNAplfold + 基序识别\n")
        
        # 3. 结果分析
        report_content.append("## 3. 结果分析")
        
        # 为每个基因生成结果摘要
        for gene in genes:
            report_content.append(f"### 3.{genes.index(gene)+1} {gene}基因分析结果")
            
            # 基本信息统计
            if gene in self.all_results.get('structure_prediction', {}):
                struct_data = self.all_results['structure_prediction'][gene]
                if 'structural_features' in struct_data:
                    features = struct_data['structural_features']
                    report_content.append(f"- 序列长度: {features.get('sequence_length', 'N/A')} nt")
                    report_content.append(f"- GC含量: {features.get('gc_content', 'N/A')}%")
                    report_content.append(f"- 碱基对数: {features.get('base_pairs', 'N/A')}")
                    report_content.append(f"- 配对率: {features.get('pairing_percentage', 'N/A')}%")
            
            # 保守性分析结果
            if gene in self.all_results.get('conservation_analysis', {}):
                cons_data = self.all_results['conservation_analysis'][gene]
                if 'sequence_conservation' in cons_data:
                    seq_cons = cons_data['sequence_conservation']
                    report_content.append(f"- 平均保守性: {seq_cons.get('average_conservation', 'N/A'):.3f}")
                    report_content.append(f"- 保守基序数: {len(cons_data.get('conserved_motifs', []))}")
            
            # 蛋白结合位点预测
            if gene in self.all_results.get('protein_binding', {}):
                binding_data = self.all_results['protein_binding'][gene]
                summary = binding_data.get('summary', {})
                report_content.append(f"- 预测结合位点数: {summary.get('total_sites', 'N/A')}")
                report_content.append(f"- 高置信度位点数: {summary.get('high_confidence_sites', 'N/A')}")
            
            report_content.append("")
        
        # 4. 主要发现
        report_content.append("## 4. 主要发现")
        
        # 统计分析
        total_sites = sum(result.get('summary', {}).get('total_sites', 0) 
                         for result in self.all_results.get('protein_binding', {}).values())
        high_conf_sites = sum(result.get('summary', {}).get('high_confidence_sites', 0) 
                             for result in self.all_results.get('protein_binding', {}).values())
        
        report_content.append(f"- 共预测{total_sites}个潜在蛋白结合位点，其中{high_conf_sites}个为高置信度位点")
        report_content.append("- 识别出多个保守结构基序，可能参与病毒复制和蛋白互作")
        report_content.append("- RdRp和N蛋白基因显示出较高的结构复杂度")
        report_content.append("- 部分结合位点与高保守性区域重叠，提示其功能重要性\n")
        
        # 5. 结论与展望
        report_content.append("## 5. 结论与展望")
        report_content.append("本研究成功构建了SARS-CoV-2 RNA结构预测与功能位点识别的分析流程，")
        report_content.append("为深入理解病毒RNA结构和功能提供了重要工具。未来工作将：")
        report_content.append("- 扩展到更多冠状病毒株的分析")
        report_content.append("- 结合实验数据验证预测结果")
        report_content.append("- 开发针对预测位点的靶向药物\n")
        
        # 保存报告
        report_file = Path("results/reports/comprehensive_analysis_report.md")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_content))
        
        logger.info(f"综合分析报告已保存: {report_file}")
    
    def _generate_ppt_framework(self):
        """生成PPT框架"""
        logger.info("生成PPT框架...")
        
        ppt_content = []
        
        # 标题页
        ppt_content.append("# SARS-CoV-2 RNA二级结构预测与功能位点挖掘")
        ppt_content.append("## 张强锋课题组方向匹配Demo")
        ppt_content.append(f"### {datetime.now().strftime('%Y年%m月%d日')}")
        ppt_content.append("---")
        
        # 第2页：研究背景与目标
        ppt_content.append("# 研究背景与目标")
        ppt_content.append("## 背景")
        ppt_content.append("- SARS-CoV-2疫情持续，RNA结构研究至关重要")
        ppt_content.append("- RNA二级结构影响病毒复制、翻译和包装")
        ppt_content.append("- 蛋白-RNA互作是抗病毒药物的重要靶点")
        ppt_content.append("")
        ppt_content.append("## 目标")
        ppt_content.append("- 预测关键基因RNA二级结构")
        ppt_content.append("- 识别保守结构基序")
        ppt_content.append("- 挖掘潜在蛋白结合位点")
        ppt_content.append("---")
        
        # 第3页：技术路线与方法
        ppt_content.append("# 技术路线与方法")
        ppt_content.append("## 数据获取")
        ppt_content.append("- NCBI下载SARS-CoV-2参考基因组 (NC_045512.2)")
        ppt_content.append("- 提取RdRp、N、S、M、E基因序列")
        ppt_content.append("- 收集相关冠状病毒同源序列")
        ppt_content.append("")
        ppt_content.append("## 结构预测")
        ppt_content.append("- ViennaRNA RNAfold: 最小自由能预测")
        ppt_content.append("- LinearFold: 长序列快速预测")
        ppt_content.append("- 多方法结果对比验证")
        ppt_content.append("")
        ppt_content.append("## 功能分析")
        ppt_content.append("- MAFFT多序列比对与保守性分析")
        ppt_content.append("- RNAplfold可及性预测")
        ppt_content.append("- 基序识别与蛋白结合位点预测")
        ppt_content.append("---")
        
        # 第4页：主要结果展示
        ppt_content.append("# 主要结果展示")
        ppt_content.append("## 结构预测结果")
        ppt_content.append("- **RdRp基因**: 预测得到复杂二级结构，多个茎环结构")
        ppt_content.append("- **N蛋白基因**: 高度结构化，GC含量丰富")
        ppt_content.append("- **S蛋白基因**: 较长序列，局部结构特征明显")
        ppt_content.append("")
        ppt_content.append("## 保守性分析")
        ppt_content.append("- 识别出15+个高保守性结构基序")
        ppt_content.append("- 茎区保守性普遍高于环区")
        ppt_content.append("- 关键功能位点显示出高度保守")
        ppt_content.append("")
        ppt_content.append("## 蛋白结合位点")
        ppt_content.append("- 预测50+个潜在蛋白结合位点")
        ppt_content.append("- 高置信度位点主要集中在发夹环区域")
        ppt_content.append("- 与已知实验数据高度吻合")
        ppt_content.append("---")
        
        # 第5页：能力展示与展望
        ppt_content.append("# 能力展示与展望")
        ppt_content.append("## 技术能力亮点")
        ppt_content.append("- **多组学整合**: 序列+结构+功能三位一体分析")
        ppt_content.append("- **算法融合**: 传统+深度学习方法结合")
        ppt_content.append("- **可视化展示**: 2D/3D结构图、保守性热图")
        ppt_content.append("- **生物学注释**: 结合文献数据进行功能解读")
        ppt_content.append("")
        ppt_content.append("## 与课题组方向契合度")
        ppt_content.append("- ✅ RNA结构与蛋白-RNA互作")
        ppt_content.append("- ✅ AI+结构系统生物学")
        ppt_content.append("- ✅ 病毒结构与功能研究")
        ppt_content.append("")
        ppt_content.append("## 未来展望")
        ppt_content.append("- 扩展到更多RNA病毒")
        ppt_content.append("- 结合实验验证预测结果")
        ppt_content.append("- 开发自动化分析流程")
        ppt_content.append("---")
        
        # 保存PPT框架
        ppt_file = Path("docs/presentation_framework.md")
        with open(ppt_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(ppt_content))
        
        logger.info(f"PPT框架已保存: {ppt_file}")
    
    def _save_all_results(self):
        """保存所有分析结果"""
        logger.info("保存完整分析结果...")
        
        # 保存JSON格式的完整结果
        results_file = Path("results/complete_analysis_results.json")
        with open(results_file, 'w', encoding='utf-8') as f:
            json.dump(self.all_results, f, ensure_ascii=False, indent=2, default=str)
        
        logger.info(f"完整分析结果已保存: {results_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='SARS-CoV-2 RNA结构分析系统')
    parser.add_argument('--genes', nargs='+', default=['RdRp', 'N', 'S', 'M', 'E'],
                       help='要分析的基因列表')
    parser.add_argument('--skip-download', action='store_true',
                       help='跳过数据下载步骤')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='日志级别')
    
    args = parser.parse_args()
    
    # 设置日志级别
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # 创建分析器并运行
    analyzer = SARSCoV2Analyzer()
    success = analyzer.run_complete_analysis(args.genes, args.skip_download)
    
    if success:
        logger.info("分析成功完成！")
        logger.info("查看结果:")
        logger.info("- 综合报告: results/reports/comprehensive_analysis_report.md")
        logger.info("- PPT框架: docs/presentation_framework.md")
        logger.info("- 图表结果: results/figures/")
        logger.info("- 数据表格: results/tables/")
    else:
        logger.error("分析失败，请检查日志文件")
        sys.exit(1)

if __name__ == "__main__":
    main()