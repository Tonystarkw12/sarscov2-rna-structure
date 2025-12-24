#!/usr/bin/env python3
"""
RNA结构可视化模块
使用forgi、matplotlib、py3Dmol等工具进行2D/3D结构可视化
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# 尝试导入forgi
try:
    import forgi.visual.mplotlib as fvm
    import forgi.graph.bulge_graph as fgb
    FORGI_AVAILABLE = True
except ImportError:
    FORGI_AVAILABLE = False
    logging.warning("forgi未安装，将使用简化可视化")

# 尝试导入py3Dmol
try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False
    logging.warning("py3Dmol未安装，3D可视化功能受限")

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RNAStructureVisualizer:
    """RNA结构可视化器"""
    
    def __init__(self):
        self.figures_dir = Path("results/figures")
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.structures_dir = Path("data/structures")
        
        # 设置颜色方案
        self.colors = {
            'stem': '#2E86AB',
            'loop': '#A23B72',
            'bulge': '#F18F01',
            'multibranch': '#C73E1D',
            'unpaired': '#6C757D'
        }
    
    def plot_structure_2d_forgi(self, gene_name: str, structure: str, sequence: str) -> bool:
        """
        使用forgi绘制2D RNA结构图
        Args:
            gene_name: 基因名称
            structure: 点括号格式结构
            sequence: RNA序列
        Returns:
            是否成功绘制
        """
        if not FORGI_AVAILABLE:
            logger.warning("forgi不可用，跳过2D结构绘制")
            return False
        
        try:
            logger.info(f"绘制{gene_name}的2D RNA结构...")
            
            # 创建dot-bracket格式文件
            dot_bracket = f">{gene_name}\n{sequence}\n{structure}\n"
            
            # 加载RNA结构
            bg = fgb.BulgeGraph.from_dot_bracket(dot_bracket)
            
            # 创建图形
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # 绘制结构
            fvm.plot_rna(bg, ax=ax, 
                         text_kwargs={"fontsize": 8, "color": "black"},
                         backbone_kwargs={"color": "black", "linewidth": 2})
            
            # 设置标题和样式
            ax.set_title(f"{gene_name} RNA二级结构 (长度: {len(sequence)}nt)", 
                        fontsize=14, fontweight='bold')
            ax.set_aspect('equal')
            
            # 保存图形
            output_file = self.figures_dir / f"{gene_name}_2d_structure.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"2D结构图已保存: {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"绘制{gene_name}的2D结构失败: {e}")
            return False
    
    def plot_structure_simple(self, gene_name: str, structure: str, sequence: str):
        """
        简化的RNA结构可视化（不依赖forgi）
        Args:
            gene_name: 基因名称
            structure: 点括号格式结构
            sequence: RNA序列
        """
        logger.info(f"绘制{gene_name}的简化结构图...")
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
        
        # 上图：序列和结构对应图
        positions = list(range(len(sequence)))
        colors = []
        
        for char in structure:
            if char == '(':
                colors.append(self.colors['stem'])
            elif char == ')':
                colors.append(self.colors['stem'])
            elif char == '.':
                colors.append(self.colors['unpaired'])
            else:
                colors.append(self.colors['loop'])
        
        # 绘制序列-结构对应
        for i, (pos, char, color) in enumerate(zip(positions, sequence, colors)):
            ax1.bar(i, 1, color=color, alpha=0.8, width=1.0)
            if i % 10 == 0:  # 每10个位置标记一个
                ax1.text(i, 0.5, char, ha='center', va='center', fontsize=8)
        
        ax1.set_xlim(-1, len(sequence))
        ax1.set_ylim(0, 1.2)
        ax1.set_title(f"{gene_name} 序列-结构对应图", fontsize=12, fontweight='bold')
        ax1.set_xlabel("核苷酸位置")
        ax1.set_ylabel("结构状态")
        
        # 添加图例
        legend_elements = [
            plt.Rectangle((0,0),1,1, color=self.colors['stem'], label='茎区'),
            plt.Rectangle((0,0),1,1, color=self.colors['unpaired'], label='未配对'),
            plt.Rectangle((0,0),1,1, color=self.colors['loop'], label='环区')
        ]
        ax1.legend(handles=legend_elements, loc='upper right')
        
        # 下图：结构弧线图
        ax2.set_xlim(-1, len(sequence))
        ax2.set_ylim(-1, len(sequence))
        
        # 绘制配对弧线
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')' and stack:
                start = stack.pop()
                # 绘制弧线
                x = np.linspace(start, i, 50)
                y = np.sqrt(max(0, (i-start)**2/4 - (x - (start+i)/2)**2))
                ax2.plot(x, y, color=self.colors['stem'], alpha=0.6, linewidth=1)
        
        # 绘制序列轴线
        ax2.plot([0, len(sequence)], [0, 0], 'k-', linewidth=2)
        ax2.set_title(f"{gene_name} RNA二级结构弧线图", fontsize=12, fontweight='bold')
        ax2.set_xlabel("核苷酸位置")
        ax2.set_ylabel("配对高度")
        
        plt.tight_layout()
        
        # 保存图形
        output_file = self.figures_dir / f"{gene_name}_simple_structure.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"简化结构图已保存: {output_file}")
    
    def plot_conservation_heatmap(self, conservation_data: Dict, gene_name: str):
        """
        绘制保守性热图
        Args:
            conservation_data: 保守性数据
            gene_name: 基因名称
        """
        logger.info(f"绘制{gene_name}的保守性热图...")
        
        # 准备数据
        if isinstance(conservation_data, dict) and 'scores' in conservation_data:
            scores = conservation_data['scores']
        else:
            logger.warning("保守性数据格式不正确")
            return
        
        # 创建热图数据
        window_size = 50  # 滑动窗口大小
        if len(scores) > window_size:
            # 使用滑动窗口平滑
            smoothed_scores = []
            for i in range(len(scores) - window_size + 1):
                window_avg = np.mean(scores[i:i+window_size])
                smoothed_scores.append(window_avg)
            
            positions = list(range(window_size//2, len(scores) - window_size//2))
            plot_scores = smoothed_scores
        else:
            positions = list(range(len(scores)))
            plot_scores = scores
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # 绘制热图
        heatmap_data = np.array(plot_scores).reshape(1, -1)
        im = ax.imshow(heatmap_data, cmap='RdYlBu', aspect='auto', 
                      vmin=0, vmax=1, interpolation='nearest')
        
        # 设置坐标轴
        ax.set_xlabel("核苷酸位置", fontsize=12)
        ax.set_ylabel("保守性", fontsize=12)
        ax.set_title(f"{gene_name} 序列保守性分析", fontsize=14, fontweight='bold')
        
        # 设置x轴刻度
        step = max(1, len(positions) // 10)
        ax.set_xticks(range(0, len(positions), step))
        ax.set_xticklabels([positions[i] for i in range(0, len(positions), step)])
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('保守性分数', rotation=270, labelpad=15)
        
        # 标记高保守性区域
        threshold = np.percentile(plot_scores, 80)  # 前20%为高保守性
        high_conservation = [i for i, score in enumerate(plot_scores) if score > threshold]
        
        if high_conservation:
            for pos in high_conservation[::max(1, len(high_conservation)//5)]:  # 最多标记5个位置
                ax.axvline(x=pos, color='red', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        
        # 保存图形
        output_file = self.figures_dir / f"{gene_name}_conservation_heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"保守性热图已保存: {output_file}")
    
    def plot_structure_features(self, features_data: Dict, gene_name: str):
        """
        绘制结构特征统计图
        Args:
            features_data: 结构特征数据
            gene_name: 基因名称
        """
        logger.info(f"绘制{gene_name}的结构特征图...")
        
        # 创建子图
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. 基本统计信息
        basic_stats = {
            '序列长度': features_data.get('sequence_length', 0),
            '碱基对数': features_data.get('base_pairs', 0),
            '茎区数': features_data.get('stem_count', 0),
            '环区数': features_data.get('loop_count', 0)
        }
        
        bars1 = ax1.bar(basic_stats.keys(), basic_stats.values(), 
                       color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
        ax1.set_title('基本结构统计', fontweight='bold')
        ax1.set_ylabel('数量')
        
        # 添加数值标签
        for bar, value in zip(bars1, basic_stats.values()):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(value), ha='center', va='bottom')
        
        # 2. 配对率和GC含量
        composition_data = {
            '配对率': features_data.get('pairing_percentage', 0),
            'GC含量': features_data.get('gc_content', 0),
            '未配对率': 100 - features_data.get('pairing_percentage', 0)
        }
        
        colors2 = ['#2E86AB', '#A23B72', '#F18F01']
        wedges, texts, autotexts = ax2.pie(composition_data.values(), 
                                          labels=composition_data.keys(),
                                          colors=colors2, autopct='%1.1f%%',
                                          startangle=90)
        ax2.set_title('序列组成分析', fontweight='bold')
        
        # 3. 结构基序分布
        motifs = features_data.get('motifs', [])
        if motifs:
            motif_types = {}
            for motif in motifs:
                motif_type = motif.get('type', 'unknown')
                motif_types[motif_type] = motif_types.get(motif_type, 0) + 1
            
            bars3 = ax3.bar(motif_types.keys(), motif_types.values(),
                           color='#C73E1D')
            ax3.set_title('结构基序分布', fontweight='bold')
            ax3.set_ylabel('数量')
            
            # 添加数值标签
            for bar, value in zip(bars3, motif_types.values()):
                ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                        str(value), ha='center', va='bottom')
        else:
            ax3.text(0.5, 0.5, '未发现明显基序', ha='center', va='center',
                    transform=ax3.transAxes, fontsize=12)
            ax3.set_title('结构基序分布', fontweight='bold')
        
        # 4. 长度分布（如果有多个基因的数据）
        ax4.text(0.5, 0.5, f'序列长度:\n{features_data.get("sequence_length", 0)} nt\n\n'
                          f'结构复杂度:\n{features_data.get("base_pairs", 0)} 个碱基对',
                ha='center', va='center', transform=ax4.transAxes, fontsize=12,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        ax4.set_title('序列信息', fontweight='bold')
        ax4.axis('off')
        
        plt.suptitle(f'{gene_name} RNA结构特征分析', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # 保存图形
        output_file = self.figures_dir / f"{gene_name}_structure_features.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"结构特征图已保存: {output_file}")
    
    def create_interactive_plot(self, gene_name: str, sequence: str, structure: str):
        """
        创建交互式结构图
        Args:
            gene_name: 基因名称
            sequence: RNA序列
            structure: 点括号格式结构
        """
        logger.info(f"创建{gene_name}的交互式结构图...")
        
        # 创建交互式图形
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('序列-结构对应', '配对弧线图'),
            vertical_spacing=0.1
        )
        
        positions = list(range(len(sequence)))
        colors = []
        
        for char in structure:
            if char == '(':
                colors.append(self.colors['stem'])
            elif char == ')':
                colors.append(self.colors['stem'])
            elif char == '.':
                colors.append(self.colors['unpaired'])
            else:
                colors.append(self.colors['loop'])
        
        # 序列-结构对应条形图
        fig.add_trace(
            go.Bar(
                x=positions,
                y=[1] * len(sequence),
                marker_color=colors,
                name='结构状态',
                hovertemplate='位置: %{x}<br>碱基: %{text}<extra></extra>',
                text=list(sequence)
            ),
            row=1, col=1
        )
        
        # 配对弧线图
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')' and stack:
                start = stack.pop()
                x_arc = np.linspace(start, i, 50)
                y_arc = np.sqrt(max(0, (i-start)**2/4 - (x_arc - (start+i)/2)**2))
                
                fig.add_trace(
                    go.Scatter(
                        x=x_arc,
                        y=y_arc,
                        mode='lines',
                        line=dict(color=self.colors['stem'], width=2),
                        name=f'配对 {start}-{i}',
                        hovertemplate=f'配对: {start}-{i}<extra></extra>'
                    ),
                    row=2, col=1
                )
        
        # 更新布局
        fig.update_layout(
            title=f'{gene_name} RNA二级结构 (交互式)',
            height=800,
            showlegend=False
        )
        
        fig.update_xaxes(title_text="核苷酸位置", row=1, col=1)
        fig.update_xaxes(title_text="核苷酸位置", row=2, col=1)
        fig.update_yaxes(title_text="结构状态", row=1, col=1, showticklabels=False)
        fig.update_yaxes(title_text="配对高度", row=2, col=1)
        
        # 保存交互式图形
        output_file = self.figures_dir / f"{gene_name}_interactive_structure.html"
        fig.write_html(output_file)
        
        logger.info(f"交互式结构图已保存: {output_file}")
    
    def generate_summary_report(self, all_predictions: Dict):
        """
        生成综合分析报告图
        Args:
            all_predictions: 所有预测结果
        """
        logger.info("生成综合分析报告...")
        
        # 创建大图
        fig = plt.figure(figsize=(20, 15))
        
        # 收集所有基因的数据
        genes = list(all_predictions.keys())
        
        # 1. 序列长度对比
        ax1 = plt.subplot(3, 3, 1)
        lengths = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                lengths.append(all_predictions[gene]['structural_features'].get('sequence_length', 0))
        
        bars1 = ax1.bar(genes, lengths, color='#1f77b4')
        ax1.set_title('各基因序列长度对比', fontweight='bold')
        ax1.set_ylabel('长度 (nt)')
        plt.xticks(rotation=45)
        
        # 添加数值标签
        for bar, value in zip(bars1, lengths):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(lengths)*0.01,
                    str(value), ha='center', va='bottom')
        
        # 2. GC含量对比
        ax2 = plt.subplot(3, 3, 2)
        gc_contents = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                gc_contents.append(all_predictions[gene]['structural_features'].get('gc_content', 0))
        
        bars2 = ax2.bar(genes, gc_contents, color='#ff7f0e')
        ax2.set_title('各基因GC含量对比', fontweight='bold')
        ax2.set_ylabel('GC含量 (%)')
        plt.xticks(rotation=45)
        
        # 3. 配对率对比
        ax3 = plt.subplot(3, 3, 3)
        pairing_rates = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                pairing_rates.append(all_predictions[gene]['structural_features'].get('pairing_percentage', 0))
        
        bars3 = ax3.bar(genes, pairing_rates, color='#2ca02c')
        ax3.set_title('各基因配对率对比', fontweight='bold')
        ax3.set_ylabel('配对率 (%)')
        plt.xticks(rotation=45)
        
        # 4. 碱基对数量
        ax4 = plt.subplot(3, 3, 4)
        base_pairs = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                base_pairs.append(all_predictions[gene]['structural_features'].get('base_pairs', 0))
        
        bars4 = ax4.bar(genes, base_pairs, color='#d62728')
        ax4.set_title('各基因碱基对数量', fontweight='bold')
        ax4.set_ylabel('碱基对数量')
        plt.xticks(rotation=45)
        
        # 5. 茎区数量
        ax5 = plt.subplot(3, 3, 5)
        stem_counts = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                stem_counts.append(all_predictions[gene]['structural_features'].get('stem_count', 0))
        
        bars5 = ax5.bar(genes, stem_counts, color='#9467bd')
        ax5.set_title('各基因茎区数量', fontweight='bold')
        ax5.set_ylabel('茎区数量')
        plt.xticks(rotation=45)
        
        # 6. 能量对比
        ax6 = plt.subplot(3, 3, 6)
        energies = []
        gene_names_energy = []
        for gene in genes:
            if 'rnafold' in all_predictions[gene] and all_predictions[gene]['rnafold'].get('energy'):
                energies.append(all_predictions[gene]['rnafold']['energy'])
                gene_names_energy.append(gene)
        
        if energies:
            bars6 = ax6.bar(gene_names_energy, energies, color='#8c564b')
            ax6.set_title('各基因预测自由能', fontweight='bold')
            ax6.set_ylabel('自由能 (kcal/mol)')
            plt.xticks(rotation=45)
            ax6.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        # 7. 综合散点图
        ax7 = plt.subplot(3, 3, 7)
        if len(lengths) == len(pairing_rates):
            scatter = ax7.scatter(lengths, pairing_rates, 
                                 s=100, c=gc_contents, 
                                 cmap='viridis', alpha=0.7)
            ax7.set_xlabel('序列长度 (nt)')
            ax7.set_ylabel('配对率 (%)')
            ax7.set_title('序列长度 vs 配对率', fontweight='bold')
            plt.colorbar(scatter, ax=ax7, label='GC含量 (%)')
            
            # 添加基因标签
            for i, gene in enumerate(genes):
                ax7.annotate(gene, (lengths[i], pairing_rates[i]), 
                           xytext=(5, 5), textcoords='offset points')
        
        # 8. 结构复杂度指标
        ax8 = plt.subplot(3, 3, 8)
        complexity_scores = []
        for gene in genes:
            if 'structural_features' in all_predictions[gene]:
                features = all_predictions[gene]['structural_features']
                # 简单的复杂度评分：茎区数 * 配对率 / 100
                complexity = features.get('stem_count', 0) * features.get('pairing_percentage', 0) / 100
                complexity_scores.append(complexity)
        
        bars8 = ax8.bar(genes, complexity_scores, color='#e377c2')
        ax8.set_title('各基因结构复杂度', fontweight='bold')
        ax8.set_ylabel('复杂度评分')
        plt.xticks(rotation=45)
        
        # 9. 总结统计表
        ax9 = plt.subplot(3, 3, 9)
        ax9.axis('off')
        
        # 创建统计表格
        summary_text = "分析总结:\n\n"
        summary_text += f"分析基因数: {len(genes)}\n"
        summary_text += f"平均序列长度: {np.mean(lengths):.1f} nt\n"
        summary_text += f"平均GC含量: {np.mean(gc_contents):.1f}%\n"
        summary_text += f"平均配对率: {np.mean(pairing_rates):.1f}%\n"
        summary_text += f"总碱基对数: {sum(base_pairs)}\n\n"
        summary_text += "主要发现:\n"
        
        # 找出最长序列
        max_len_gene = genes[np.argmax(lengths)]
        summary_text += f"• 最长序列: {max_len_gene} ({max(lengths)} nt)\n"
        
        # 找出最高GC含量
        max_gc_gene = genes[np.argmax(gc_contents)]
        summary_text += f"• 最高GC含量: {max_gc_gene} ({max(gc_contents):.1f}%)\n"
        
        # 找出最高配对率
        max_pairing_gene = genes[np.argmax(pairing_rates)]
        summary_text += f"• 最高配对率: {max_pairing_gene} ({max(pairing_rates):.1f}%)"
        
        ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes, 
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
        
        plt.suptitle('SARS-CoV-2 RNA结构综合分析报告', fontsize=18, fontweight='bold')
        plt.tight_layout()
        
        # 保存综合报告图
        output_file = self.figures_dir / "comprehensive_analysis_report.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"综合分析报告已保存: {output_file}")

def main():
    """主函数"""
    logger.info("开始RNA结构可视化...")
    
    # 创建可视化器
    visualizer = RNAStructureVisualizer()
    
    # 读取预测结果
    structures_dir = Path("data/structures")
    all_predictions = {}
    
    for json_file in structures_dir.glob("*_prediction.json"):
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                prediction = json.load(f)
                gene_name = json_file.stem.replace('_prediction', '')
                all_predictions[gene_name] = prediction
        except Exception as e:
            logger.error(f"读取{json_file}失败: {e}")
    
    # 为每个基因生成可视化
    for gene_name, prediction in all_predictions.items():
        # 读取序列文件
        sequence_file = Path(f"data/sequences/sars_cov2_{gene_name.lower()}.fasta")
        if sequence_file.exists():
            from Bio import SeqIO
            record = SeqIO.read(sequence_file, "fasta")
            sequence = str(record.seq)
            
            # 获取结构信息
            if 'rnafold' in prediction:
                structure = prediction['rnafold']['structure']
                
                # 生成各种可视化
                visualizer.plot_structure_simple(gene_name, structure, sequence)
                
                if FORGI_AVAILABLE:
                    visualizer.plot_structure_2d_forgi(gene_name, structure, sequence)
                
                visualizer.create_interactive_plot(gene_name, sequence, structure)
            
            # 绘制结构特征
            if 'structural_features' in prediction:
                visualizer.plot_structure_features(prediction['structural_features'], gene_name)
    
    # 生成综合报告
    if all_predictions:
        visualizer.generate_summary_report(all_predictions)
    
    logger.info("RNA结构可视化完成！")

if __name__ == "__main__":
    main()
