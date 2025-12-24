#!/usr/bin/env python3
"""
数据获取模块
用于从NCBI下载SARS-CoV-2及相关冠状病毒的基因组序列
"""

import os
import requests
from Bio import SeqIO
from Bio import Entrez
import pandas as pd
from pathlib import Path
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SequenceDownloader:
    """序列下载器类"""
    
    def __init__(self, email="your_email@example.com"):
        """
        初始化下载器
        Args:
            email: NCBI Entrez API需要的邮箱地址
        """
        self.email = email
        Entrez.email = email
        self.data_dir = Path("data/sequences")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # SARS-CoV-2关键基因坐标 (基于NC_045512.2)
        self.gene_coordinates = {
            "RdRp": {"start": 26244, "end": 29533, "description": "RNA依赖的RNA聚合酶"},
            "N": {"start": 28274, "end": 29533, "description": "核衣壳蛋白"},
            "S": {"start": 21563, "end": 25384, "description": "刺突蛋白"},
            "M": {"start": 26523, "end": 27191, "description": "膜蛋白"},
            "E": {"start": 26245, "end": 26472, "description": "包膜蛋白"}
        }
    
    def download_sars_cov2_genome(self):
        """下载SARS-CoV-2参考基因组(NC_045512.2)"""
        logger.info("开始下载SARS-CoV-2参考基因组...")
        
        try:
            # 从NCBI下载SARS-CoV-2基因组
            handle = Entrez.efetch(
                db="nucleotide",
                id="NC_045512.2",
                rettype="fasta",
                retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            # 保存完整基因组
            genome_file = self.data_dir / "sars_cov2_genome.fasta"
            SeqIO.write(record, genome_file, "fasta")
            logger.info(f"SARS-CoV-2基因组已保存到: {genome_file}")
            
            return str(record.seq)
            
        except Exception as e:
            logger.error(f"下载SARS-CoV-2基因组失败: {e}")
            return None
    
    def extract_gene_sequences(self, genome_seq):
        """从基因组序列中提取关键基因序列"""
        logger.info("提取关键基因序列...")
        
        gene_sequences = {}
        
        for gene_name, coords in self.gene_coordinates.items():
            start, end = coords["start"], coords["end"]
            gene_seq = genome_seq[start-1:end]  # 转换为0-based索引
            
            # 保存基因序列
            gene_record = {
                "id": f"NC_045512.2_{gene_name}",
                "description": f"SARS-CoV-2 {gene_name} gene ({coords['description']})",
                "sequence": str(gene_seq)
            }
            
            gene_file = self.data_dir / f"sars_cov2_{gene_name.lower()}.fasta"
            with open(gene_file, 'w') as f:
                f.write(f">{gene_record['id']} {gene_record['description']}\n")
                f.write(f"{gene_record['sequence']}\n")
            
            gene_sequences[gene_name] = gene_record
            logger.info(f"{gene_name}基因序列已保存，长度: {len(gene_seq)}nt")
        
        return gene_sequences
    
    def download_related_coronavirus_sequences(self):
        """下载相关冠状病毒序列用于保守性分析"""
        logger.info("下载相关冠状病毒序列...")
        
        # 相关冠状病毒的Accession号
        related_viruses = {
            "SARS": "NC_004718.3",  # SARS-CoV
            "MERS": "NC_019843.3",   # MERS-CoV
            "SARS_CoV_2_Omicron": "OM061194.1",  # Omicron变异株
        }
        
        sequences = {}
        
        for virus_name, accession in related_viruses.items():
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="fasta",
                    retmode="text"
                )
                record = SeqIO.read(handle, "fasta")
                handle.close()
                
                # 保存序列
                virus_file = self.data_dir / f"{virus_name.lower()}_genome.fasta"
                SeqIO.write(record, virus_file, "fasta")
                
                sequences[virus_name] = {
                    "accession": accession,
                    "sequence": str(record.seq),
                    "description": record.description
                }
                
                logger.info(f"{virus_name}序列已下载，长度: {len(record.seq)}nt")
                
            except Exception as e:
                logger.error(f"下载{virus_name}序列失败: {e}")
        
        return sequences
    
    def create_sequence_summary(self):
        """创建序列信息汇总表"""
        logger.info("创建序列信息汇总...")
        
        summary_data = []
        
        # 收集所有序列信息
        sequence_files = list(self.data_dir.glob("*.fasta"))
        
        for seq_file in sequence_files:
            try:
                record = SeqIO.read(seq_file, "fasta")
                summary_data.append({
                    "文件名": seq_file.name,
                    "序列ID": record.id,
                    "描述": record.description,
                    "长度(nt)": len(record.seq),
                    "GC含量(%)": round((record.seq.count('G') + record.seq.count('C')) / len(record.seq) * 100, 2)
                })
            except Exception as e:
                logger.warning(f"读取文件{seq_file}失败: {e}")
        
        # 保存汇总表
        summary_df = pd.DataFrame(summary_data)
        summary_file = Path("results/tables/sequence_summary.csv")
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_file, index=False, encoding='utf-8-sig')
        
        logger.info(f"序列信息汇总已保存到: {summary_file}")
        return summary_df

def main():
    """主函数"""
    logger.info("开始数据获取流程...")
    
    # 创建下载器实例
    downloader = SequenceDownloader(email="researcher@example.com")
    
    # 下载SARS-CoV-2基因组
    genome_seq = downloader.download_sars_cov2_genome()
    if not genome_seq:
        logger.error("无法获取SARS-CoV-2基因组，程序退出")
        return
    
    # 提取关键基因序列
    gene_sequences = downloader.extract_gene_sequences(genome_seq)
    
    # 下载相关冠状病毒序列
    related_sequences = downloader.download_related_coronavirus_sequences()
    
    # 创建序列汇总
    summary_df = downloader.create_sequence_summary()
    
    logger.info("数据获取完成！")
    logger.info(f"共提取{len(gene_sequences)}个关键基因序列")
    logger.info(f"共下载{len(related_sequences)}个相关冠状病毒序列")
    
    return gene_sequences, related_sequences

if __name__ == "__main__":
    main()