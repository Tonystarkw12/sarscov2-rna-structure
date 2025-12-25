#!/usr/bin/env python3
"""
快速测试脚本
用于验证各个模块的基本功能
"""

import sys
import os
from pathlib import Path
import logging

# 添加src目录到Python路径
sys.path.append(str(Path(__file__).parent / "src"))

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_basic_functionality():
    """测试基本功能"""
    logger.info("开始基本功能测试...")
    
    try:
        # 测试导入
        logger.info("测试模块导入...")
        from data_download import SequenceDownloader
        from structure_prediction import RNAStructurePredictor
        from visualization import RNAStructureVisualizer
        from conservation_analysis import ConservationAnalyzer
        from protein_binding_prediction import ProteinBindingPredictor
        logger.info("✅ 所有模块导入成功")
        
        # 测试目录创建
        logger.info("测试目录结构...")
        directories = [
            "data/sequences",
            "data/structures", 
            "data/alignments",
            "results/figures",
            "results/tables",
            "results/reports"
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
        
        logger.info("✅ 目录结构创建成功")
        
        # 创建示例序列文件用于测试
        logger.info("创建示例序列...")
        test_sequence = "AUGGCUAAUCGUAGCUAUCGUAGCUAGCUAGCUAGCUAGCUA"
        
        sequence_file = Path("data/sequences/test_gene.fasta")
        with open(sequence_file, 'w') as f:
            f.write(">test_gene Test sequence for SARS-CoV-2 analysis\n")
            f.write(f"{test_sequence}\n")
        
        logger.info("✅ 示例序列创建成功")
        
        # 测试结构预测（简化版）
        logger.info("测试结构预测...")
        predictor = RNAStructurePredictor()
        
        # 模拟预测结果
        mock_structure = "(((...)))...(((...)))"
        features = predictor.analyze_structural_features(test_sequence, mock_structure)
        
        logger.info(f"✅ 结构特征分析完成: {features.get('stem_count')}个茎区, {features.get('base_pairs')}个碱基对")
        
        # 测试可视化
        logger.info("测试可视化...")
        visualizer = RNAStructureVisualizer()
        visualizer.plot_structure_simple("TestGene", mock_structure, test_sequence)
        visualizer.plot_structure_features(features, "TestGene")
        
        logger.info("✅ 基本可视化测试完成")
        
        # 清理测试文件
        if sequence_file.exists():
            sequence_file.unlink()
        
        logger.info("🎉 所有基本功能测试通过！")
        return True
        
    except Exception as e:
        logger.error(f"❌ 测试失败: {e}")
        return False

def check_dependencies():
    """检查依赖"""
    logger.info("检查系统依赖...")
    
    required_packages = [
        'Bio', 'pandas', 'numpy', 'matplotlib',
        'seaborn', 'plotly', 'tqdm'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            logger.info(f"✅ {package}")
        except ImportError:
            missing_packages.append(package)
            logger.warning(f"❌ {package} 未安装")
    
    if missing_packages:
        logger.warning(f"缺少依赖包: {', '.join(missing_packages)}")
        logger.info("请运行: pip install -r requirements.txt")
        return False
    
    logger.info("✅ 所有Python依赖检查通过")
    
    # 检查系统工具
    import shutil
    
    system_tools = {
        'RNAfold': 'RNAfold',
        'RNAplfold': 'RNAplfold', 
        'mafft': 'mafft'
    }
    
    missing_tools = []
    
    for name, command in system_tools.items():
        if shutil.which(command):
            logger.info(f"✅ {name}")
        else:
            missing_tools.append(name)
            logger.warning(f"❌ {name} 未找到")
    
    if missing_tools:
        logger.warning(f"缺少系统工具: {', '.join(missing_tools)}")
        logger.info("某些功能可能受限，但基本分析仍可进行")
    
    return len(missing_tools) == 0

def main():
    """主函数"""
    logger.info("SARS-CoV-2 RNA结构分析系统 - 快速测试")
    logger.info("="*50)
    
    # 检查依赖
    deps_ok = check_dependencies()
    
    print()
    
    # 基本功能测试
    func_ok = test_basic_functionality()
    
    print()
    logger.info("="*50)
    
    if deps_ok and func_ok:
        logger.info("🎉 所有测试通过！系统准备就绪。")
        logger.info("运行完整分析: python src/main.py")
        return 0
    else:
        logger.warning("⚠️  部分测试未通过，请检查依赖安装。")
        if not deps_ok:
            logger.info("安装依赖: pip install -r requirements.txt")
        return 1

if __name__ == "__main__":
    sys.exit(main())