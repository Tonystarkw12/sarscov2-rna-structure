#!/bin/bash

# SARS-CoV-2 RNA结构预测环境配置脚本
# 适用于Ubuntu/Debian系统

echo "开始配置SARS-CoV-2 RNA结构预测环境..."

# 更新包管理器
sudo apt update

# 安装Python和pip
sudo apt install -y python3 python3-pip python3-venv

# 创建虚拟环境
python3 -m venv venv
source venv/bin/activate

# 升级pip
pip install --upgrade pip

# 安装Python依赖
pip install -r requirements.txt

# 安装ViennaRNA
sudo apt install -y vienna-rna

# 安装MAFFT
sudo apt install -y mafft

# 安装其他系统依赖
sudo apt install -y build-essential

# 创建数据目录
mkdir -p data/sequences data/structures results/figures results/tables

echo "环境配置完成！"
echo "激活虚拟环境: source venv/bin/activate"
echo "运行主脚本: python src/main.py"