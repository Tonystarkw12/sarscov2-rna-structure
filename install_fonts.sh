#!/bin/bash
# 安装中文字体以修复matplotlib中文显示问题

echo "==================================================="
echo "安装中文字体以修复图片显示问题"
echo "==================================================="

# 检测系统类型
if [ -f /etc/debian_version ]; then
    # Debian/Ubuntu系统
    echo "检测到 Debian/Ubuntu 系统"
    echo "安装文泉驿字体..."
    sudo apt update
    sudo apt install -y fonts-wqy-microhei fonts-wqy-zenhei fonts-noto-cjk

    # 清除字体缓存
    echo "清除字体缓存..."
    fc-cache -fv

    echo "✓ 字体安装完成!"
    echo ""
    echo "安装的字体包括:"
    echo "  - WenQuanYi Micro Hei (文泉驿微米黑)"
    echo "  - WenQuanYi Zen Hei (文泉驿正黑)"
    echo "  - Noto Sans CJK (思源黑体)"
    echo ""
    echo "请重新运行Python程序以使用新字体"

elif [ -f /etc/redhat-release ]; then
    # RedHat/CentOS系统
    echo "检测到 RedHat/CentOS 系统"
    sudo yum install -y wqy-microhei-fonts wqy-zenhei-fonts google-noto-sans-cjk-fonts
    fc-cache -fv
    echo "✓ 字体安装完成!"

elif [ "$(uname)" == "Darwin" ]; then
    # macOS系统
    echo "检测到 macOS 系统"
    echo "macOS通常已经包含中文字体(PingFang等)"
    echo "如果仍有问题,请手动安装字体"
else
    echo "未识别的系统类型,请手动安装字体"
    echo "推荐安装: WenQuanYi Micro Hei 或 Noto Sans CJK"
fi

echo "==================================================="
