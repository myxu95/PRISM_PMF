# 修改 tests/test_basic.py
#!/usr/bin/env python3
"""基础功能测试"""

import sys
import os
from pathlib import Path

# 添加项目根目录到Python路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_import():
    """测试模块导入"""
    try:
        from pmf_package import pmf  # ← 使用包导入
        print("✓ PMF模块导入成功")
        print(f"✓ PMF模块路径: {pmf.__file__}")
        return True
    except Exception as e:
        print(f"✗ PMF模块导入失败: {e}")
        print(f"当前工作目录: {os.getcwd()}")
        print(f"Python路径: {sys.path}")
        return False

def test_setup():
    """测试配置设置"""
    from pmf_package import pmf  # ← 更新导入方式
    
    config = {
        # ... 配置保持不变
    }
    
    try:
        pmf.setup(config, "./pmf_output", "./md_results")
        print("✓ PMF配置设置成功")
        
        status = pmf.get_status()
        print(f"✓ 状态查询成功: {status}")
        return True
    except Exception as e:
        print(f"✗ PMF配置设置失败: {e}")
        return False

# ... 其余代码保持不变