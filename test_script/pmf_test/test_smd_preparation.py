#!/usr/bin/env python3
"""SMD准备模块测试"""

import sys
import os
sys.path.append('./test_scripts')
import pmf

def test_smd_preparation():
    """测试SMD准备功能"""
    print("\n=== 测试SMD准备模块 ===")
    
    # 设置配置
    config = {
        'reference_group': 'Protein',
        'moving_group': 'LIG',
        'smd': {'pull_rate': 0.01, 'pull_k': 1000.0, 'pull_direction': [-1, 0, 0]},
        'distance': {'start': 0.5, 'end': 1.0},
        'simulation': {'dt': 0.002, 'temperature': 300.0}
    }
    
    try:
        # 配置PMF系统
        pmf.setup(config, "./pmf_output", "./md_results")
        print("✓ PMF系统配置完成")
        
        # 运行SMD准备
        results = pmf.smd_preparation()
        print("✓ SMD准备模块执行完成")
        
        # 验证结果
        expected_files = [
            './pmf_output/pmf/solv_ions.gro',
            './pmf_output/pmf/topol.top',
            './pmf_output/mdps/smd.mdp',
            './pmf_output/pmf/create_index.sh',
            './pmf_output/pmf/smd_preparation.sh'
        ]
        
        missing_files = []
        for file_path in expected_files:
            if os.path.exists(file_path):
                print(f"✓ 文件存在: {file_path}")
            else:
                missing_files.append(file_path)
                print(f"✗ 文件缺失: {file_path}")
        
        if not missing_files:
            print("✓ 所有预期文件都已生成")
            return True
        else:
            print(f"✗ 缺失文件: {missing_files}")
            return False
            
    except Exception as e:
        print(f"✗ SMD准备模块测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_smd_preparation()
    if success:
        print("\n✓ SMD准备模块测试通过")
    else:
        print("\n✗ SMD准备模块测试失败")
        sys.exit(1)
