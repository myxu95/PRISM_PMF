import sys
import os
import shutil
from pathlib import Path
sys.path.append('./test_scripts')
import pmf  # ← 正确的导入位置

def create_mock_smd_results():
    """创建模拟的SMD结果文件"""
    print("创建模拟的SMD结果文件...")
    
    smd_dir = Path("./pmf_output/pmf/smd")
    smd_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建模拟的结果文件
    files_to_create = {
        'smd.gro': '''Test SMD result
     4
    1ALA      N    1   1.500   1.000   1.000
    1ALA     CA    2   1.600   1.000   1.000
    2LIG      C1   3   3.000   2.000   2.000
    2LIG      C2   4   3.100   2.000   2.000
   10.00000  10.00000  10.00000''',
        
        'smd_pullf.xvg': '''# Pull force data
@ title "Pull Force"
@ xaxis label "Time (ps)"
@ yaxis label "Force (kJ/mol/nm)"
0.0    0.0
1.0    10.5
2.0    25.3
3.0    45.1
4.0    32.8
5.0    15.2''',

        'smd_pullx.xvg': '''# Pull distance data
@ title "Pull Distance"
@ xaxis label "Time (ps)"
@ yaxis label "Distance (nm)"
0.0    0.5
1.0    0.6
2.0    0.7
3.0    0.8
4.0    0.9
5.0    1.0'''
    }
    
    for filename, content in files_to_create.items():
        file_path = smd_dir / filename
        with open(file_path, 'w') as f:
            f.write(content)
        print(f"✓ 创建文件: {file_path}")
    
    # 创建轨迹帧目录和文件
    frames_dir = Path("./pmf_output/pmf/trajectory_frames")
    frames_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建几个测试帧
    for i in range(5):
        frame_content = f'''Frame {i}
     4
    1ALA      N    1   {1.5 + i*0.1:.3f}   1.000   1.000
    1ALA     CA    2   {1.6 + i*0.1:.3f}   1.000   1.000
    2LIG      C1   3   {3.0 + i*0.1:.3f}   2.000   2.000
    2LIG      C2   4   {3.1 + i*0.1:.3f}   2.000   2.000
   10.00000  10.00000  10.00000'''
        
        frame_file = frames_dir / f"frame_{i}.gro"
        with open(frame_file, 'w') as f:
            f.write(frame_content)
    
    print(f"✓ 创建了5个轨迹帧文件")

def test_smd_analysis():
    """测试SMD结果分析功能"""
    print("\n=== 测试SMD分析功能 ===")
    
    try:
        # 设置配置
        config = {
            'reference_group': 'Protein',
            'moving_group': 'LIG',
            'smd': {'pull_rate': 0.01, 'pull_k': 1000.0},
            'distance': {'start': 0.5, 'end': 1.0}
        }
        
        pmf.setup(config, "./pmf_output", "./md_results")
        
        # 创建模拟结果
        create_mock_smd_results()
        
        # 测试结果获取功能
        results = pmf._get_smd_results()  # ← 直接使用pmf，不再import
        print("✓ SMD结果获取成功")
        
        # 测试分析功能
        pmf._analyze_smd_results()
        print("✓ SMD结果分析完成")
        
        # 验证分析图表是否生成
        expected_plots = [
            './pmf_output/pmf/smd/force_profile.png',
            './pmf_output/pmf/smd/distance_profile.png'
        ]
        
        for plot_file in expected_plots:
            if os.path.exists(plot_file):
                print(f"✓ 分析图表生成: {plot_file}")
            else:
                print(f"⚠ 分析图表未生成: {plot_file}")
        
        return True
        
    except Exception as e:
        print(f"✗ SMD分析测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_smd_analysis()
    if success:
        print("\n✓ SMD模拟模块测试通过")
    else:
        print("\n✗ SMD模拟模块测试失败")
        sys.exit(1)