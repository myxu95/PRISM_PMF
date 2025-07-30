import sys
sys.path.append('./test_scripts')
import pmf

# 配置PMF参数
config = {
    'reference_group': 'Protein',
    'moving_group': 'LIG',
    'smd': {
        'pull_rate': 0.005,
        'pull_k': 1000.0
    },
    'distance': {
        'start': 0.3,
        'end': 2.0
    },
    'umbrella': {
        'production_time_ps': 22000
    },
    'analysis': {
        'energy_unit': 'kCal'
    }
}

# 设置PMF系统
pmf.setup(config, "./pmf_output", "./md_results")

# 顺序执行各模块
print("1. SMD准备...")
prep_results = pmf.smd_preparation()

print("2. SMD模拟...")
smd_results = pmf.smd()

print("3. 伞形采样设置...")
umbrella_results = pmf.umbrella_sampling(auto_run=False)

print("4. PMF分析...")
analysis_results = pmf.pmf_analysis()

print(f"结合能: {analysis_results['binding_energy']['value']:.2f} kcal/mol")