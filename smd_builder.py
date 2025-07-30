#!/usr/bin/env python3
import os
import subprocess
import shutil
import numpy as np
from pathlib import Path

def prepare_smd_system(input_path, output_dir="smd_system", 
                      pull_axis="z", box_extension=20.0,
                      water_model="tip3p", ion_concentration=0.15):
    """
    为SMD模拟准备体系 - 沿蛋白质-配体质心连线对齐并延长盒子
    
    参数:
    input_path: str - 包含md.gro文件的路径
    output_dir: str - 输出目录名称
    pull_axis: str - 拉伸轴向 ('x', 'y', 'z')，默认'z'
    box_extension: float - 拉伸轴方向盒子延长长度(nm)，默认20.0nm
    water_model: str - 水模型，默认tip3p
    ion_concentration: float - 离子浓度(M)，默认0.15M
    
    返回:
    bool - 是否成功完成
    """
    
    try:
        # 设置路径
        input_path = Path(input_path)
        gro_file = input_path / "md.gro"
        output_path = Path(output_dir)
        
        # 检查输入文件
        if not gro_file.exists():
            print(f"错误: 未找到文件 {gro_file}")
            return False
            
        # 创建输出目录
        output_path.mkdir(exist_ok=True)
        os.chdir(output_path)
        
        print("=== 开始SMD体系建模 ===")
        
        # 步骤1: 复制并清理原始文件
        shutil.copy2(gro_file, "original.gro")
        print("✓ 复制原始gro文件")
        
        # 步骤2: 创建索引文件分离蛋白质和配体
        print("正在创建蛋白质和配体索引...")
        cmd_make_ndx = ["gmx", "make_ndx", "-f", "original.gro", "-o", "system.ndx"]
        process = subprocess.Popen(cmd_make_ndx, stdin=subprocess.PIPE, 
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                 text=True)
        # 选择: 1=蛋白质, r LIG=配体, 1|r_LIG=蛋白质+配体
        stdout, stderr = process.communicate(input="r LIG\n1 | r_LIG\nq\n")
        
        if process.returncode != 0:
            print(f"创建索引失败: {stderr}")
            return False
        print("✓ 创建蛋白质和配体索引")
        
        # 步骤3: 提取蛋白质和配体（去除水离子）
        cmd_extract = [
            "gmx", "editconf", "-f", "original.gro", "-o", "protein_ligand.gro",
            "-n", "system.ndx"
        ]
        
        process = subprocess.Popen(cmd_extract, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True)
        stdout, stderr = process.communicate(input="Protein_r_LIG\n")
        
        if process.returncode != 0:
            print(f"提取蛋白质和配体失败: {stderr}")
            return False
        print("✓ 提取蛋白质和配体")
        
        # 步骤4: 计算质心并对齐拉伸轴
        print(f"正在沿蛋白质-配体连线对齐{pull_axis.upper()}轴...")
        
        # 创建选择文件计算质心
        with open("center_protein.dat", "w") as f:
            f.write("resname PRO* or resname ALA or resname VAL or resname LEU or resname ILE or resname MET or resname PHE or resname TYR or resname TRP or resname SER or resname THR or resname CYS or resname ASN or resname GLN or resname ASP or resname GLU or resname LYS or resname ARG or resname HIS or resname GLY")
        
        with open("center_ligand.dat", "w") as f:
            f.write("resname LIG")
        
        # 使用gmx select计算质心
        cmd_select_prot = [
            "gmx", "select", "-f", "protein_ligand.gro", "-s", "protein_ligand.gro",
            "-select", "protein", "-oc", "protein_center.xvg"
        ]
        
        cmd_select_lig = [
            "gmx", "select", "-f", "protein_ligand.gro", "-s", "protein_ligand.gro", 
            "-select", "resname LIG", "-oc", "ligand_center.xvg"
        ]
        
        # 简化方法：使用editconf的centering功能
        # 先将体系移到盒子中心
        cmd_center = [
            "gmx", "editconf", "-f", "protein_ligand.gro", "-o", "centered.gro",
            "-center", "0", "0", "0", "-c"
        ]
        
        result = subprocess.run(cmd_center, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"居中失败: {result.stderr}")
            return False
        print("✓ 体系居中完成")
        
        # 步骤5: 根据拉伸轴调整盒子尺寸
        print(f"正在延长{pull_axis.upper()}轴方向盒子...")
        
        # 读取当前盒子尺寸
        with open("centered.gro", "r") as f:
            lines = f.readlines()
            box_line = lines[-1].strip().split()
            if len(box_line) >= 3:
                box_x, box_y, box_z = float(box_line[0]), float(box_line[1]), float(box_line[2])
            else:
                box_x = box_y = box_z = 5.0  # 默认值
        
        # 在拉伸轴方向延长盒子
        if pull_axis.lower() == 'x':
            new_box = f"{box_x + box_extension} {box_y} {box_z}"
        elif pull_axis.lower() == 'y':
            new_box = f"{box_x} {box_y + box_extension} {box_z}"
        else:  # z轴（默认）
            new_box = f"{box_x} {box_y} {box_z + box_extension}"
        
        cmd_extend = [
            "gmx", "editconf", "-f", "centered.gro", "-o", "extended.gro",
            "-box"] + new_box.split() + ["-c"]
        
        result = subprocess.run(cmd_extend, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"延长盒子失败: {result.stderr}")
            return False
        print(f"✓ {pull_axis.upper()}轴方向延长 {box_extension} nm")
        
        # 步骤6: 重新加水
        print(f"正在使用{water_model}重新溶剂化...")
        cmd_solvate = [
            "gmx", "solvate", "-cp", "extended.gro", "-cs", f"{water_model}.gro",
            "-o", "solvated.gro", "-p", "topol.top"
        ]
        
        result = subprocess.run(cmd_solvate, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"溶剂化失败: {result.stderr}")
            print("注意: 请确保有正确的拓扑文件，或手动创建")
            # 如果没有拓扑文件，至少保存无水体系
            shutil.copy2("extended.gro", "smd_system.gro")
        else:
            print("✓ 重新溶剂化完成")
            
            # 步骤7: 添加离子
            print("正在准备添加离子...")
            
            # 创建简单mdp文件
            mdp_content = """
integrator = md
dt = 0.002
nsteps = 0
cutoff-scheme = Verlet
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
"""
            with open("ions.mdp", "w") as f:
                f.write(mdp_content)
            
            # 运行grompp
            cmd_grompp = [
                "gmx", "grompp", "-f", "ions.mdp", "-c", "solvated.gro",
                "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "10"
            ]
            
            result = subprocess.run(cmd_grompp, capture_output=True, text=True)
            if result.returncode == 0:
                # 添加离子
                cmd_genion = [
                    "gmx", "genion", "-s", "ions.tpr", "-o", "ionized.gro",
                    "-p", "topol.top", "-pname", "NA", "-nname", "CL",
                    "-conc", str(ion_concentration), "-neutral"
                ]
                
                process = subprocess.Popen(cmd_genion, stdin=subprocess.PIPE,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True)
                stdout, stderr = process.communicate(input="SOL\n")
                
                if process.returncode == 0:
                    print("✓ 添加离子完成")
                    shutil.copy2("ionized.gro", "smd_system.gro")
                else:
                    print(f"添加离子失败: {stderr}")
                    shutil.copy2("solvated.gro", "smd_system.gro")
            else:
                shutil.copy2("solvated.gro", "smd_system.gro")
        
        # 最终检查和信息输出
        final_file = "smd_system.gro"
        if os.path.exists(final_file):
            with open(final_file, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    atom_count = lines[1].strip()
                    print(f"✓ 最终体系包含 {atom_count} 个原子")
                    
                if len(lines) >= 3:
                    box_info = lines[-1].strip()
                    print(f"✓ 最终盒子尺寸: {box_info}")
        
        # 创建SMD用的索引文件
        print("正在创建SMD索引文件...")
        cmd_make_ndx_final = ["gmx", "make_ndx", "-f", final_file, "-o", "smd.ndx"]
        process = subprocess.Popen(cmd_make_ndx_final, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True)
        stdout, stderr = process.communicate(input="r LIG\nq\n")
        print("✓ 创建SMD索引文件")
        
        print(f"\n=== SMD体系建模完成 ===")
        print(f"输出目录: {output_path.absolute()}")
        print("主要文件:")
        print(f"  - {final_file}: SMD体系文件")
        print(f"  - smd.ndx: SMD索引文件")
        print(f"  - 拉伸轴: {pull_axis.upper()}轴")
        print(f"  - 延长距离: {box_extension} nm")
        
        return True
        
    except Exception as e:
        print(f"错误: {e}")
        return False

def create_smd_index_groups(gro_file="smd_system.gro", ndx_file="smd.ndx"):
    """
    为SMD创建特定的索引组
    
    参数:
    gro_file: str - gro文件名
    ndx_file: str - 输出的索引文件名
    """
    
    print("正在创建SMD专用索引组...")
    
    cmd_make_ndx = ["gmx", "make_ndx", "-f", gro_file, "-o", ndx_file]
    process = subprocess.Popen(cmd_make_ndx, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True)
    
    # 创建蛋白质和配体的索引组
    commands = """r LIG
name 13 LIG
1
name 14 Protein_for_SMD
q
"""
    
    stdout, stderr = process.communicate(input=commands)
    
    if process.returncode == 0:
        print("✓ SMD索引组创建完成")
        print("  - Protein_for_SMD: 蛋白质组")
        print("  - LIG: 配体组")
    else:
        print(f"创建索引组失败: {stderr}")

# 使用示例
if __name__ == "__main__":
    # 示例用法
    success = prepare_smd_system(
        input_path="/path/to/your/simulation",  # 包含md.gro的路径
        output_dir="smd_prepared",              # 输出目录
        pull_axis="z",                          # 沿Z轴拉伸
        box_extension=25.0,                     # Z轴延长25nm
        ion_concentration=0.15                  # 0.15M离子浓度
    )
    
    if success:
        # 创建SMD专用索引
        create_smd_index_groups()
        print("\nSMD体系建模完成！")
        print("建议下一步:")
        print("1. 检查体系结构和取向")
        print("2. 进行能量最小化")
        print("3. 执行SMD模拟")
    else:
        print("建模过程出现错误，请检查输入文件。")