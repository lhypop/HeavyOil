from pathlib import Path
import os
import numpy as np
import pandas as pd
import subprocess
from tqdm import tqdm
from collections import Counter
from concurrent.futures import ThreadPoolExecutor


class EnergyAnalysis():
    def __init__(self, topol, trajectoy, dt,clusters, parent_file):
        self.topol = topol
        self.trajectory = trajectoy
        self.dt = dt
        self.clusters = clusters
        self.parent_file = parent_file
        self.target_file_root = self.parent_file/"sub_clusters"
        self.target_file = None

        self.mk_dir(target_file=self.target_file_root)

    def cal_energy(self):
        n_frames = len(self.clusters)
        n_cols = 12
        energy_array = np.empty((n_frames, n_cols))

        for i, (frame, cluster_info) in tqdm(
            enumerate(self.clusters.items()),
            total=n_frames,
            desc="Energy analysis"
        ):
            self.target_file = self.target_file_root / f"frame_{frame}"
            Path(self.target_file).mkdir(exist_ok=True)
            energy_array[i] = self.extract_cluster_trjconv(frame, cluster_info)

        energy_frames = pd.DataFrame(
            energy_array,
            columns=[
                "AllClusterE", "AllClusterV", "AllClusterC", 
                "clusterE", "clusterV", "clusterC",
                "coheE", "coheV", "coheC", 
                "incuE", "incuV", "incuC"
            ]
        )
        energy_frames["frame"] = list(self.clusters.keys())
        
        return energy_frames
    
    def mk_dir(self,target_file):
        if not Path(target_file).exists():
            Path(target_file).mkdir()
            conf_path = Path(target_file)/"conf"
        else:
            def clear_folder(folder: str | Path) -> None:
                folder = Path(folder)
                if not folder.exists() or not folder.is_dir():
                    return
                for item in folder.iterdir():
                    if item.is_file() or item.is_symlink():
                        item.unlink()   # 删除文件或符号链接
                    elif item.is_dir():
                        clear_folder(item)  # 递归清空子文件夹
                        item.rmdir()        # 删除空文件夹
            clear_folder(target_file)
        
    def py_subprocess(self,command,inter_content,parent_file): 

        process = subprocess.Popen(command, 
                                stdin=subprocess.PIPE, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE,
                                cwd= parent_file,
                                text=True)
        if inter_content != None:
            stdout, stderr = process.communicate(input= inter_content)
        else:
            stdout, stderr = process.communicate()
        return stdout, stderr
    
    def read_xvg(self,file_path):
        """
        read data from xvg format
        parameters：
        - file_path: .xvg file path
        return：
        - data: data as numpy array
        """
        data = []
        with open(file_path, 'r') as file:
            for line in file:
                # skip the line strated with # or @
                if line.startswith('#') or line.startswith('@'):
                    continue
                # divide the data
                parts = line.strip().split()
                if parts:  # not none
                    data.append([float(x) for x in parts])
        # transpose inte np.array
        return np.array(data)
    
    def generate_modified_ndx_cluster(self,frame,cluster_info):

        import re
        # 获取已有 group 数量
        def get_existing_group_count(gro_file,ndx_file,parent_file):
            cmd = ['gmx', 'make_ndx', 
                '-f', gro_file, 
                '-o', ndx_file]
            
            stdout, stderr = self.py_subprocess(cmd, 'l\nq\n', parent_file)

            lines = stdout.splitlines()
            group_lines = []
            for line in lines:
                # 使用正则判断行是否以数字开头且只含一个冒号
                if re.match(r'\s*\d+\s+\S.*?:[^:]*$', line):
                    group_lines.append(line)
            return len(group_lines)
        

        # output the structures files
        ndx_file = f"sub_clusters_{frame}_cluster.ndx"
        
        if Path(self.target_file/ndx_file).exists():
            Path(self.target_file/ndx_file).unlink()
        else:
            pass

        n_existing_groups = get_existing_group_count(str(self.topol), self.target_file/ndx_file, self.parent_file)
        Path(self.target_file/ndx_file).unlink()

        make_ndx_command = ['gmx', 'make_ndx', 
                            '-f', str(self.topol), 
                            '-o', self.target_file/ndx_file]
        try:
            input_commands = ""
            research_list = []


            #  residues in all cluster
            all_residues = [res for cluster in cluster_info for res in cluster]
            all_residues = sorted(list(set(all_residues)), key=lambda res: res.resid) ### 重要的

            res_ids = ' '.join(str(res.resid) for res in all_residues)
            input_commands += f"ri {res_ids}\n"
            input_commands += f"name {n_existing_groups + 0} allcluster\n"
            research_list.append("allcluster")

            #  residues in each cluster
            for cluster_i, cluster_residues in enumerate(cluster_info):
                cluster_residues = sorted(list(set(cluster_residues)), key=lambda res: res.resid)
                res_ids = ' '.join(str(res.resid) for res in cluster_residues)
                input_commands += f"ri {res_ids}\n"
                input_commands += f"name {n_existing_groups + 1 + cluster_i} cluster_{cluster_i}\n"
                research_list.append(f"cluster_{cluster_i}")

            #  molecules
            for res_i, res in enumerate(all_residues):
                input_commands += f"ri {res.resid}\n"
                input_commands += f"name {n_existing_groups + 1 + len(cluster_info)+ res_i} res_{res.resid}\n"
                research_list.append(f"res_{res.resid}")
                
            input_commands += "q\n"

            _, _ = self.py_subprocess(make_ndx_command, input_commands,self.parent_file)

            # if stderr:
            #     print(f"Warning from make_ndx:{stderr}")

        except Exception as e:
            print(f"Failed to create cluster .ndx file: {e}")
            raise

        return ndx_file,research_list
    
    def check_file(self,Path,stderr):

        absolute_file = Path
        if absolute_file.exists():
            pass
            # print(f"sucessfully generate {str(Path)} file")  
        else:
            raise stderr
            # print(stderr)
            # print(f"Failed to generate {str(Path)} file") 

    def energy_cal(self,frame, ndx_file, item):

        # xtc extraction of cluster 
        trjconv_command = ['gmx', 'trjconv', 
                '-s', str(self.topol), 
                '-f', str(self.trajectory),
                '-n', str(self.target_file/ndx_file),
                '-o', str(self.target_file/f'conf{frame}_{item}.gro') ,
                '-dump', str(int(frame)*self.dt)]
        input_commands = f"{item}\n"
        stdout, stderr = self.py_subprocess(trjconv_command, input_commands, self.parent_file)
        self.check_file(self.target_file/f'conf{frame}_{item}.gro', stderr)


        
        # tpr extraction of cluster 
        trjconv_command = ['gmx', 'convert-tpr', 
                '-s', str(self.topol), 
                '-n', str(self.target_file/ndx_file),
                '-o', str(self.target_file/f"conf{frame}_{item}.tpr")]
        input_commands = f"{item}\n"
        stdout, stderr= self.py_subprocess(trjconv_command, input_commands, self.parent_file)
        self.check_file(self.target_file/f'conf{frame}_{item}.tpr', stderr)
        # energy rerun of trajectory 
        trjconv_command = ['gmx', 'mdrun', 
                '-v', 
                '-deffnm', f"conf{frame}_{item}", 
                '-rerun', f'conf{frame}_{item}.gro',
                '-e', f"conf{frame}_{item}.edr",
                ]
        input_commands = None
        stdout, stderr = self.py_subprocess(trjconv_command, input_commands, self.target_file)
        self.check_file(self.target_file/f'conf{frame}_{item}.edr', stderr)

        # energy extraction & consider Disper.-corr
        try:
            trjconv_command = ['gmx', 'energy', 
                            '-f', f"conf{frame}_{item}.edr",
                            '-o', f"conf{frame}_{item}.xvg"]
            parameters = ['LJ-(SR)','Disper.-corr.','Coulomb-(SR)','Coul.-recip']
            input_commands = "\n".join(parameters) + "\n"
            stdout, stderr = self.py_subprocess(trjconv_command, input_commands, self.target_file)
        except Exception as e:
            print(f"Warning: Disper.-corr. not found, fallback to simpler energy terms. ({e})")
            trjconv_command = ['gmx', 'energy', 
                            '-f', f"conf{frame}_{item}.edr",
                            '-o', f"conf{frame}_{item}.xvg"]
            parameters = ['LJ-(SR)','Coulomb-(SR)','Coul.-recip']
            input_commands = "\n".join(parameters) + "\n"
            stdout, stderr = self.py_subprocess(trjconv_command, input_commands, self.target_file)
            
        self.check_file(self.target_file/f'conf{frame}_{item}.xvg', stderr)
        data = self.read_xvg(str(self.target_file/f"conf{frame}_{item}.xvg"))


        if len(data) == 5:
            total_cluster_VDW = data[:,1] +data[:,2]
            total_cluster_Coul = data[:,3]+data[:,4]
            total_cluster_edr = data[:,1] +data[:,2]+ data[:,3]+data[:,4]
        else:
            # print("The simulation tpr file don't consider the Disper.-corr")
            total_cluster_VDW = data[:,1] 
            total_cluster_Coul = data[:,2]+data[:,3]
            total_cluster_edr = data[:,1] +data[:,2]+ data[:,3]

        return [item, total_cluster_edr,total_cluster_VDW,total_cluster_Coul]

    def extract_cluster_trjconv(self,frame,cluster_info):
        #load the structure file 

        
        ndx_file,research_list = self.generate_modified_ndx_cluster(frame,cluster_info)

        energy_items = {}

        from functools import partial
        # 把 frame 和 ndx_file 固定住
        func = partial(self.energy_cal, frame, ndx_file)

        with ThreadPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
            results = list(executor.map(func, research_list))

        for item in results:
            energy_items[item[0]] = [item[1], item[2],item[3]]

        # Energy of all cluster
        AllClusterE, AllClusterV, AllClusterC = energy_items["allcluster"]

        # Energy of cohesion == pi-pi interaction
        clusterE, clusterV, clusterC = 0, 0, 0
        coheE, coheV, coheC = 0, 0, 0
        for cluster_i, cluster_residues in enumerate(cluster_info):
            clusteriE, clusteriV, clusteriC = energy_items[f"cluster_{cluster_i}"]
            clusterE += clusteriE
            clusterV += clusteriV
            clusterC += clusteriC

            allresE, allresV, allresC = 0, 0, 0
            for res in cluster_residues:
                resE, resV, resC = energy_items[f"res_{res.resid}"]
                allresE += resE
                allresV += resV
                allresC += resC

            coheiE, coheiV, coheiC = clusteriE - allresE, clusteriV - allresV, clusteriC - allresC
            coheE += coheiE
            coheV += coheiV
            coheC += coheiC

        incuE, incuV, incuC = AllClusterE - clusterE, AllClusterV - clusterV, AllClusterC - clusterC

        # print(AllClusterE, AllClusterV, AllClusterC, clusterE, clusterV, clusterC, coheE, coheV, coheC, incuE, incuV, incuC)
        return      AllClusterE[0], AllClusterV[0], AllClusterC[0], \
                    clusterE[0], clusterV[0], clusterC[0], \
                    coheE[0], coheV[0], coheC[0],\
                    incuE[0], incuV[0], incuC[0]