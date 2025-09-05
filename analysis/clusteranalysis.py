import numpy as np
import pandas as pd
import fnmatch
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter

class ClusterAnalysis:
    """Analyzes molecular cluster properties and statistics."""
    
    def __init__(self,datafile):
        self.datafile = datafile
        
    def get_cluster_statistics(self, summarys: list, cluster_time) -> pd.DataFrame:
        """
        Calculate descriptive statistics for cluster sizes per frame.
        
        Args:
            cluster_time_array: list of lists or 2D array, cluster sizes per frame
            total_number: total number of particles/nodes
            cluster_time: list/array of frame times
        
        Returns:
            DataFrame with descriptive statistics per frame
        """
        statis_cluster = []
        
        for cluster_sizes, time in zip(summarys, cluster_time):
            sizes = list(cluster_sizes)

            # print(f"Total number of molecular in clusters is : {np.sum(sizes)}")
            size_series = pd.Series(sizes, name='cluster_size')
            stats = size_series.describe().to_dict()
            stats['time'] = time
            statis_cluster.append(stats)
        
        return pd.DataFrame(statis_cluster)
    
    def cal_Csize_probability(self, summarys: list):

        flat_list = [x for row in summarys for x in row]
        counts = Counter(flat_list)
        n_frame = len(summarys)

        values = np.array(list(counts.keys()), dtype=float)
        freqs = np.array(list(counts.values()), dtype=float)

        f_counts = freqs / n_frame
        p_counts = values * f_counts / sum(summarys[0])

        sort_idx = np.argsort(values)
        values_sorted = values[sort_idx]
        p_counts_sorted = p_counts[sort_idx]

        return (values_sorted, p_counts_sorted)
    
    def get_components_ratio(self,components,cluster_summary,res_counts):
        
        # ====== Part 1: 簇内组分比例 ======
        all_ratios = []

        for clusters in cluster_summary:
            cluster_resnames = []
            for cluster in clusters:
                for res in cluster:
                    cluster_resnames.append(res.resname)

            total = len(cluster_resnames)
            if total == 0:
                ratios = {p: 0 for p in components}
            else:
                ratios = {
                    p: sum(1 for name in cluster_resnames if fnmatch.fnmatch(name, p)) / total
                    for p in components
                }
            all_ratios.append(ratios)

        df_clusters = pd.DataFrame(all_ratios)

        stats_clusters = pd.DataFrame({
            'mean': df_clusters.mean(),
            'std': df_clusters.std()
        })

        # ====== Part 2: 系统整体参与率（逐帧） ======
        frame_ratios = []
        for clusters in cluster_summary:  # 每一帧的聚类结果
            clustered_resnames = [res.resname for cluster in clusters for res in cluster]
            clustered_counts = Counter(clustered_resnames)

            ratios = {}
            for pattern in components:
                total_count = sum(count for resname, count in res_counts.items()
                                if fnmatch.fnmatch(resname, pattern))
                cluster_count = sum(count for resname, count in clustered_counts.items()
                                    if fnmatch.fnmatch(resname, pattern))
                ratio = cluster_count / total_count if total_count > 0 else 0
                ratios[pattern] = ratio
            frame_ratios.append(ratios)

        df_participation_frames = pd.DataFrame(frame_ratios)

        stats_participation = pd.DataFrame({
            'mean': df_participation_frames.mean(),
            'std': df_participation_frames.std()
        })
        return stats_clusters, stats_participation
            
    def unwrap_positions(self, positions, box):
        """
        基于第一个原子作为参考点展开原子坐标，确保坐标连续性。
        
        Parameters
        ----------
        positions : np.ndarray
            原子坐标数组，形状为 (N, 3)，可能包含周期性跳跃
        box : np.ndarray
            盒子尺寸数组，形状为 (3,)，表示 [Lx, Ly, Lz]
            
        Returns
        -------
        np.ndarray
            展开后的原子坐标，形状为 (N, 3)
        """
        box = np.asarray(box)
        unwrapped = np.empty_like(positions)
        unwrapped[0] = positions[0].copy()
        
        for i in range(1, len(positions)):
            delta = positions[i] - positions[i-1]
            delta_corrected = delta - np.round(delta / box) * box
            unwrapped[i] = unwrapped[i-1] + delta_corrected
        
        return unwrapped

    def gyration_tensor(self, cluster_set, box):
        S_set = []
        Rgs = []
        principal_axes_set = []
        anisotropy_set = []
        axis_ratios_set = []
        Volume_set = []
        density_set = []
        box = box * 10
        for res_set in cluster_set:
            positions = np.vstack([res.atoms.positions for res in res_set])
            masses = np.concatenate([res.atoms.masses for res in res_set])
            M = np.sum(masses)

            positions = self.unwrap_positions(positions, box)
            r_cm = np.sum(positions * masses[:, None], axis=0) / M
            delta = positions - r_cm

            S = (delta.T * masses) @ delta / M
            Rg = np.sqrt(np.trace(S))
            eigvals = np.linalg.eigvalsh(S)
            principal_axes = np.sqrt(eigvals)  # [短, 中, 长]

            # 各向异性参数 kappa^2
            l1, l2, l3 = eigvals
            kappa2 = 1-3*(l1*l2 + l2*l3 + l3*l1)/(l1+l2+l3)**2

            # 长宽比
            S_len, M_len, L_len = principal_axes
            ratios = [L_len/M_len, M_len/S_len, L_len/S_len]

            #团簇密度
            Volume = 4/3 * np.pi * 5**(3/2) * np.sqrt(l1*l2*l3)
            density = M / Volume

            S_set.append(S)
            Rgs.append(Rg)
            principal_axes_set.append(principal_axes)
            anisotropy_set.append(kappa2)
            axis_ratios_set.append(ratios)
            Volume_set.append(Volume)
            density_set.append(density)

        return S_set, Rgs, principal_axes_set, anisotropy_set, axis_ratios_set,Volume_set, density_set


    # def save_outcome(summarys, time, data_file):

    #     lines_list = []
    #     header = f"{'times':<10} {'count':<10} {'mean':<10} {'std':<10} {'min':<10} {'25%':<10} {'50%':<10} {'75%':<10} {'max':<10}\n"
    #     for j in range(len(summarys)):
    #         line = f"{time[j]:<10} {summarys[j]['count']:<10} {summarys[j]['mean']:<10.2f} {summarys[j]['std']:<10.2f} " \
    #                f"{summarys[j]['min']:<10} {summarys[j]['25%']:<10} {summarys[j]['50%']:<10} {summarys[j]['75%']:<10} " \
    #                f"{summarys[j]['max']:<10}\n"
    #         print(line)
    #         lines_list.append(line)

    #     with open(data_file, "w") as f:
    #         f.writelines(lines_list)

    # def plot_MaxSize_and_AverSize():
