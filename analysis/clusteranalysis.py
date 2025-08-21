import numpy as np
import pandas as pd
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

        print("values_sorted:", values_sorted)
        print("p_counts_sorted:", p_counts_sorted)

    def unwrap_positions_center(self, positions, box):
        """
        基于质心展开原子坐标
        positions: (N,3) numpy array
        box: np.array([Lx,Ly,Lz])
        """
        r_cm = np.mean(positions, axis=0)  # 先粗略计算质心
        delta = positions - r_cm
        # 沿每个维度映射到离质心最近的周期盒子
        for dim in range(3):
            delta[:, dim] -= np.round(delta[:, dim] / box[dim]) * box[dim]
        positions_unwrapped = r_cm + delta
        return positions_unwrapped

    def gyration_tensor(self, cluster_set, box):
        S_set = []
        Rgs = []
        principal_axes_set = []
        anisotropy_set = []
        axis_ratios_set = []
        density_set = []

        for res_set in cluster_set:
            positions = np.vstack([res.atoms.positions for res in res_set])
            masses = np.concatenate([res.atoms.masses for res in res_set])
            M = np.sum(masses)

            positions = self.unwrap_positions_center(positions, box)
            r_cm = np.sum(positions * masses[:, None], axis=0) / M
            delta = positions - r_cm

            S = (delta.T * masses) @ delta / M
            Rg = np.sqrt(np.trace(S))
            eigvals = np.linalg.eigvalsh(S)
            principal_axes = np.sqrt(eigvals)  # [短, 中, 长]

            # 各向异性参数 kappa^2
            l1, l2, l3 = eigvals
            print(l1,l2,l3)
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
            density_set.append(density)

        return S_set, Rgs, principal_axes_set, anisotropy_set, axis_ratios_set, density_set

    def plot_cluster_shape_distribution(self,
                                        Rg_list, 
                                        principal_axes_list, 
                                        anisotropy_list, 
                                        axis_ratios_list, 
                                        density_list,
                                        top_prefix):
        """
        绘制簇的回转半径、主轴、各向异性和长宽比概率密度图。

        参数：
            Rg_list: list, 所有簇的回转半径
            principal_axes_list: list of arrays, 每个簇的三条主轴 [S, M, L]
            anisotropy_list: list, 各向异性指标 kappa^2
            axis_ratios_list: list of arrays, 每个簇的长宽比 [L/M, M/S, L/S]
            top_prefix: str or Path, 输出文件前缀
        """
        # 转为数组
        Rg_array = np.array(Rg_list)
        principal_axes_array = np.vstack(principal_axes_list)
        anisotropy_array = np.array(anisotropy_list)
        axis_ratios_array = np.vstack(axis_ratios_list)
        density_array = np.array(density_list)


        # ---------- 回转半径 ----------
        plt.figure(figsize=(6,4))
        sns.kdeplot(Rg_array, fill=True, bw_adjust=0.5)
        plt.xlabel("Radius of Gyration (Rg)")
        plt.ylabel("Probability Density")
        plt.title("Probability Density of Radius of Gyration")
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_Rg_based_S.png")
        plt.close()

        # ---------- 主轴一维概率密度 ----------
        short_axes = principal_axes_array[:,0]
        middle_axes = principal_axes_array[:,1]
        long_axes = principal_axes_array[:,2]

        plt.figure(figsize=(6,4))
        sns.kdeplot(short_axes, fill=True, label='Short axis', bw_adjust=0.5)
        sns.kdeplot(middle_axes, fill=True, label='Middle axis', bw_adjust=0.5)
        sns.kdeplot(long_axes, fill=True, label='Long axis', bw_adjust=0.5)
        plt.xlabel("Principal Axis Length")
        plt.ylabel("Probability Density")
        plt.title("Probability Density of Cluster Principal Axes")
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_Axis_1D.png")
        plt.close()

        # ---------- 主轴二维投影（长轴 vs 中轴） ----------
        # middle_axes 和 long_axes 是你的数据数组
        corr = np.corrcoef(middle_axes, long_axes)[0,1]  # 计算皮尔逊相关系数
        a, b = np.polyfit(middle_axes, long_axes, 1)     # 拟合直线

        plt.figure(figsize=(6,5))
        # 散点图
        plt.scatter(middle_axes, long_axes, color='blue', alpha=0.5, s=20, label='Data Points')

        # 拟合直线
        x_fit = np.array([min(middle_axes), max(middle_axes)])
        y_fit = a * x_fit + b
        plt.plot(x_fit, y_fit, color='red', lw=2, label=f'Fit Line (r={corr:.2f})')

        plt.xlabel("Middle Axis Length")
        plt.ylabel("Long Axis Length")
        plt.title("Scatter plot with Fit Line")
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_Axis_2D.png")
        plt.close()

        # ---------- 各向异性 ----------
        plt.figure(figsize=(6,4))
        sns.kdeplot(anisotropy_array, fill=True, bw_adjust=0.5)
        plt.xlabel("Shape Anisotropy (kappa^2)")
        plt.ylabel("Probability Density")
        plt.title("Probability Density of Cluster Shape Anisotropy")
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_Anisotropy.png")
        plt.close()

        # ---------- 长宽比 L/M, M/S, L/S ----------
        ratio_labels = ['L/M', 'M/S', 'L/S']
        plt.figure(figsize=(6,4))
        for i, label in enumerate(ratio_labels):
            sns.kdeplot(axis_ratios_array[:,i], fill=True, label=label, bw_adjust=0.5)
        plt.xlabel("Axis Ratios")
        plt.ylabel("Probability Density")
        plt.title("Probability Density of Cluster Axis Ratios")
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_AxisRatios.png")
        plt.close()

        # ---------- 团簇密度 ----------
        plt.figure(figsize=(6,4))
        sns.kdeplot(density_array, fill=True, bw_adjust=0.5)
        plt.xlabel("Cluster density")
        plt.ylabel("Probability Density")
        plt.title("Probability Density of Cluster density")
        plt.grid(True)
        plt.savefig(f"{Path(top_prefix).stem}_Cdensity.png")
        plt.close()


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
