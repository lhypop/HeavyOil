import numpy as np
import pandas as pd
import fnmatch
import seaborn as sns
from scipy import stats
from pathlib import Path
from collections import Counter

class ClusterAnalysis:
    """Analyzes molecular cluster properties and statistics."""
    
    def __init__(self,traject ,datafile, frame_clusters):
        self.traject = traject
        self.datafile = datafile

        self.summarys = [
                    [int(len(cluster)) for cluster in merged_clusters] 
                    for _, merged_clusters in frame_clusters.items()
                    ]
        self.cluster_summary = [ merged_clusters for _, merged_clusters in frame_clusters.items()]
        self.frames = [frame for frame, _ in frame_clusters.items()]
        self.cluster_time = [frame * self.traject.trajectory.dt for frame in self.frames]
        
    def cal_Cnums_overtime(self: dict) -> pd.DataFrame:
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
        
        for cluster_sizes, time in zip(self.summarys, self.cluster_time):
            sizes = list(cluster_sizes)

            # print(f"Total number of molecular in clusters is : {np.sum(sizes)}")
            size_series = pd.Series(sizes, name='cluster_size')
            cstats = size_series.describe().to_dict()
            cstats['time'] = time
            statis_cluster.append(cstats)

        statis_cluster = pd.DataFrame(statis_cluster)
        return statis_cluster
    
    def cal_Csize_probability(self):

        flat_list = [i for item in self.summarys for i in item]           
        flat_array = np.array(flat_list)

        n_frame = len(self.summarys)
        counts = Counter(flat_array)

        values = np.array(list(counts.keys()), dtype=float)
        freqs = np.array(list(counts.values()), dtype=float)

        f_counts = freqs / n_frame
        p_counts = values * f_counts / sum(self.summarys[0])

        sort_idx = np.argsort(values)
        values_sorted = np.array(values[sort_idx],dtype=float)
        p_counts_sorted = np.array(p_counts[sort_idx],dtype=float)

        return {"size": values_sorted, "p": p_counts_sorted}
    
    def cal_components_ratio(self, components, res_counts):
        # Part 1: 簇内组分比例
        all_ratios = []
        for clusters in self.cluster_summary:
            cluster_resnames = [res.resname for cluster in clusters for res in cluster]
            total = len(cluster_resnames)
            
            ratios = {
                key: sum(1 for name in cluster_resnames if name in pattern_list)
                    / total if total > 0 else 0
                for key, pattern_list in components.items()
            }

            all_ratios.append(ratios)
        df_clusters = pd.DataFrame(all_ratios)

        # Part 2: 系统整体参与率（逐帧）
        frame_ratios = []
        for clusters in self.cluster_summary:
            clustered_counts = Counter(res.resname for cluster in clusters for res in cluster)
            ratios = {}
            for key, pattern_list in components.items():
                total_count = sum(res_counts.get(name, 0) for name in pattern_list)
                cluster_count = sum(clustered_counts.get(name, 0) for name in pattern_list)
                ratios[key] = cluster_count / total_count if total_count > 0 else 0
            frame_ratios.append(ratios)
        df_participation = pd.DataFrame(frame_ratios)

        return df_clusters, df_participation
    

    
    # 工具函数
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

    def cal_shape_describe(self):
        """
        Calculate cluster shape descriptors for each frame and each cluster.
        Returns a dictionary of lists.
        """
        Rg_results_all = []
        principal_axes_all = []
        anisotropy_all = []
        axis_ratios_all = []
        volume_all = []
        density_all = []

        for index, _ in enumerate(self.frames):
            box = self.traject.dimensions[:3] / 10
            cluster_set = self.cluster_summary[index]

            for res_set in cluster_set:
                positions = np.vstack([res.atoms.positions for res in res_set]) /10
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
                kappa2 = 1 - 3 * (l1*l2 + l2*l3 + l3*l1) / (l1 + l2 + l3)**2

                # 长宽比
                S_len, M_len, L_len = principal_axes
                ratios = [L_len/M_len, M_len/S_len, L_len/S_len]

                # 团簇体积（椭球近似）
                Volume = 4/3 * np.pi * 5 **(3/2) * np.prod(principal_axes)
                density = M / Volume

                Rg_results_all.append(Rg)
                principal_axes_all.append(principal_axes)
                anisotropy_all.append(kappa2)
                axis_ratios_all.append(ratios)
                volume_all.append(Volume)
                density_all.append(density)

        bins = 50

        Rg_results_all = np.array(Rg_results_all)
        hist_vals_Rg, edges_Rg = np.histogram(Rg_results_all, bins=bins, density=True)
        center_Rg =  0.5 * (edges_Rg[1:] + edges_Rg[:-1])

        principal_axes_all = np.vstack(principal_axes_all)
        hist_vals_PAS, edges_PAS = np.histogram(principal_axes_all[:,0], bins=bins, density=True)
        center_PAS =  0.5 * (edges_PAS[1:] + edges_PAS[:-1])   
        hist_vals_PAM, edges_PAM = np.histogram(principal_axes_all[:,1], bins=bins, density=True)
        center_PAM =  0.5 * (edges_PAM[1:] + edges_PAM[:-1])   
        hist_vals_PAL, edges_PAL = np.histogram(principal_axes_all[:,2], bins=bins, density=True)
        center_PAL =  0.5 * (edges_PAL[1:] + edges_PAL[:-1])   

        anisotropy_all = np.array(anisotropy_all)
        hist_vals_AS, edges_AS = np.histogram(anisotropy_all, bins=bins, density=True)
        center_AS =  0.5 * (edges_AS[1:] + edges_AS[:-1])   

        axis_ratios_all = np.vstack(axis_ratios_all)
        hist_vals_ALM, edges_ALM = np.histogram(axis_ratios_all[:,0], bins=bins, density=True)
        center_ALM =  0.5 * (edges_ALM[1:] + edges_ALM[:-1])   
        hist_vals_AMS, edges_AMS = np.histogram(axis_ratios_all[:,1], bins=bins, density=True)
        center_AMS =  0.5 * (edges_AMS[1:] + edges_AMS[:-1])   
        hist_vals_ALS, edges_ALS = np.histogram(axis_ratios_all[:,2], bins=bins, density=True)
        center_ALS =  0.5 * (edges_ALS[1:] + edges_ALS[:-1])   

        volume_all = np.array(volume_all)
        hist_vals_VO, edges_VO = np.histogram(volume_all, bins=bins, density=True)
        center_VO =  0.5 * (edges_VO[1:] + edges_VO[:-1])   

        density_all = np.array(density_all)
        hist_vals_DE, edges_DE = np.histogram(density_all, bins=bins, density=True)
        center_DE =  0.5 * (edges_DE[1:] + edges_DE[:-1])  

        shape_features = {
            "Rg"    : {"cen": center_Rg , "P": hist_vals_Rg , "values":Rg_results_all},
            "PAS"   : {"cen": center_PAS, "P": hist_vals_PAS, "values":principal_axes_all[:,0]},
            "PAM"   : {"cen": center_PAM, "P": hist_vals_PAM, "values":principal_axes_all[:,1]},
            "PAL"   : {"cen": center_PAL, "P": hist_vals_PAL, "values":principal_axes_all[:,2]},
            "AS"    : {"cen": center_AS , "P": hist_vals_AS , "values":anisotropy_all},
            "ALM"   : {"cen": center_ALM, "P": hist_vals_ALM, "values":axis_ratios_all[:,0]},
            "AMS"   : {"cen": center_AMS, "P": hist_vals_AMS, "values":axis_ratios_all[:,1]},
            "ALS"   : {"cen": center_ALS, "P": hist_vals_ALS, "values":axis_ratios_all[:,2]},
            "VO"    : {"cen": center_VO , "P": hist_vals_VO , "values":volume_all},
            "DE"    : {"cen": center_DE , "P": hist_vals_DE , "values":density_all},
        }


        return shape_features
    
    

    def cal_distribution_analysis(self, dr=0.05):
        from scipy.spatial.distance import squareform
        from scipy.ndimage import gaussian_filter1d
        import matplotlib.pyplot as plt

        bins_total = None
        g_r_total = None
        N_frames = len(self.frames)
        dis_scale = min([np.min(self.traject.trajectory[frame].dimensions[:3]) for frame in self.frames]) / 2 / 10

        nnd_all = []
        for index, specific_frame in enumerate(self.frames):
            self.traject.trajectory[specific_frame] 
            box_size = self.traject.dimensions[:3] / 10
            cluster_set = self.cluster_summary[index]
            
            bins = np.arange(0, dis_scale + dr, dr)
            edges = bins

            
            cluster_rcm = []
            for cluster in cluster_set:
                positions = np.vstack([res.atoms.positions for res in cluster]) / 10
                masses = np.concatenate([res.atoms.masses for res in cluster])
                M = np.sum(masses)

                positions = self.unwrap_positions(positions, box_size)
                r_cm = np.sum(positions * masses[:, None], axis=0) / M
                cluster_rcm.append(r_cm)

            coords = np.array(cluster_rcm)

            delta = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            delta = np.abs(delta)
            delta = np.where(delta > 0.5 * box_size, box_size - delta, delta)
            dis_matrix = np.sqrt(np.sum(delta**2, axis=-1))

            # 最近邻距离
            dis_matrix[dis_matrix == 0] = np.inf
            nnd_frame = np.min(dis_matrix, axis=1)
            nnd_all.extend(nnd_frame.tolist())

            # RDF
            distances = squareform(dis_matrix, checks=False)
            hist, _ = np.histogram(distances, bins=bins)

            V = np.prod(box_size)
            rho = len(cluster_rcm) / V
            N_clusters = len(cluster_rcm)

            shell_volumes = 4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
            g_r_frame = hist / (N_clusters * rho * shell_volumes)

            if g_r_total is None:
                g_r_total = np.zeros_like(g_r_frame)
                bins_total = edges
            
            g_r_total += g_r_frame

        # --- RDF 部分 ---
        g_r_avg = g_r_total / N_frames
        r = 0.5 * (bins_total[1:] + bins_total[:-1])
        g_r_avg_smooth = gaussian_filter1d(g_r_avg, sigma=1)
        # plt.plot(r, g_r_avg_smooth)
        # plt.xlabel("r")
        # plt.ylabel("g(r)")
        # plt.show()        


        nnd_all = np.array(nnd_all)
        hist_vals_nll, edges_nll = np.histogram(nnd_all, bins=bins, density=True)
        center_nll =  0.5 * (edges_nll[1:] + edges_nll[:-1])  

        rdf_stats = {
            "r" : r,
            "gr":g_r_avg_smooth,
        }

        nnd_stats = {
            "cen": center_nll , "p": hist_vals_nll , "values":nnd_all
        }

        return rdf_stats, nnd_stats

