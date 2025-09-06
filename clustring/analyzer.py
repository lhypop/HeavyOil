import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict
from collections import Counter
from typing import List,Tuple
from tqdm import tqdm

from ..edbscan.FastPeriodicDBSCAN import FastPeriodicDBSCAN
from .coordinate import coordinate
from .postprocess import postprocess
from ..analysis.clusteranalysis import ClusterAnalysis
from ..visual.generate import GMXClusterTool
from ..visual.translate import clusterTranslator
from ..visual.ploter import MyPlot
from ..utils.preprocess import Preprocessor

class ClusterAnalyzer:
    def __init__(self, args: argparse.Namespace):
        """Initialize trajectory analyzer with command line arguments.
        
        Args:
            args: Parsed command line arguments containing:
                - gro: Topology file path (.gro)
                - trj: Trajectory file path (.xtc, optional)
                - residueselect: Residue selection string ('all' or comma-separated names)
                - begin: Start frame index (default: 0)
                - end: End frame index (default: last frame) 
                - interval: Frame sampling interval (default: 1)
                - config: Path to config file (optional)
        """
        # Store raw arguments for reference
        self.args = args
        
        # Required topology file (GRO format)
        self.top = args.gro
        
        # Initialize preprocessing pipeline - handles:
        # 1. Trajectory loading (GRO/XTC)
        # 2. Residue selection parsing
        # 3. Frame range validation
        # 4. Configuration loading
        preprocessor = Preprocessor(args)
        
        # Extract preprocessed parameters
        params = preprocessor.parameters
        self.frames = params['frames']       # Selected frame indices (range object)
        self.residues = params['residues']   # List of residue names to analyze
        self.eps = params['eps']             # DBSCAN epsilon parameter (Å)
        self.min = params['min_samples']  # DBSCAN min_samples parameter
        self.traject = params['universe']    # MDAnalysis Universe object
        
        # Set up temporary directory for intermediate files
        # Default location: ./temp relative to current working directory
        self.temp = Path.cwd() / "temp"
        self.temp.mkdir(exist_ok=True)  # Create if doesn't exist
        
        # Initialize box dimensions - will be set during analysis
        # based on trajectory dimensions
        self.box = None  # Will store [lx, ly, lz] box vectors
    
    def _DBSCAN_clustering(self,
                coords: np.ndarray,
                aromatic_residues: List[str]) -> List[Tuple[int, str]]:
        """
        Perform improved DBSCAN clustering with periodic boundary conditions (PBC).
        
        Parameters
        ----------
        coords : np.ndarray
            Array of coordinates with shape (n_samples, n_features).
        aromatic_residues : List[str]
            Residue names or IDs associated with each point.        
        Returns
        -------
        List[Tuple[int, str]]
            List of (cluster_label, residue) tuples. Noise points are excluded.
        """
        model = FastPeriodicDBSCAN(eps=self.eps, 
                                   min_samples=self.min, 
                                   box_size=self.box)
        cluster_labels = model.fit_predict(coords)

        clustered_residues = [
            (label, res)
            for label, res in zip(cluster_labels, aromatic_residues)
            if label != -1  # -1 means noise
        ]

        return clustered_residues
    
    def analyze_aromatic_clusters(self):

        frame_clusters = defaultdict(list)

        summarys = []
        cluster_summary = []
        cluster_time = []

        coordinator = coordinate(self.temp, self.residues)

        for specific_frame in tqdm(self.frames, desc="Processing"):
            # Get aromatic ring centers and corresponding residues
            self.traject.trajectory[specific_frame] 
            gro_ring_centers, aromatic_residue, self.box = coordinator._get_center_of_aroring(self.traject)
            cluster_labels = self._DBSCAN_clustering(gro_ring_centers, aromatic_residue)

            merged_clusters = postprocess(cluster_labels)._postprocess(self.args.infoselect)
            # Store cluster members for current frame
            for cluster in merged_clusters:
                frame_clusters[specific_frame].append(set(cluster))


            cluster_sizes = [int(len(cluster)) for cluster in merged_clusters]
            summarys.append(cluster_sizes)
            cluster_summary.append(merged_clusters)
            cluster_time.append(self.traject.trajectory.time)


        ######visualization########

        all_cluster_list = frame_clusters.get(self.frames[-1])
        clusters_residues_atoms_position = []

        for _, cluster in enumerate(all_cluster_list):
            cluster_atoms_position = []
            for residue in cluster:
                cluster_atoms_position.append(residue.atoms.positions / 10)
            cluster_atoms_position = np.vstack(cluster_atoms_position)
            clusters_residues_atoms_position.append(cluster_atoms_position)


        # Generate output files
        output_base = Path.cwd() / Path(self.top).stem
        move_tcl_path = output_base.with_name(f"{output_base.name}_move.tcl")
        clusterTranslator(all_cluster_list, self.box, move_tcl_path)
    
        
        processor = GMXClusterTool(
            molname = str(output_base)
        )
        processor.process_clusters(all_cluster_list)

        data_file = Path(f"cluster_nums_{Path(self.top).stem}.txt")
        ClusterInfoAnalyzer = ClusterAnalysis(data_file)

        # # ------数量特征模块------
        # cluster_size = ClusterInfoAnalyzer.cal_Csize_probability(summarys)
        # cluster_describe = ClusterInfoAnalyzer.get_cluster_statistics(summarys, cluster_time)
        # components = ["sa*", "ar*", "re*", "as*"]
        # res_counts = Counter(self.traject.residues.resnames)
        # stats_clusters, stats_participation = ClusterInfoAnalyzer.get_components_ratio(components,cluster_summary,res_counts)

        # # ------形状特征模块------
        # Rg_results_all,principal_axes_all,anisotropy_all,axis_ratios_all,volume_all,density_all = [],[],[],[],[],[]
        # for index, specific_frame in enumerate(self.frames):
        #     self.traject.trajectory[specific_frame] 
        #     self.box = coordinator._get_box_dimension(self.traject)
        #     cluster_set = cluster_summary[index]

        #     S_set, Rgs, principal_axes_set, anisotropy_set, axis_ratios_set, volume_set, density_set = \
        #         ClusterInfoAnalyzer.gyration_tensor(cluster_set, self.box)

        #     Rg_results_all.extend(Rgs)
        #     principal_axes_all.extend(principal_axes_set)
        #     anisotropy_all.extend(anisotropy_set)
        #     axis_ratios_all.extend(axis_ratios_set)
        #     volume_all.extend(volume_set)
        #     density_all.extend(density_set)

        # Rg_array = np.array(Rg_results_all)
        # principal_axes_array = np.vstack(principal_axes_all)
        # anisotropy_array = np.array(anisotropy_all)
        # axis_ratios_array = np.vstack(axis_ratios_all)
        # volume_array = np.array(volume_all)
        # density_array = np.array(density_all)


        # ------分布特征模块------
        import matplotlib.pyplot as plt
        bins_total = None
        g_r_total = None
        N_frames = len(self.frames)
        ## 1. 以最小的盒子大小生成区间

        for index, specific_frame in enumerate(self.frames):
            self.traject.trajectory[specific_frame] 
            self.box = coordinator._get_box_dimension(self.traject)
            cluster_set = cluster_summary[index]

            hist, edges, n_clusters = ClusterInfoAnalyzer.distribution_analysis(cluster_set,self.box)
            # 当前帧体积和密度
            V = np.prod(self.box)
            rho = n_clusters / V

            # RDF归一化
            shell_volumes = 4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
            g_r_frame = hist / (rho * shell_volumes)

            if g_r_total is None:
                g_r_total = np.zeros_like(g_r_frame)
                bins_total = edges

            g_r_total += g_r_frame  # 累积每帧归一化结果

        # 平均
        g_r_avg = g_r_total / N_frames
        r = 0.5 * (bins_total[1:] + bins_total[:-1])
        ## 2. 平滑曲线
        # 绘图
        plt.plot(r, g_r_avg)
        plt.xlabel("r")
        plt.ylabel("g(r)")
        plt.show()
                
        # ------氢键计算模块------
        #氢键数量，离子键数量，分布


        # ------能量计算模块------
        # 能量分解计算


        # data = {
        #     "cluster_size" : cluster_size,
        #     "cluster_describe" : cluster_describe,
        #     "stats_clusters": stats_clusters,
        #     "stats_participation" : stats_participation,
        #     "Rg_array": Rg_array,
        #     "principal_axes_array" : principal_axes_array,
        #     "anisotropy_array": anisotropy_array,
        #     "axis_ratios_array": axis_ratios_array,
        #     "volume_array" : volume_array,
        #     "density_array": density_array

        # }
        # if self.args.info:
        #     self.plot_data(self.top,data)

        # return data

    @staticmethod
    def plot_data(data):
        
        cluster_size = data["cluster_size"] 
        cluster_describe = data["cluster_describe"] 
        stats_clusters = data["stats_clusters"] 
        stats_participation = data["stats_participation"] 
        Rg_array = data["Rg_array"] 
        principal_axes_array = data["principal_axes_array"] 
        anisotropy_array = data["anisotropy_array"] 
        axis_ratios_array = data["axis_ratios_array"] 
        volume_array = data["volume_array"]
        density_array = data["density_array"] 



        ######cluster_analysis#####
        ploter = MyPlot(nrows=2, ncols=2,figsize=(6,5), dpi= 1080)
        colors = ploter.colors
        #团簇尺寸分布
        ax = 0
        ploter.bar(cluster_size[0], cluster_size[1] ,label= "cluster_size", ax = ax)
        ploter.set_axis(ax = ax, xlabel= "Cluster Size", ylabel= "Probability")

        #团簇信息
        ax = 1
        x2 = cluster_describe['time'].values / 1000
        x2 = x2.reshape(1,-1)
        counts = cluster_describe[['count']].values.T
        ploter.step(x2,counts,ax = ax)
        ploter.set_axis(ax = ax, xlabel= "Time (ns)", ylabel= "Counts")

        y_groups = cluster_describe[['mean','std']].values.T
        labels = ['mean','std']

        ax =  2
        for y, label,color in zip(y_groups, labels, colors):
            ploter.step(x2, y, label=label, color=color,ax=ax) 
        ploter.set_axis(ax = ax, xlabel= "Time (ns)", ylabel= "Cluster Size",legend_loc=1,)

        ax = 3
        x = stats_clusters.index.tolist()
        y1 = stats_clusters['mean'].values
        std1 = stats_clusters['std'].values
        y2 = stats_participation['mean'].values
        std2 = stats_participation['std'].values
        bar_width = 0.4
        x_pos = np.arange(len(x))

        ploter.bar(x_pos - bar_width/2, y1, yerr = std1, 
                   bar_width = bar_width,capsize=4, error_kw={'elinewidth':1, 'ecolor':'black'},
                   color=colors[0],label="Clusters",ax=3)
        ploter.bar(x_pos + bar_width/2, y2, yerr = std2, 
                   bar_width = bar_width,capsize=4, error_kw={'elinewidth':1, 'ecolor':'black'},
                   color=colors[1],label="Participation",ax=3)

        ploter.set_axis(ax = ax, xlabel= "Components", ylabel= "Ratio",)

        # ploter.show(show_if=False, save_path=f"{Path(top_path).stem}_cluster.png")
        ploter.show(show_if=False)


        #计算形状特征            
        ploter = MyPlot(nrows=2, ncols=2, figsize=(6,5), dpi= 1080)
        colors = ploter.colors
        ax = 0
        ploter.probility_densit(Rg_array, ax = ax)
        ploter.set_axis(ax = ax, xlabel= r"Rg ($\AA$)", ylabel= "Probability Density", 
                        # title="Radius of Gyration",
                        grid=True,
                        # xlim=(0,20),
                        )
        
        short_axes = principal_axes_array[:,0]
        middle_axes = principal_axes_array[:,1]
        long_axes = principal_axes_array[:,2]

        ax = 1
        ploter.probility_densit(short_axes, color = colors[0], label='Short axis', ax = ax)
        ploter.probility_densit(middle_axes, color = colors[1], label='Middle axis', ax = ax)
        ploter.probility_densit(long_axes, color = colors[2], label='Long axis', ax = ax)
        ploter.set_axis(ax = ax, xlabel= "Axis", ylabel= "Probability Density",
                        # title=r"Cluster Principal Axes (Rg ($\AA$))",
                        grid=True,
                        legend_loc=1,flegend_size=10,
                        # xlim=(5,15),
                        )
        # ax = 2
        # corr = np.corrcoef(middle_axes, long_axes)[0,1]  # 计算皮尔逊相关系数
        # a, b = np.polyfit(middle_axes, long_axes, 1)     # 拟合直线
        # x_fit = np.array([min(middle_axes), max(middle_axes)])
        # y_fit = a * x_fit + b
        # ploter.line(x_fit, y_fit,color="red", lw=2, label=f'Fit Line (r={corr:.2f})',ax=ax)
        # ploter.scatter(middle_axes, long_axes, alpha=0.5, label='Data Points',ax=ax)
        # ploter.set_axis(ax = ax, xlabel= "Middle Axis Length", ylabel= "Long Axis Length",
        #                 title="Scatter plot with Fit Line",grid=True,)

        # ax = 2
        # ratio_labels = ['L/M', 'M/S', 'L/S']
        # for i, label in enumerate(ratio_labels):
        #     ploter.probility_densit(axis_ratios_array[:,i],label=label, color=colors[i],ax=ax)
        # ploter.set_axis(ax = ax, xlabel= "Axis Ratios", ylabel= "Probability Density",
        #                 title="Cluster Axis Ratios",grid=True,
        #                 # xlim=(0.9,2.0),
        #                 ) 

        ax = 2
        ploter.probility_densit(anisotropy_array, ax = ax)    
        ploter.set_axis(ax = ax, xlabel= r"Shape Anisotropy ($\kappa^2$)", ylabel= "Probability Density",
                        # title="Cluster Shape Anisotropy",
                        grid=True,
                        # xlim=(0,0.6),
                        )   

        ax = 3
        ploter.probility_densit(density_array,ax=ax)
        ploter.set_axis(ax = ax, xlabel= "Cluster density", ylabel= "Probability Density",
                        # title="Cluster density",
                        grid=True,
                        # xlim=(0.25,1.5),
                        ) 
        
        # ax = 4
        # ploter.probility_densit(volume_array,ax=ax)
        # ploter.set_axis(ax = ax, xlabel= "Cluster volume", ylabel= "Probability Density",
        #                 title="Cluster Volume",grid=True,
        #                 # xlim=(500,2500),
        #                 ) 
        # ploter.show(show_if=False, save_path=f"{Path(top_path).stem}_Cinfo.png")
        ploter.show(show_if=False)
