import os
import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict
from collections import Counter
from typing import List,Tuple
from tqdm import tqdm
import MDAnalysis as mda

from ..edbscan.FastPeriodicDBSCAN import FastPeriodicDBSCAN
from .coordinate import coordinate
from .postprocess import postprocess
from ..analysis.clusteranalysis import ClusterAnalysis
from ..analysis.energyanalysis import EnergyAnalysis
from ..analysis.interactanalysis import Interactanalysis

from ..visual.generate import GMXClusterTool
from ..visual.translate import clusterTranslator
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
        self.xtc = args.trj
        
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
        self.components = params['components']
        
        # Set up temporary directory for intermediate files
        # Default location: ./temp relative to current working directory
        self.temp = Path.cwd() / "temp"
        self.temp.mkdir(exist_ok=True)  # Create if doesn't exist
    
    
    def _DBSCAN_clustering(self,
                coords: np.ndarray,
                aromatic_residues: List[str],
                box) -> List[Tuple[int, str]]:
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
                                   box_size=box)
        cluster_labels = model.fit_predict(coords)

        clustered_residues = [
            (label, res)
            for label, res in zip(cluster_labels, aromatic_residues)
            if label != -1  # -1 means noise
        ]

        return clustered_residues
    

    def analyzing(self, specific_frame, topfile, xtcfile, temp_file, residues):
        # traject = mda.Universe(topfile, xtcfile)
        traject = mda.Universe(topfile, xtcfile)
        traject.trajectory[specific_frame] 

        temp = temp_file / str(specific_frame)
        temp.mkdir(exist_ok=True)
        coordinator = coordinate(temp, residues)

        # Get aromatic ring centers and corresponding residues
        gro_ring_centers, aromatic_residue, box = coordinator._get_center_of_aroring(traject)

        cluster_labels = self._DBSCAN_clustering(gro_ring_centers, aromatic_residue, box)
        merged_clusters = postprocess(cluster_labels)._postprocess(self.args.infoselect)

        return specific_frame, merged_clusters
    
    
    # def analyze_aromatic_clusters(self):

    #     frame_clusters = defaultdict(list)

    #     from concurrent.futures import ProcessPoolExecutor

    #     from functools import partial
    #     func = partial(
    #         self.analyzing,
    #         topfile=self.args.gro,
    #         xtcfile=self.args.trj,
    #         temp_file=self.temp,
    #         residues=self.residues
    #     )

    #     # 多进程执行
    #     with ProcessPoolExecutor(max_workers=1) as executor:
    #         results = list(tqdm(executor.map(func, self.frames), total=len(self.frames), desc="cluster processing"))

    #     # 收集结果
    #     for frame, merged_clusters in results:
    #         for cluster in merged_clusters:
    #             frame_clusters[frame].append(set(cluster))


    #     ######visualization########
    #     all_cluster_list = frame_clusters.get(self.frames[-1])

    #     # Generate output files
    #     output_base = Path.cwd() / Path(self.top).stem
    #     move_tcl_path = output_base.with_name(f"{output_base.name}_move.tcl")
    #     clusterTranslator(all_cluster_list, self.traject.trajectory[self.frames[-1]].dimensions[:3], move_tcl_path)
    
        
    #     processor = GMXClusterTool(
    #         molname = str(output_base)
    #     )
    #     processor.process_clusters(all_cluster_list)

    #     if self.args.info:

    #         features = self.extract_characteristic(frame_clusters)

    #         self.features_set(features, filename=f"{output_base.name}_features.json")

    #         return features

    def analyze_aromatic_clusters(self):

        frame_clusters = defaultdict(list)

        # 顺序执行（去掉多进程）
        for specific_frame in tqdm(self.frames, desc="cluster processing"):
            frame, merged_clusters = self.analyzing(
                specific_frame,
                topfile=self.args.gro,
                xtcfile=self.args.trj,
                temp_file=self.temp,
                residues=self.residues
            )
            for cluster in merged_clusters:
                frame_clusters[frame].append(set(cluster))

        ###### visualization ########
        all_cluster_list = frame_clusters.get(self.frames[-1])

        # Generate output files
        output_base = Path.cwd() / Path(self.top).stem
        
        move_tcl_path = output_base.with_name(f"{output_base.name}_move.tcl")
        clusterTranslator(
            all_cluster_list,
            self.traject.trajectory[self.frames[-1]].dimensions[:3],
            move_tcl_path
        )

        processor = GMXClusterTool(
            molname=str(output_base)
        )
        processor.process_clusters(all_cluster_list)


        if self.args.info:
            features = self.extract_characteristic(frame_clusters)
            return features

    def extract_characteristic(self, frame_clusters):
        """
        从聚集体轨迹中提取多维特征信息，包括：
        - 数量特征
        - 时间特征
        - 组分特征
        - 形状特征
        - 分布特征
        - 杂原子与氢键信息
        - 能量特征
        """
        import json

        # ====== 初始化 ======
        data_file = Path(f"cluster_nums_{Path(self.top).stem}.txt")
        ClusterInfoAnalyzer = ClusterAnalysis(self.traject, data_file, frame_clusters)

        Features_init = {}
        Features_set = {}

        # ====== 辅助函数 ======
        def get_stats(cen, p):
            """根据分布中心 cen 和概率密度 p 计算加权平均与方差"""
            bin_width = cen[1] - cen[0] if len(cen) > 1 else 1.0
            mean_val = np.sum(p * cen * bin_width)
            mean_square = np.sum(p * (cen**2) * bin_width)
            var_cal = mean_square - mean_val**2
            return mean_val, var_cal

        # ====== 数量特征模块 ======
        if "All" in self.args.info or "SizeInfo" in self.args.info:
            size_dis = ClusterInfoAnalyzer.cal_Csize_probability()
            Features_init["size_dis"] = size_dis

            size, p = size_dis["size"], size_dis["p"]
            avg_size = np.sum(size * p)
            var_size = np.sum(p * (size - avg_size) ** 2)
            Features_set.update({
                "avg_size": avg_size,
                "var_size": var_size,
                "max_size": size[-1],
            })

        # ====== 时间特征模块 ======
        if "All" in self.args.info or "TimeInfo" in self.args.info:
            time_stats = ClusterInfoAnalyzer.cal_Cnums_overtime()
            Features_init["time_stats"] = time_stats

            counts = time_stats['count'].values
            sizes = time_stats['mean'].values
            Features_set.update({
                "Tavg_num": np.mean(counts),
                "Tvar_num": np.std(counts),
                "Tavg_size": np.mean(sizes),
                "Tvar_size": np.std(sizes),
            })

        # ====== 组分特征模块 ======
        if "All" in self.args.info or "Components" in self.args.info:
            if not self.components:
                raise ValueError("You need to specify the four components (sat, aro, res, asp).")

            res_counts = Counter(self.traject.residues.resnames)
            df_clusters, df_participation = ClusterInfoAnalyzer.cal_components_ratio(self.components, res_counts)
            Features_init["df_clusters"] = df_clusters
            Features_init["df_participation"] = df_participation

            clu_com = df_clusters.mean().tolist()
            par_com = df_participation.mean().tolist()

            Features_set.update({
                "Clust_sat": clu_com[0], "Parti_sat": par_com[0],
                "Clust_aro": clu_com[1], "Parti_aro": par_com[1],
                "Clust_res": clu_com[2], "Parti_res": par_com[2],
                "Clust_asp": clu_com[3], "Parti_asp": par_com[3],
            })

        # ====== 形状特征模块 ======
        if "All" in self.args.info or "Shape" in self.args.info:
            shape_features = ClusterInfoAnalyzer.cal_shape_describe()
            Features_init["shape_features"] = shape_features

            sf = shape_features  # 简写
            feature_keys = ["Rg", "PAS", "PAM", "PAL", "AS", "ALM", "AMS", "ALS", "VO", "DE"]

            for key in feature_keys:
                mean_val, var_val = get_stats(sf[key]["cen"], sf[key]["P"])
                Features_set[f"{key}_mean"] = mean_val
                Features_set[f"{key}_var"] = var_val

        # ====== 分布特征模块 ======
        if "All" in self.args.info or "Distribution" in self.args.info:
            rdf_stats, nnd_stats = ClusterInfoAnalyzer.cal_distribution_analysis()
            Features_init["rdf_stats"] = rdf_stats
            Features_init["nnd_stats"] = nnd_stats

            rdf_means = []
            for i in range(5):
                mask = (rdf_stats["r"] > i) & (rdf_stats["r"] < i + 1)
                rdf_means.append(np.mean(rdf_stats["gr"][mask]))

            nnd_mean, nnd_var = get_stats(nnd_stats["cen"], nnd_stats["p"])

            Features_set.update({
                "rdf_01": rdf_means[0], "rdf_12": rdf_means[1],
                "rdf_23": rdf_means[2], "rdf_34": rdf_means[3],
                "rdf_45": rdf_means[4],
                "nnd_mean": nnd_mean, "nnd_var": nnd_var,
            })

        # ====== 杂原子与氢键特征 ======
        if "All" in self.args.info or "Heteratom" in self.args.info:
            InteractAnalyst = Interactanalysis(self.traject, frame_clusters)
            Heter_info = InteractAnalyst.Heter_statistic_analysis()

            Hydro_select = " or ".join([f"resname {res}" for res in self.residues])
            Hydro_info = InteractAnalyst.Hydro_statistic_analysis(select_info=Hydro_select)

            Features_init["Heter_info"] = Heter_info
            Features_init["Hydro_info"] = Hydro_info

            Heter_mean, Heter_var = Heter_info.mean(), Heter_info.var()
            Hydro_mean, Hydro_var = Hydro_info.mean(), Hydro_info.var()

            Features_set.update({
                "I_heter_mean": Heter_mean[1], "I_heter_var": Heter_var[1],
                "O_heter_mean": Heter_mean[2], "O_heter_var": Heter_var[2],
                "S_heter_mean": Heter_mean[3], "S_heter_var": Heter_var[3],

                "I_hydro_mean": Hydro_mean[1], "I_hydro_var": Hydro_var[1],
                "O_hydro_mean": Hydro_mean[2], "O_hydro_var": Hydro_var[2],
                "S_hydro_mean": Hydro_mean[3], "S_hydro_var": Hydro_var[3],
            })

        # ====== 能量特征模块 ======
        if "All" in self.args.info or "Energy" in self.args.info:
            EAnalysis = EnergyAnalysis(
                self.top,
                self.xtc,
                self.traject.trajectory.dt,
                frame_clusters,
                Path.cwd(),
            )
            Energy_info = EAnalysis.cal_energy()
            Features_init["Energy_info"] = Energy_info

            E_desc = Energy_info.describe()
            E_mean = E_desc.loc["mean"]

            Features_set.update({
                "All_Energy": E_mean[0], "All_VDW": E_mean[1], "All_Coulom": E_mean[2],
                "Clu_Energy": E_mean[3], "Clu_VDW": E_mean[4], "Clu_Coulom": E_mean[5],
                "Cohe_Energy": E_mean[6], "Cohe_VDW": E_mean[7], "Cohe_Coulom": E_mean[8],
                "inclu_Energy": E_mean[9], "inclu_VDW": E_mean[10], "inclu_Coulom": E_mean[11],
            })

        # ====== 保存特征 ======
        output_file = f"{Path(self.top).stem}_features.json"
        with open(output_file, "w") as f:
            json.dump(Features_set, f, indent=4)

        print(f"[Feature Extraction Completed] Results saved to {output_file}")
        return Features_init
