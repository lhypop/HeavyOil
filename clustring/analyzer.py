import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import List,Tuple
from tqdm import tqdm

from ..edbscan.FastPeriodicDBSCAN import FastPeriodicDBSCAN
from .coordinate import coordinate
from .postprocess import postprocess
from ..analysis.clusteranalysis import ClusterAnalysis
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

            merged_clusters = postprocess(cluster_labels)._postprocss(self.args.infoselect)
            # Store cluster members for current frame
            for cluster in merged_clusters:
                frame_clusters[specific_frame].append(set(cluster))


            cluster_sizes = [int(len(cluster)) for cluster in merged_clusters]
            summarys.append(cluster_sizes)
            cluster_summary.append(merged_clusters)
            cluster_time.append(self.traject.trajectory.time)


        ######cluster_analysis#####
        if self.args.info:
            data_file = Path(f"cluster_nums_{Path(self.top).stem}.txt")
            ClusterInfoAnalyzer = ClusterAnalysis(data_file)

            #限制一下，执行这个命令，需要sellect and all
            pd2 = ClusterInfoAnalyzer.cal_Csize_probability(summarys)
            print(pd2)

            #可选part or all
            pd1 = ClusterInfoAnalyzer.get_cluster_statistics(summarys, cluster_time)
            print(pd1)

            #计算回转半径part
            import matplotlib.pyplot as plt
            import seaborn as sns

            Rg_results = []

            for index, specific_frame in tqdm(enumerate(self.frames), 
                                            total=len(self.frames), 
                                            desc="Calculating the radius of gyration"):
                self.traject.trajectory[specific_frame]  
                cluster_set = cluster_summary[index]

                for set_idx, res_set in enumerate(cluster_set):

                    positions = []
                    masses = []
                    for res in res_set:
                        positions.append(res.atoms.positions)  # shape (N_atoms, 3)
                        masses.append(res.atoms.masses)       # shape (N_atoms,)
                    
                    positions = np.vstack(positions)  # 合并成 (total_atoms, 3)
                    masses = np.concatenate(masses)   # 合并成 (total_atoms,)
                    
                    # 计算质心
                    r_cm = np.sum(positions * masses[:, None], axis=0) / masses.sum()
                    
                    # 计算回转半径
                    Rg = np.sqrt(np.sum(masses[:, None] * (positions - r_cm)**2) / masses.sum())
                    
                    Rg_results.append(Rg)

            Rg_array = np.array(Rg_results)
            # 1. 直方图归一化为概率密度
            plt.hist(Rg_array, bins=10, density=True, alpha=0.6, color='skyblue', label='Histogram')

            # 2. KDE 平滑曲线
            sns.kdeplot(Rg_array, bw_adjust=0.5, color='red', label='KDE')

            plt.xlabel("Radius of Gyration (Rg)")
            plt.ylabel("Probability Density")
            plt.title("Rg Probability Density Distribution")
            plt.legend()
            plt.savefig(f"{Path(self.top).stem}.png")
            # plt.show()

            



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
