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
        cluster_time = []

        for specific_frame in tqdm(self.frames, desc="Processing"):
            # Get aromatic ring centers and corresponding residues
            self.traject.trajectory[specific_frame] 
            gro_ring_centers, aromatic_residue, self.box = coordinate(self.traject, self.temp)._get_center_of_aroring()
            cluster_labels = self._DBSCAN_clustering(gro_ring_centers, aromatic_residue)

            tolnum, merged_clusters = postprocess(cluster_labels)._postprocss(self.args.infoselect)

            # Store cluster members for current frame
            for cluster in merged_clusters:
                frame_clusters[specific_frame].append(set(cluster))

            if self.args.info:
                info_cluster = ClusterAnalysis().get_cluster_statistics(merged_clusters,tolnum)
                summarys.append(info_cluster)
                cluster_time.append(self.traject.trajectory.time)

        if self.args.info:
            data_file = Path.cwd() / f"data/cluster_nums_{Path(self.top).stem}.txt"
            ClusterAnalysis().save_outcome(summarys, cluster_time, data_file)
            

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
