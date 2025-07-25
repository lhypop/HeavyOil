import os
import re
import subprocess
from collections import Counter
from typing import List, Set

class GMXClusterTool:
    """Handles processing and visualization of aromatic residue clusters in GROMACS.
    
    Features:
    - Generates NDX index files for cluster selections
    - Creates VMD visualization scripts
    - Extracts cluster coordinates
    - Prepares topology files for extracted clusters
    """

    def __init__(self, molname: str):
        """Initialize processor with base molecule name.
        
        Args:
            molname: Base name for input/output files (without extensions)
        """
        self.molname = molname
        self._validate_input_files()

    def _validate_input_files(self) -> None:
        """Check required input files exist."""
        required_files = [f"{self.molname}.gro"]
        for f in required_files:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Required input file missing: {f}")

    def process_clusters(self, filtered_clusters: List[Set]) -> None:
        """Main pipeline for processing residue clusters.
        
        Args:
            filtered_clusters: List of residue sets representing clusters
        """
        sorted_residues = self._prepare_residue_list(filtered_clusters)
        
        # Generate output files
        self._create_index_file(sorted_residues)
        self._generate_vmd_scripts(filtered_clusters)
        self._extract_cluster_coordinates()
        self._generate_topology_file(sorted_residues)

    def _prepare_residue_list(self, clusters: List[Set]) -> List:
        """Prepare sorted unique residue list from clusters."""
        residues = []
        for cluster in clusters:
            residues.extend(cluster)
        return sorted(list(set(residues)), key=lambda res: res.resid)

    def _create_index_file(self, residues: List) -> None:
        """Generate GROMACS index file for clusters."""
        self._cleanup_file(f"{self.molname}.ndx")
        
        n_groups = self._get_existing_group_count()
        selection = (
            f'ri {" ".join(str(res.resid) for res in residues)}\n'
            f'name {n_groups} aro_zone\n'
            'q\n'
        )
        
        self._run_gmx_command(
            ['gmx', 'make_ndx', '-f', f"{self.molname}.gro", '-o', f"{self.molname}.ndx"],
            input=selection
        )

    def _get_existing_group_count(self) -> int:
        """Count existing groups in GRO file."""
        result = self._run_gmx_command(
            ['gmx', 'make_ndx', '-f', f"{self.molname}.gro", '-o', f"{self.molname}.ndx"],
            input='l\nq\n'
        )
        
        return len([
            line for line in result.stdout.splitlines() 
            if re.match(r'\s*\d+\s+\S.*?:[^:]*$', line)
        ])

    def _generate_vmd_scripts(self, clusters: List[Set]) -> None:
        """Generate VMD visualization scripts."""
        # Main visualization script
        self._cleanup_file(f"{self.molname}.tcl")
        with open(f"{self.molname}.tcl", "w") as f:
            for rep_id, cluster in enumerate(clusters, 1):
                color_id = (rep_id - 1) % 33
                f.write(
                    f"mol selection resid {' '.join(str(res.resid) for res in cluster)}\n"
                    f"mol addrep 0\n"
                    f"mol modstyle {rep_id} 0 VDW 1.000000 12.000000\n"
                    f"mol modcolor {rep_id} 0 ColorID {color_id}\n"
                    f"mol modmaterial {rep_id} 0 Opaque\n"
                )

        # Line representation script
        self._cleanup_file(f"{self.molname}_line.tcl")
        with open(f"{self.molname}_line.tcl", "w") as f:
            residues = [res.resid for cluster in clusters for res in cluster]
            f.write(
                f"mol selection resid {' '.join(map(str, residues))}\n"
                "mol addrep 0\n"
            )

    def _extract_cluster_coordinates(self) -> None:
        """Extract cluster coordinates to new GRO file."""
        self._cleanup_file(f"{self.molname}_aro_extract.gro")
        self._run_gmx_command(
            ['gmx', 'trjconv', '-f', f"{self.molname}.gro", 
             '-s', f"{self.molname}.gro", '-n', f"{self.molname}.ndx",
             '-o', f"{self.molname}_aro_extract.gro"],
            input="aro_zone\n"
        )

    def _generate_topology_file(self, residues: List) -> None:
        """Generate topology file for extracted clusters."""
        self._cleanup_file(f"{self.molname}_aro_extract_reindex.top")
        res_counts = Counter(res.resname for res in residues)
        
        with open(f"{self.molname}_aro_extract_reindex.top", "w") as f:
            f.write('#include "../topol/charmm36-jul2021.ff/forcefield.itp"\n\n')
            for resname in res_counts:
                f.write(f'#include "../topol/itp/{resname}.itp"\n')
            f.write('\n[system]\nSARA\n\n[molecules]\n; Compound   nmols\n')
            for resname, count in res_counts.items():
                f.write(f'{resname}\t{count}\n')

    def _run_gmx_command(self, command: List[str], input: str = "") -> subprocess.CompletedProcess:
        """Execute GROMACS command with error handling."""
        try:
            return subprocess.run(
                command,
                input=input,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"GROMACS command failed: {e.stderr}")
            raise

    def _cleanup_file(self, filename: str) -> None:
        """Remove file if it exists."""
        if os.path.exists(filename):
            os.remove(filename)