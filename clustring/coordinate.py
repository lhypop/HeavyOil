import numpy as np
import MDAnalysis as mda
from rdkit import Chem
from pathlib import Path
from typing import List
from collections import Counter

from ..utils.read_molecule_silent import MoleculeReader

class coordinate:
    def __init__(self, traj, temp):
        self.traj = traj
        self.temp = temp
        self.box = self.traj.dimensions[:3] / 10 

    def _atoms_in_aromatics(self,mol2_file: Path) -> List[List[int]]:
            """
            Identify atom indices of aromatic rings in a molecule file.

            Parameters
            ----------
            mol2_file : Path
                Path to a .mol2 file containing the molecule.
            Returns
            -------
            List[List[int]]
                A list of rings, each ring is a list of atom indices.
            """
            mol = Chem.MolFromMol2File(str(mol2_file), sanitize=False)
            if mol is None:
                raise ValueError(f"Failed to load molecule from {mol2_file}")

            ring_info = mol.GetRingInfo()
            ring_atoms = ring_info.AtomRings()

            aromatic_rings = []

            for ring in ring_atoms:
                atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                symbols = [atom.GetSymbol() for atom in atoms]
                degrees = [atom.GetTotalDegree() for atom in atoms]

                if len(ring) == 6:
                    # For 6-membered rings like benzene/pyridine
                    expected_degrees = [3 if s == 'C' else 2 for s in symbols]
                elif len(ring) == 5:
                    # For 5-membered rings like pyrrole, imidazole
                    expected_degrees = [3 if s in ('C', 'N') else 2 for s in symbols]
                else:
                    continue  # Skip non-aromatic-size rings

                if degrees == expected_degrees:
                    aromatic_rings.append(ring)

            return aromatic_rings

    def _get_center_of_aroring(self):
        """
        Extract geometric centers of aromatic rings in a specified frame.

        Parameters
        ----------
        frame : int
            Frame index to analyze. Use -1 for the last frame.
        """
        residues = self.traj.residues

        ring_centers = []
        aromatic_residues = []

        for _, res in enumerate(residues):
            temp_pdb = self.temp / f"{res.resname}.pdb"
            mol2_file = temp_pdb.with_suffix(".mol2")

            # If mol2 already exists, skip conversion
            if not mol2_file.exists():
                res.atoms.write(temp_pdb)
                reader = MoleculeReader()
                _ = reader.read_molecule(temp_pdb, output_path=mol2_file)

            aromatic_rings = self._atoms_in_aromatics(mol2_file)

            if len(aromatic_rings) != 0:
                atoms_coord = res.atoms.positions

                for ring in aromatic_rings:
                    #distance unit change from nm to angstrom
                    ring_atom_coords = np.array([atoms_coord[atom_idx] for atom_idx in ring])/10
                    ring_center = ring_atom_coords.mean(axis=0)

                    ring_centers.append(ring_center)
                    aromatic_residues.append(res)

        gro_ring_centers = np.array(np.vstack(ring_centers))

        return gro_ring_centers,aromatic_residues,self.box