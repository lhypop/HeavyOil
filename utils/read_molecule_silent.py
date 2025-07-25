from contextlib import redirect_stderr, nullcontext
from openbabel import pybel
from pathlib import Path
from typing import Optional, Set, Union
import os
import tempfile
import shutil
import logging

class MoleculeReader:
    """Handles molecular file reading and conversion with error suppression."""
    
    # Valid chemical elements (up to atomic number 20)
    VALID_ELEMENTS = {
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca'
    }

    def __init__(self, suppress_warnings: bool = True):
        """
        Args:
            suppress_warnings: Whether to suppress OpenBabel warnings
        """
        self.suppress_warnings = suppress_warnings

    def fix_pdb_element_column(self, pdb_in: Path, backup: bool = True) -> None:
        """Corrects element symbols in PDB files (columns 77-78).
        
        Args:
            pdb_in: Path to input PDB file
            backup: Whether to create a backup file
            
        Raises:
            ValueError: For invalid atom names or elements
            OSError: For file operation failures
        """
        # Validate input file
        if not pdb_in.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_in}")

        # Create backup if requested
        if backup:
            backup_path = pdb_in.with_suffix('.pdb.bak')
            try:
                shutil.copy2(pdb_in, backup_path)
            except OSError as e:
                raise OSError(f"Backup failed: {e}") from e

        # Process in temporary file
        temp_dir = pdb_in.parent
        try:
            with tempfile.NamedTemporaryFile(
                mode='w',
                dir=str(temp_dir),
                prefix=f"{pdb_in.stem}_",
                suffix='.tmp',
                delete=False
            ) as tmp_file:
                temp_path = Path(tmp_file.name)
                
                try:
                    # Process each line
                    with open(pdb_in) as src_file:
                        for line in src_file:
                            if line.startswith(("ATOM", "HETATM")):
                                line = self._process_pdb_line(line)
                            tmp_file.write(line)
                    
                    # Ensure write completion
                    tmp_file.flush()
                    os.fsync(tmp_file.fileno())
                    
                    # Atomic replacement
                    shutil.move(str(temp_path), str(pdb_in))
                    
                except Exception:
                    temp_path.unlink(missing_ok=True)
                    if backup:
                        try:
                            shutil.move(str(backup_path), str(pdb_in))
                        except OSError:
                            pass
                    raise
                    
        except OSError as e:
            raise OSError(f"Temp file creation failed: {e}") from e

    def _process_pdb_line(self, line: str) -> str:
        """Process a single PDB line to correct element column."""
        atom_name = line[12:16].strip()
        element = ''.join(c for c in atom_name if c.isalpha())
        
        if not element:
            raise ValueError(f"Invalid atom name: {atom_name}")
        
        # Normalize element symbol
        element = element.capitalize()
        if len(element) > 1 and element[:2] in self.VALID_ELEMENTS:
            element = element[:2]
        else:
            element = element[0]

        if element not in self.VALID_ELEMENTS:
            raise ValueError(f"Invalid element {element} from {atom_name}")

        return line[:76] + f"{element:>2}" + line[78:]

    def read_molecule(
        self,
        file_path: Union[str, Path],
        output_format: str = "mol2",
        output_path: Optional[Union[str, Path]] = None
    ) -> Optional[pybel.Molecule]:
        """Read a molecule file with optional format conversion.
        
        Args:
            file_path: Input file path
            output_format: Desired output format
            output_path: Optional output path
            
        Returns:
            pybel.Molecule or None if reading fails
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

        # Auto-fix PDB files
        if file_path.suffix.lower() == '.pdb':
            self.fix_pdb_element_column(file_path)

        # Prepare output
        output_path = Path(output_path) if output_path else None
        if output_path and output_path.exists():
            logging.warning(f"Overwriting existing file: {output_path}")

        # Context manager for warning suppression
        err_context = nullcontext() if not self.suppress_warnings else redirect_stderr(open(os.devnull, 'w'))

        try:
            with err_context:
                # Read molecule
                mol = next(pybel.readfile(file_path.suffix[1:], str(file_path)), None)
                
                # Write output if requested
                if mol and output_path:
                    mol.write(output_format, str(output_path), overwrite=True)
                
                return mol
                
        except StopIteration:
            return None
        except Exception as e:
            logging.error(f"Error reading {file_path}: {str(e)}")
            raise