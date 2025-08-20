from contextlib import redirect_stderr, nullcontext
from openbabel import pybel
from pathlib import Path
import os
import tempfile
import shutil

class MoleculeReader:
    """Simplified molecule reader with PDB element correction and format conversion."""
    
    VALID_ELEMENTS = {'H', 'C', 'N', 'O', 'S'}

    def __init__(self, suppress_warnings: bool = True):
        self.suppress_warnings = suppress_warnings

    def fix_pdb_element_column(self, pdb_path: Path) -> None:
        """Fix element symbols in PDB files (columns 77-78)."""
        temp_dir = pdb_path.parent
        with tempfile.NamedTemporaryFile(mode='w', dir=str(temp_dir), delete=False) as tmp_file:
            temp_path = Path(tmp_file.name)
            for line in open(pdb_path):
                if line.startswith(("ATOM", "HETATM")):
                    atom_name = line[12:16].strip()
                    element = ''.join(c for c in atom_name if c.isalpha()).capitalize()
                    if len(element) > 1 and element[:2] in self.VALID_ELEMENTS:
                        element = element[:2]
                    else:
                        element = element[0]
                    line = line[:76] + f"{element:>2}" + line[78:]
                tmp_file.write(line)
        shutil.move(str(temp_path), str(pdb_path))

    def read_molecule(self, file_path: Path, output_path: Path = None):
        """Read molecule file and optionally write to another format."""
        file_path = Path(file_path)
        if file_path.suffix.lower() == '.pdb':
            self.fix_pdb_element_column(file_path)

        err_context = nullcontext() if not self.suppress_warnings else redirect_stderr(open(os.devnull, 'w'))
        with err_context:
            mol = next(pybel.readfile("pdb", str(file_path)), None)
            if mol and output_path:
                mol.write("mol2", str(output_path), overwrite=True)
            return mol
