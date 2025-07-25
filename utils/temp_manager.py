"""
temp.py

Module for managing temporary directories, including specifying paths, 
creating folders, and clearing their contents.

Author: Haiyang Li
"""

from pathlib import Path
import shutil

class TempMangager:
    """
    Manage temporary directory: specify, create, clear contents.
    """

    def __init__(self, temp_dir: str):
        self.temp_dir = Path(temp_dir)

    def ensure_dir(self):
        """Ensure the temp directory exists."""
        if not self.temp_dir.exists():
            self.temp_dir.mkdir(parents=True, exist_ok=True)
    
    def clear(self):
        """Clear all files and subdirectories inside the temp directory."""
        if not self.temp_dir.exists():
            return
        
        for item in self.temp_dir.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)

    def prepare(self):
        """
        Ensure the temp directory exists and is empty.
        """

        self.ensure_dir()
        self.clear()
