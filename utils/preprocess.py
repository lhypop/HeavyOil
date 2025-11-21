import yaml
import numpy as np
from pathlib import Path
import MDAnalysis as mda

class Preprocessor:
    def __init__(self, args):
        """Initialize with preprocessing results directly available"""
        self.args = args
        self._setup()
    
    @property
    def parameters(self):
        """Return all analysis parameters in a structured way"""
        return {
            'residues': self.residues,
            'frames': self.frames,
            'eps': self.eps,
            'min_samples': self.min,
            'universe': self.traject,  # MDAnalysis universe if needed
            'components':self.components
        }
    
    def _setup(self):
        """Configure preprocessing (internal use only)"""
        # Initialize trajectory
        self.traject = mda.Universe(self.args.gro, self.args.trj or self.args.gro)
        
        # Process selections
        self.residues = self._get_residues()
        self.frames = self._get_frames()
        self.eps, self.min = self._get_config()
        self.components = self._get_components()
        
        print(f"Preprocessing complete - {len(self.residues)} residues, "
              f"{len(self.frames)} frames")
         

    def _get_residues(self):
        """Get selected residues"""
        if self.args.residueselect == 'all':
            self.traject.trajectory[0]
            return list(set(self.traject.residues.resnames))
            
        names = [x.strip() for x in self.args.residueselect.split(",") if x.strip()]
        resnames = [self.traject.select_atoms(f"resname {n}").residues.resnames for n in names]
        return list(set(np.concatenate(resnames).tolist()))

    def _get_frames(self):
        """Validate and get frame range"""
        total = len(self.traject.trajectory)
        begin = self.args.begin or 0
        end = self.args.end or total
        interval = self.args.interval or 1
        
        if not (0 <= begin < end <= total):
            raise ValueError(f"Invalid frames: {begin}-{end} (max: {total})")
            
        return range(begin, end, interval)
    
    def _get_components(self):
        components = {
            'sat':[],
            'aro':[],
            'res':[],
            'asp':[]
        }

        sat = self.args.sat
        aro = self.args.aro
        res = self.args.res
        asp = self.args.asp

        if all([sat,aro,res,asp]):
            sat_list = [x.strip() for x in sat.split(",") if x.strip()]
            aro_list = [x.strip() for x in aro.split(",") if x.strip()]
            res_list = [x.strip() for x in res.split(",") if x.strip()]
            asp_list = [x.strip() for x in asp.split(",") if x.strip()]
            
            all_lists = [sat_list, aro_list, res_list, asp_list]
            keys = list(components.keys())

            for i, coms in enumerate(all_lists):
                for res_name in coms:
                    sel_resnames = self.traject.select_atoms(f"resname {res_name}").residues.resnames
                    if len(sel_resnames) > 0:
                        components[keys[i]].extend(sel_resnames)
                    else:
                        print(f"Warning: {res_name} specified for {keys[i]} but not found in trajectory.")
                        
        elif any([sat,aro,res,asp]):
            raise ValueError("You must specify all four component groups: --sat, --aro, --res, and --asp.")
                    
        for k in components:
            components[k] = list(set(components[k]))

        return components



    def _get_config(self):
        """Load config parameters"""
        config_path = Path(__file__).parent.parent / "config/config.yaml"
        config = self._load_config(config_path)  # Implement YAML loading
        return config['eps'], config['min_samples']
    
    def _load_config(self,config_path) -> dict:
        """
        Load configuration parameters from the YAML config file.
        """
        if not config_path.exists():
            print(f"Looking for config file at: {config_path.resolve()}")
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config