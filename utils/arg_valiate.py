import argparse
from pathlib import Path
from typing import Optional

class InputValidator:
    """Validates molecular dynamics simulation input parameters.
    
    This class provides comprehensive validation for:
    - File existence and formats
    - Parameter value ranges
    - Logical consistency between parameters
    """

    def __init__(self):
        """Initialize validator with input arguments.
        """
        self.args = None
        self.top = False
        self.trj = False

    def valiate_args(self, args: argparse.Namespace) -> None:
        """Execute complete validation pipeline.
        
        Raises:
            FileNotFoundError: For missing files
            ValueError: For invalid values/combinations
            TypeError: For incorrect data types
        """
        self.args = args
        self._validate_topology()
        self._validate_trajectory()
        self._validate_frames()
        self._validate_dependencies()

    def _validate_topology(self) -> None:
        """Validate topology file parameters."""
        gro_path = Path(self.args.gro)
        if not gro_path.is_file():
            raise FileNotFoundError(f"Topology file not found: {gro_path}")
        if gro_path.suffix != '.gro':
            raise ValueError(f"Invalid topology format. Expected .gro, got {gro_path.suffix}")
        self.top = True

    def _validate_trajectory(self) -> None:
        """Validate trajectory file parameters if provided."""
        if self.args.trj:
            trj_path = Path(self.args.trj)
            if not trj_path.is_file():
                raise FileNotFoundError(f"Trajectory file not found: {trj_path}")
            if trj_path.suffix != '.xtc':
                raise ValueError(f"Invalid trajectory format. Expected .xtc, got {trj_path.suffix}")
        else:
            self.trj = False


    def _validate_frames(self) -> None:
        """Validate frame selection parameters."""
        self._validate_individual_frame('begin', minimum=0)
        self._validate_individual_frame('end')
        self._validate_interval('interval')

    def _validate_individual_frame(self, param: str, minimum: Optional[int] = None) -> None:
        """Validate a single frame parameter.
        
        Args:
            param: Parameter name ('begin' or 'end')
            minimum: Minimum allowed value
        """
        if not hasattr(self.args, param) or getattr(self.args, param) is None:
            return

        value = getattr(self.args, param)
        if not isinstance(value, int):
            raise TypeError(f"{param} frame must be integer. Got {type(value)}")
        if minimum is not None and value < minimum:
            raise ValueError(f"{param} frame cannot be below {minimum}. Got {value}")

    def _validate_interval(self,param) -> None:
        """Validate frame interval parameter."""
        if not hasattr(self.args, param) or getattr(self.args, param) is None:
            return

        if not isinstance(self.args.interval, int):
            raise TypeError(f"Interval must be integer. Got {type(self.args.interval)}")
        if self.args.interval <= 0:
            raise ValueError(f"Interval must be positive. Got {self.args.interval}")


    def _validate_dependencies(self) -> None:
        """Validate relationships between parameters."""
        has_frames = (hasattr(self.args, 'begin') and self.args.begin is not None) or \
                    (hasattr(self.args, 'end') and self.args.end is not None)
        
        if has_frames and (not hasattr(self.args, 'trj') or not self.args.trj):
            raise ValueError("Frame selection requires trajectory file")

        if (hasattr(self.args, 'begin') and self.args.begin is not None and 
            hasattr(self.args, 'end') and self.args.end is not None and 
            self.args.end <= self.args.begin):
            raise ValueError(f"End frame ({self.args.end}) must exceed begin frame ({self.args.begin})")