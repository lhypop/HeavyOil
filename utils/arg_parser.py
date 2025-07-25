# arg_parser.py
import argparse

class ArgParser:
    """Command line argument parser for trajectory analysis tool."""
    
    def __init__(self):
        self.parser = self._create_parser()
    
    def _create_parser(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Molecular Dynamics Trajectory Analysis Tool",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Examples:\n"
                   "  python main.py -s system.gro -t traj.xtc\n"
                   "  python main.py -s system.gro --info"
        )
        
        self._add_input_arguments(parser)
        self._add_frame_arguments(parser)
        self._add_clusterinfo_arguments(parser)
        self._add_analysis_arguments(parser)
        return parser
    
    def _add_input_arguments(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_argument_group('Input files')
        group.add_argument(
            "-s", "--gro", 
            required=True,
            help="Topology file (GRO format)"
        )
        group.add_argument(
            "-t", "--trj",
            help="Trajectory file (XTC format)"
        )
    
    def _add_frame_arguments(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_argument_group('Frame selection')
        group.add_argument(
            "-b", "--begin",
            type=int,
            default=0,
            help="Start frame index (0-based)"
        )
        group.add_argument(
            "-e", "--end",
            type=int,
            help="End frame index (exclusive)"
        )
        group.add_argument(
            "-i", "--interval",
            type=int,
            default=1,
            help="Frame sampling interval"
        )
    
    def _add_clusterinfo_arguments(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_argument_group('Cluster information options')
        group.add_argument(
            "--info",
            action="store_true",
            help="Perform info analysis"
        )
        group.add_argument(
        "--infoselect",
        choices=['part', 'all'],
        default='part',
        help="Analysis scope: 'part'=clustered, 'all'=entire system"
        )

    def _add_analysis_arguments(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_argument_group('Analysis options')
        group.add_argument(
            "--residueselect",
            default='all',
            help="Residue selection ('all' or specify names (e.g. 'ASP1,ASN2' or 'AS*'))"
        )
    
    def parse_args(self) -> argparse.Namespace:
        """Parse and validate arguments"""
        args = self.parser.parse_args()
        self._validate_args(args)
        return args
    
    def _validate_args(self, args: argparse.Namespace) -> None:
        if args.end is not None and args.begin >= args.end:
            self.parser.error("End frame must be > start frame")
        if args.interval <= 0:
            self.parser.error("Interval must be > 0")