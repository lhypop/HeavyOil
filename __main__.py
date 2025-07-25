# main.py
import sys

from .utils.arg_parser import ArgParser
from .utils.arg_valiate import InputValidator
from .clustring.analyzer import ClusterAnalyzer


def main():
    """
    Main entry point for trajectory analysis application.
    
    Handles:
    - Argument parsing
    - Mode selection
    - Error handling
    """
    try:
        parser = ArgParser()
        args = parser.parse_args()

        valiator = InputValidator()
        valiator.valiate_args(args)
        
        analyzer = ClusterAnalyzer(args)
        analyzer.analyze_aromatic_clusters()
                  
    except Exception as e:
        print(f"Fatal Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()