# HeavyOil
Molecular dynamics simulation analysis scripts of heavy oil systems.

---

## Version History

### origin
- Basic functionalities for recognizing structures.
- Initial output of structural information.

### 0.0.1
- Improved frame selection and residue filtering.  
- Added analysis of cluster number and shape characteristics.
- Enabled visualization of analysis results.

### 0.0.2 (Maintenance/0.0.2, under development)
- Planned: Advanced cluster information analysis.  


---

## Overview
**HeavyOil** is a command-line tool for analyzing molecular dynamics (MD) trajectories, with a focus on cluster analysis and structural information extraction.  
It is designed to work with **GROMACS** output files (`.gro` and `.xtc`) and provides flexible frame selection and residue-based analysis.

---

## Features
- 📂 Load and analyze MD trajectories (`.gro` + `.xtc`).
- 🎯 Flexible frame selection (begin, end, interval).
- 🔍 Cluster-based information analysis (`--info`).
- 🌐 Support for different analysis scopes:  
  - **part**: clustered molecules only  
  - **all**: the entire system  
- 🧬 Residue-based filtering with wildcard support (e.g., `AS*`).
- 🖥️ Simple command-line interface, easy integration into workflows.

---

## Installation
Clone this repository and install dependencies:

```bash
git clone https://github.com/yourusername/HeavyOil.git

cd HeavyOil

pip install -r requirements.txt
```

## Usage
Display help message:
```
python -m HeavyOil -h
```

## Command-Line Options
```
Input files:
  -s GRO, --gro GRO          Topology file (.gro)
  -t TRJ, --trj TRJ          Trajectory file (.xtc)

Frame selection:
  -b BEGIN, --begin BEGIN    Start frame index (0-based)
  -e END, --end END          End frame index (exclusive)
  -i INTERVAL, --interval    Frame sampling interval

Cluster information:
  --info                     Perform cluster information analysis
  --infoselect {part,all}    'part' = clustered molecules only  
                             'all' = entire system  

Analysis options:
  --residueselect RESIDUESELECT  
                             Residue selection (default: all)  
                             Example: 'ASP1,ASN2' or 'AS*'
```

# Samples
Run a basic analysis with GRO:
```
python main.py -s system.gro
```
Perform cluster information analysis on the entire system:
```
python main.py -s system.gro -t traj.xtc --info --infoselect all

```
Select specific residues (e.g., ASP1 and ASN2 or AS*) for analysis:
```
python main.py -s system.gro -t traj.xtc --residueselect "ASP1,ASN2"
```
Analyze a specific time window (frames 0–100, every 10th frame):
```
python main.py -s system.gro -t traj.xtc -b 0 -e 100 -i 10
```