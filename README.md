# QdoRNAR - Genetic Circuit Modeling for Flavonoid-Responsive Systems

This repository contains a mathematical modeling framework for simulating genetic circuits, specifically focusing on the QdoR regulatory system responsive to flavonoids (particularly quercetin). The model uses mechanistic ordinary differential equations (ODEs) to predict the behavior of biosensor circuits in bacterial cells.

## Overview

The QdoRNAR project implements a deterministic mathematical model to simulate the dynamics of a genetic regulatory network where:
- **QdoR** is a repressor protein that regulates gene expression
- **Quercetin** (a flavonoid) acts as an inducer by binding to QdoR
- **GFP** (Green Fluorescent Protein) is used as a biomarker to measure circuit activity

The model accounts for:
- Protein formation and degradation kinetics
- Gene-repressor binding dynamics
- Cell volume growth over time
- Multiple promoter strengths (Anderson promoter library)

## Project Structure

```
QdoRNAR/
├── QDOR-DFFM-Versão final/
│   └── QDOR-DFFM-Versão final/
│       ├── Functions/           # Mathematical functions and ODE solvers
│       │   ├── functions cubic version.sce
│       │   └── functions cubic version tests.sce
│       ├── Input/              # Experimental data files
│       │   ├── Input.csv
│       │   ├── OD600 input.xlsx
│       │   └── *.csv (various experimental datasets)
│       └── Lib/                # Main execution scripts and parameters
│           ├── loopExe.sce     # Main entry point
│           ├── Main.sce        # Core simulation and plotting
│           ├── Parameters.sce  # Model parameters
│           └── *.csv (output files)
├── README.md
└── LICENSE
```

### Key Files

- **`loopExe.sce`**: Main entry point that imports data, sets parameters, and runs simulations across different parameter ranges
- **`Main.sce`**: Core simulation script that solves ODEs and generates plots
- **`Parameters.sce`**: Defines biological parameters (kinetic constants, dissociation constants, etc.)
- **`functions cubic version.sce`**: Contains the differential equation solvers and mathematical functions

## Requirements

- **Scilab** (version 6.0 or higher recommended)
  - Download from: https://www.scilab.org/download/

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/Harveker/QdoRNAR.git
   cd QdoRNAR
   ```

2. Ensure Scilab is installed on your system

3. No additional dependencies required - all necessary files are included

## How to Use

### Running the Simulation

1. Open Scilab

2. Navigate to the main script directory:
   ```scilab
   cd("path/to/QdoRNAR/QDOR-DFFM-Versão final/QDOR-DFFM-Versão final/Lib")
   ```

3. Execute the main script:
   ```scilab
   exec("loopExe.sce", -1)
   ```

### What the Simulation Does

The simulation will:
1. Import experimental data from CSV files
2. Set up initial conditions (cell volume, initial protein concentrations)
3. Solve the system of differential equations over time (0-900 minutes)
4. Generate plots showing:
   - **Cell volume behavior** over time
   - **GFP fluorescence** for different promoter strengths
   - **QdoR repressor concentration** dynamics
5. Export results to CSV files

### Understanding the Output

After running, you'll see three figure windows:
- **Figure 0**: Cell volume growth over time
- **Figure 1**: GFP fluorescence curves for 6 different conditions
- **Figure 2**: QdoR repressor protein curves

Output CSV files will be generated in the `Lib` directory:
- `rep*.csv`: Repressor concentration data
- `flu*.csv`: Fluorescence intensity data
- `vol*.csv`: Volume data
- `time.csv`: Time points

## Model Parameters

Key biological parameters (defined in `Parameters.sce` and `loopExe.sce`):

| Parameter | Description | Value | Units |
|-----------|-------------|-------|-------|
| k_fqdoR | QdoR formation rate | 6.60 | min⁻¹ |
| k_fgfp | GFP formation rate | 11.04 | min⁻¹ |
| k_dR | QdoR degradation rate | 3×10⁻⁵ | min⁻¹ |
| k_dgfp | GFP degradation rate | 5×10⁻⁴ | min⁻¹ |
| [Q] | Quercetin concentration | 2.5×10⁻¹² | mol L⁻¹ |
| K_dQ | Q~R dissociation constant | 1×10⁻⁴ | M |
| K_dR | Promoter~R dissociation | 4×10⁻⁹ | M |
| r | Volume increase constant | 0.003 | - |

### Modifying Parameters

To adjust parameters for your experiments:

1. Open `loopExe.sce`
2. Locate the parameter section (around line 91-105)
3. Modify values in the `pz()` array
4. Save and re-run the script

Example:
```scilab
pz(5) = 5D-12    // Change quercetin concentration
pz(12) = 0.005   // Change growth rate constant
```

## Mathematical Model

The system is described by three main differential equations:

1. **Cell Volume**: 
   ```
   dV/dt = r·V·(1 - V/k)
   ```

2. **Total Repressor Concentration**:
   ```
   d[R]_total/dt = k_fR·vr·[R_gene]_free - k_dR·[R]_free - ([R]_total/V)·(dV/dt)
   ```

3. **GFP Concentration**:
   ```
   d[GFP]/dt = k_fB·[BM_gene]_free - k_dB·[GFP] - ([GFP]/V)·(dV/dt)
   ```

Free protein concentration is calculated using a cubic equation solver (Cardano-Tartaglia method) that accounts for multiple binding equilibria.

## Experimental Data

The `Input` directory contains experimental fluorescence data for different genetic constructs:
- **QPG114, QPG115, QPG116**: Different genetic circuit variants
- **QPG105, QPG110**: Additional constructs with varying promoter configurations

These datasets are used to validate and calibrate the model.

## Customization

### Running Parameter Sweeps

The `loopExe.sce` script includes a loop for parameter variation:

```scilab
h = 0.0023:0.0005:0.0063    // Range for parameter sweep
for i = h
    pz(12) = i              // Parameter to vary
    exec("Main.sce", [-1])
end
```

Modify:
- `h`: The range and step size for parameter variation
- `pz(12)`: Which parameter to vary (change the index)

### Changing Simulation Time

In `loopExe.sce`, modify:
```scilab
t = 1:900/91:900    // Time vector: start:step:end (minutes)
```

## Troubleshooting

**Issue**: Script fails to find files
- **Solution**: Ensure you run `loopExe.sce` from within the `Lib` directory, or use full paths

**Issue**: Plots are not displayed
- **Solution**: Check that your Scilab graphics system is properly configured

**Issue**: "more than 1000 iterations" message
- **Solution**: The numerical solver is having difficulty converging. Try adjusting initial conditions or parameter values

## Contributing

Contributions are welcome! Areas for improvement:
- Adding more genetic circuit variants
- Implementing stochastic simulations
- Improving parameter estimation algorithms
- Adding experimental data fitting tools

To contribute:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## References

The model is based on published research on:
- Flavonoid-responsive genetic circuits
- Mechanistic modeling of gene regulation
- QdoR repressor systems

For questions or collaboration, contact:
- Franco Endrigo
- franco.endrigo.r@gmail.com
- francoendrigo@alunos.utfpr.edu.br

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

**Note**: This is a research tool for scientific modeling. Results should be validated against experimental data before drawing biological conclusions.
