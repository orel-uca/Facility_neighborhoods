# Facility Location Optimization Project

This project implements multiple optimization algorithms for facility location problems on graphs with non-convex neighborhoods using mathematical programming and the Gurobi solver. The project focuses on three main types of location problems: **P-Median**, **P-Center** and **P-Maximal_Covering**.

## 📋 Project Description

The project addresses classic facility location problems in a two-dimensional space with Euclidean distance, where demand nodes are represented by neighborhoods that can include:

- **Polygons** (rectangular areas)
- **Disks** (circular areas)
- **Combinations** of both geometry types

### Implemented Problem Types

#### 1. **P-Median**

Minimizes the total sum of weighted distances between customers and their nearest facilities.

#### 2. **P-Center**

Minimizes the maximum distance between any customer and their nearest facility.

#### 3. **P-Maximal_Covering**

Maximizes the demand covered within a specific coverage radius.

## 🗂️ Project Structure

```
Maximal_Facility_neigh/
├── README.md
├── LICENSE
├── bucle_covering.py       # Main script for running P-Covering problems
├── bucle_pmedian.py        # Main script for running P-Median problems
├── bucle_pcenter.py        # Main script for running P-Center problems
├── p_center_v1.py          # P-Center implementation version 1
├── p_center_v2.py          # P-Center implementation version 2
├── p_covering_v1.py        # P-Maximal_Covering implementation version 1
├── p_covering_v2.py        # P-Maximal_Covering implementation version 2
├── p_median_v1.py          # P-Median implementation version 1
├── p_median_v2.py          # P-Median implementation version 2
├── Data_COV/               # Datasets for P-Maximal_Covering problems
│   ├── Demands/            # Demand weight files (.txt)
│   ├── Graphs/             # Geometric and topological data
│   └── Radii/              # Coverage radii data
└── Data_PM-PC/             # Data for P-Median and P-Center problems
    └── Graphs/             # Node and arc data files
        ├── *.csv           # Region geometry files
        ├── *_nodos.txt     # Node data files
        └── 0.07/, 0.15/, 0.20/  # Arc density subdirectories
```

## 🛠️ System Requirements

### Software Dependencies

- **Python 3.7+**
- **Gurobi Optimizer** (license required)
- **Python Libraries:**
  ```
  gurobipy
  numpy
  pandas
  ```

### Installing Dependencies

```bash
pip install gurobipy numpy pandas
```

**Note:** For Gurobi, you will need a valid license. Visit [gurobi.com](https://www.gurobi.com) for information on academic or commercial licenses.

## 🚀 System Usage

### Basic Execution

#### For P-Maximal_Covering Problems:

```bash
python bucle_covering.py
```

#### For P-Median Problems:

```bash
python bucle_pmedian.py
```

#### For P-Center Problems:

```bash
python bucle_pcenter.py
```

### Parameter Configuration

The main scripts allow configuration of various parameters:

```python
TIEMPO = 3600    # Maximum execution time (seconds)
CORTES = -1      # Cuts configuration (0: disabled, -1: default)
HEUR = 0.05      # Heuristic gap (0: disabled, 0.05: default)
PRE = -1         # Preprocessing (-1: default)
N = 2            # Problem dimension (2D)
P = [1, 2, 3, 4] # Number of facilities to locate
```

### Input Data Format

#### Region Files (.csv)

CSV files contain geometric information of regions:

- **Columns [a-b]:** Number of region data including polygons (a) and disks (b)
- **Columns [c-(end)]:**
  - **For rectangles (first columns):** Lower-left and upper-right coordinates
  - **For disks (last columns):** Center coordinates and radius
- **Last row:** Big-M parameter for mathematical programming

#### Demand Files (*_wi.txt)

Text files with demand weight values for each region.

#### Node Files (*_nodos.txt)

Files containing node information for the graph structure.

#### Arc Data

Connectivity between regions is defined through arc density directories (0.07, 0.15, 0.20).

## 📊 Dataset Examples

The project includes multiple test datasets with the following characteristics:

### P-Median and P-Center Data (Data_PM-PC/)

- **16-40 regions:** Different problem sizes (16, 20, 24, 28, 30, 32, 35, 40 nodes)
- **Arc densities:** 0.07, 0.15, 0.20 representing different connectivity levels
- **Multiple instances:** Ej1, Ej2 variations for each size
- **Geometric data:** Region definitions with polygons and disks

### P-Maximal_Covering Data (Data_COV/)

- **16-80 regions:** Broader range including larger instances (up to 80 nodes)
- **Demand patterns:** Multiple demand scenarios (Ej1, Ej2, Ej4)
- **Coverage radii:** Various radius configurations for different coverage scenarios
- **Geometric data:** Region definitions with polygons and disks

## 📝 License

This project is under the GNU GPLv3.0 License. See the `LICENSE` file for more details.

## 👨‍💻 Authors

- **Inmaculada Espejo, Raúl Páez, Justo Puerto, Antonio M. Rodríguez-Chía** (Created in February, 2021)

## 📚 References

This project implements algorithms based on research "Facility location problems on graphs with non-convex neighborhoods". ([https://doi.org/10.1016/j.cor.2023.106356](https://doi.org/10.1016/j.cor.2023.106356)).

## 🆘 Support

If you encounter any problems or have questions:

1. Review carefully the documentation in `README.md`
2. Verify that Gurobi is correctly installed and licensed
3. Ensure that data files are in the correct format
4. Open an issue on GitHub with problem details

---

**Note:** This project requires a valid Gurobi Optimizer license for complete execution.
