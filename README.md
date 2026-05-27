# Description

The code hereafter is made to read in the xsec systematic variations extracted from `nusystematics`. It is composed of several tools in order to do so.

# Input

The input is produced using `nusystematics` by processing CAFs to which a `SystWeights` TTree has been added, recording the weight variations consequent from xsec syst dials variations.
Example of such a file can be found on the dunegpvms at `/exp/dune/data/users/pgranger/nusyst_new_sum.root` or on `lxplus` at `/afs/cern.ch/work/p/pigrange/public/nusyst_new_sum.root`

# extract_genie

`extract_genie.cc` is provide in order to extract and compute the relevant information from the GENIE NtpGHepRecord that is stored in the CAF files under an easy to read flat tree. With GENIE loaded, the following commands should allow to extract these information:

```bash
genie
.L extract_genie.cc
extract_genie("nusyst_new_sum.root", "nusyst_genie_extracted.root")
```
The output is a flat ROOT TTree in `nusyst_genie_extracted.root` that can be read in parallel to the CAFs information (`nusyst_new_sum.root`).
An example of such a file is provided on the dunegpvm at `/exp/dune/data/users/pgranger/nusyst_genie_extracted.root` (matching the `nusyst_new_sum.root` file). On lxplus, this file is located at `/afs/cern.ch/work/p/pigrange/public/nusyst_genie_extracted.root`

# read_systs.ipynb

This notebook shows how to make some example plots with the available data and save them.

# QE pre/post-FSI matching

To study FSI impact on QE events, `CAFManager` now provides a particle-matching utility between pre-FSI and post-FSI hadrons.

```python
from systs import CAFManager

cafs = CAFManager("nusyst_new_sum.root", genie="nusyst_genie_extracted.root")
cafs.add_qe_fsi_matching(rel_energy_weight=0.2)

# Event-level columns (added to cafs.data):
# - preFSI_match_post_index : list[int], matched post index for each pre particle (-1 if unmatched)
# - preFSI_delta_theta      : list[float], angle(pre, post) in radians
# - preFSI_delta_E          : list[float], postE - preE in GeV
# - prefsi_n_matched        : number of matched pre-FSI particles
# - prefsi_frac_matched     : matched fraction in each QE event

# Long-format table (one row per matched particle pair in QE events)
fsi_pairs = cafs.qe_fsi_match_long_table()
```

The matching is one-to-one, constrained to equal PDG, and uses a kinematic score combining angular distance and relative energy difference.

# Running on a Remote Cluster

To run the `read_systs.ipynb` notebook and associated scripts on a remote cluster (e.g. DUNE GPVMs or CERN lxplus), follow these steps to set up the repository, dependencies, and data files.

## 1. Clone the Repository

Clone this repository on your remote cluster:
```bash
git clone <your-git-repository-url>
cd systematics
```

## 2. Set Up the Environment

You can configure the Python environment using any of the following methods. **Conda/Mamba** is recommended on shared clusters like lxplus or GPVMs.

### Option A: Conda/Mamba (Recommended for clusters)
Conda is highly stable on scientific computing nodes. To create and activate the environment:
```bash
conda env create -f environment.yml
conda activate systematics
```

### Option B: uv (Fastest Python virtual environment tool)
If `uv` is installed, you can create the virtual environment and install all dependencies in seconds:
```bash
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

### Option C: Standard Python venv & pip
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## 3. Link or Copy the ROOT Data Files

Because large `.root` files and the cached `parquet/` directory are excluded from Git via `.gitignore`, you need to set up the data files manually on the remote cluster.

You can copy them from your local machine, or if you are on **DUNE GPVM** or **CERN lxplus**, you can symlink the pre-staged files:

### On DUNE GPVMs:
```bash
# Symlink the main CAF and GENIE-extracted files
ln -s /exp/dune/data/users/pgranger/nusyst_new_sum.root merged_new_systs.root
ln -s /exp/dune/data/users/pgranger/nusyst_genie_extracted.root genie_extracted_new_cafs.root
```

### On CERN lxplus:
```bash
# Symlink the main CAF and GENIE-extracted files
ln -s /afs/cern.ch/work/p/pigrange/public/nusyst_new_sum.root merged_new_systs.root
ln -s /afs/cern.ch/work/p/pigrange/public/nusyst_genie_extracted.root genie_extracted_new_cafs.root
```

## 4. Run the Jupyter Notebook

### Via VS Code Remote-SSH (easiest)
1. Open VS Code and connect to your remote cluster via the **Remote - SSH** extension.
2. Navigate to this folder and open `read_systs.ipynb`.
3. In the top-right of the notebook, click **Select Kernel** and choose the Python environment you created in Step 2 (`systematics` conda env or `.venv` virtual environment).

### Via SSH Port Forwarding
1. Start the Jupyter server on the remote node:
   ```bash
   jupyter notebook --no-browser --port=8888
   ```
2. On your local machine, run SSH port forwarding to map the remote Jupyter port to your localhost:
   ```bash
   ssh -N -f -L localhost:8888:localhost:8888 username@remote-cluster-address
   ```
3. Open your local browser and navigate to the URL printed by the remote Jupyter command (typically `http://localhost:8888/?token=...`).