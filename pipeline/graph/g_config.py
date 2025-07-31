# Globally defined or passed in at runtime (for reproducibility)
from __future__ import annotations

import numpy as np

# Core pocket physics (rmsd_drift … ligand_drift)	7
# contact_count (our new slot 7)	+1
# PLIP global counts
#7 contact types × 2 phases (initial + final)	+14
# domain map 4

EXPECTED_U_RAW_DIM = 26
NODE_FEATURE_DIM_STATIC = 1308

#Each node $i$ is represented by a 1,309-dimensional feature vector $x_i \in \mathbb{R}^{1309}$ composed of: a 1,280-dimensional embedding from ESM-2, 
# scalar features for hydrophobicity, B-factor, solvent accessibility (SASA), and sidechain volume, a 3-dimensional one-hot encoding of secondary structure class, a 21-dimensional PLIP fingerprint spanning initial, final, 
# and $\Delta$ contact counts for seven interaction types, and a 1-dimensional energy signal replicated from the global MD descriptor. This final energy value is attached at load time.
# x_i = [ESM-2 | hydro | bfactor | sasa | ss | volume | PLIP | energy_i]
#       = 1280 + 1 + 1 + 1 + 3 + 1 + 21 + 1 = 1309
NODE_FEATURE_DIM= 1309
NEIGHBOR_CUTOFF = 10.0
CONTACT_CUTOFF = 4.0
# the engener_feature was added by dataload on the fly, the final feature is 1309
Y_SCALE = 1.0

MAX_PROTEINS = 20                      # global upper-limit
PROTEIN_LIST = [
    "AFR1",
    "ATRF_ASPFU",
    "CDR1_CANAR",
    "CDR1_CANAR_auris",
    "CDR2_CANAL",
    "CIMG_00533",
    "CIMG_06197",
    "MDR1_CRYNH",
    "MDR1_TRIRC",
    "MDR2_TRIRC",
    "PDH1_CANGA",
    "SNQ2_CANGA",
    "PDR5_YEAST",
    "CIMG_00780",
    "CIMG_01418",
    "CIMG_09093"
]
_pid_lookup = {pid: i for i, pid in enumerate(PROTEIN_LIST)}
U_RAW_FEATURES = {
    "rmsd_drift": 0,
    "energy_variance": 1,
    "gate_aperture_median": 2,
    "nbd_distance_median": 3,
    "cavity_volume": 4,
    "pocket_rmsd_90th_percentile": 5, #not used
    "ligand_drift": 6,
    "contact_count": 7
}
THRESHOLDS = {
    "rmsd_drift": (0.0, 10.0),
    "energy_variance": (0.0, 2e3),
    "gate_aperture_median": (0.0, 50.0),
    "nbd_distance_median": (0.0, 100.0),
    "cavity_volume": (0.0, 5000.0),
    "pocket_rmsd_90th_percentile": (0.0, 10.0),
    "ligand_drift": (0.0, 10.0),
    "contact_count": (0.0, 100.0)
}
DOMAIN_START = 1308
# Add to g_config.py
#INDEX DROPPING FEATURES
D_ESM   = 1280
D_PHYS  = 7
D_PLIP  = 21
D_FULL  = D_ESM + D_PHYS + D_PLIP    # 1308

FULL_IDX      = np.arange(D_FULL)
NO_PLIP_IDX   = FULL_IDX[:D_ESM + D_PHYS]                    # drop last 21
NO_ESM_IDX    = FULL_IDX[D_ESM:]                             # keep phys+plip
PHYS_ONLY_IDX = FULL_IDX[D_ESM:D_ESM + D_PHYS]               # 7 dims
NO_ESM_PLIP_IDX = FULL_IDX[D_ESM:D_ESM + D_PHYS]

DOMAIN_FLAGS = ["NBD", "TMD", "GATE", "LOOP"]   # extend if needed
FLAG_DIM      = len(DOMAIN_FLAGS)

