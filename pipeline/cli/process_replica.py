import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import torch
import esm
import MDAnalysis as mda
from pipeline.config import ESM_MODEL_NAME
from pipeline.features.esm_embed import extract_features
from pipeline.graph.builder import build_pyg_graph, write_contact_sidecar
from pipeline.sim.processor import SimulationProcessor
LOGGER = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica-dir", type=Path, required=False, default="../mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1")
    parser.add_argument("--protein-id", type=str, required=False, default="AFR1")
    parser.add_argument("--pocket-id", type=str, required=False, default="pocket1")
    parser.add_argument("--ligand", type=str, default="Milbemycin")

    args = parser.parse_args()

    cutoff = None
    OUT_ROOT = Path("../pipeline/data_graph_v3").resolve()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # ───────────────────── run the processor ─────────────────────
    processor = SimulationProcessor(
        replica_dir        = args.replica_dir,
        protein_id         = args.protein_id,
        pocket_id          = args.pocket_id,
        ligand_name        = args.ligand,
        stability_threshold=cutoff,        # ← derived above (or None)
        logger             = LOGGER,
        target_n           = 20,
    )

    sim_features = processor.run()

    # --- PASS 2 Logic ---
    if not sim_features:
        LOGGER.error(
            f"Feature extraction failed for pocket {args.pocket_id}. Aborting."
        )
        return
    ligand_drift = sim_features.get("ligand_drift")

    if ligand_drift is None or np.isnan(ligand_drift):
        LOGGER.warning(f"Skipping graph for {args.protein_id}:{args.pocket_id} — missing ligand_drift")
        return None  # Signal to caller: don't use this graph

    frames_dir = args.replica_dir / "frames"
    LOGGER.info(f"Setting up ESM models for {args.pocket_id}...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL_NAME)
    batch_converter = alphabet.get_batch_converter()
    model.to(torch.device("cuda" if torch.cuda.is_available() else "cpu"))
    model.eval()

    frame_files = sorted(frames_dir.glob("frame_*.pdb"))
    N_FRAMES = len(frame_files) or 20
    for frame_file in frame_files:
        try:
            match = re.search(r"frame_(\d+)\.pdb", frame_file.name)
            if not match:
                LOGGER.warning(f"Skipping malformed frame filename: {frame_file.name}")
                continue
            i = int(match.group(1))
            graph_path = (OUT_ROOT  / args.protein_id/ args.pocket_id/ f"{args.replica_dir.name}"/ f"frame_{i:04d}_graph.pt")
            plip_csv_path = Path(args.replica_dir) /"plip_results/all_plip_interactions_summary.csv"
            if not graph_path.exists():
                graph_path.parent.mkdir(parents=True, exist_ok=True)
                x = extract_features(frame_file, model, batch_converter, graph_path)

                frame_data = {
                    "x": x,
                    "ca_coords": [sim_features["ca_coords"][i]],
                    "plip_ref_pdb": Path(args.replica_dir) / f"{args.protein_id}_prepared_{args.ligand}_{args.pocket_id}_complex_recombined_complex_explicit_stripped_initial_frame.pdb",
                    "plip_csv_path": plip_csv_path,

                    # drift only — binary label dropped
                    "y_data": [sim_features["pocket_rmsd_per_frame"][i]],

                    "rmsd_drift": sim_features["pocket_rmsd_per_frame"][i],
                    "energy_variance": sim_features.get("energy_variance"),
                    "gate_aperture_median": sim_features.get("gate_aperture_median"),
                    "nbd_distance_median": sim_features.get("nbd_distance_median"),
                    "cavity_volume": sim_features.get("cavity_volume"),
                    "protein_id": args.protein_id,
                    "pocket_id":  args.pocket_id,
                    "replica_id": args.replica_dir.name,
                    "ligand_drift": sim_features.get("ligand_drift"),
                    "frame_idx":  i,
                    "frames_per_replica": N_FRAMES,
                }

                g = build_pyg_graph(frame_data)
                if g is None or not hasattr(g, "x") or g.x is None:
                    LOGGER.error(f"Graph is empty or invalid for frame {i}. Skipping save.")
                    continue
                torch.save(g, graph_path)
                LOGGER.info(f"Saved graph {graph_path}")
            LOGGER.info(f"Generating contact.graph.pt ")
            # Graph already built – just ensure side-car exists
            data = torch.load(graph_path, map_location="cpu")
            contact_path = graph_path.with_name(graph_path.name.replace("_graph.pt", "_graph.contact.pt"))

            if not contact_path.exists():
               # contact_path.unlink()
                #LOGGER.info(f"Delete contact.pt : {contact_path.name}")
                write_contact_sidecar(
                    data,
                    graph_path,
                    Path(args.replica_dir) /"plip_results/all_plip_interactions_summary.csv"
                )

        except Exception as e:
            LOGGER.error(
                f"Failed to build graph for frame {frame_file} of {args.pocket_id}: {e}",
                exc_info=True,
            )

    # ----- done with all frames ------------------------------------------------
    del model, batch_converter        # free big ESM tensors
    torch.cuda.empty_cache()          # release GPU blocks
    import gc, os
    gc.collect()                      # flush Python refs
    LOGGER.info("All tensors freed, exiting cleanly")

if __name__ == "__main__":
    main()

#python -m pipeline.cli.process_replica --replica-dir "/home/zhenli/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1"  --protein-id "AFR1" --pocket-id "pocket1" --ligand "Milbemycin"