import json
from pathlib import Path

DOMAIN_MAP_DIR = Path("/home/zhenli/git/valleyfevermutation/mutation_pipeline/dataset/proteins/gff/")
# coarse labels you already use in the graph builder
DOMAIN_GROUPS = {
    "NBD":  "NBD",
    "TM":   "TMD",
    "GATE": "GATE",
    "LOOP": "LOOP",
}

def parse_gff(gff_path: Path) -> dict[int, str]:
    """Parse GFF file into residue-to-domain mapping"""
    residue_map = {}

    if not gff_path.exists():
        return residue_map

    with gff_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            # Handle different column formats
            if len(parts) >= 5:
                feature, start, end = parts[2], int(parts[3]), int(parts[4])
            else:
                feature, start, end = parts[1], int(parts[2]), int(parts[3])

            # Determine coarse domain
            f_lower = feature.lower()
            if "tmhelix" in f_lower:
                domain = DOMAIN_GROUPS["TM"]
            elif feature in {"inside", "outside"}:
                domain = DOMAIN_GROUPS["LOOP"]
            elif "nbd" in f_lower:
                domain = DOMAIN_GROUPS["NBD"]
            elif "gate" in f_lower:
                domain = DOMAIN_GROUPS["GATE"]
            else:
                continue

            # Add residues to map
            for resid in range(start, end + 1):
                residue_map[resid] = domain

    return residue_map

def generate_all_domain_maps(gff_dir: Path, output_dir: Path):
    """Generate domain maps for all proteins and save as JSON files"""
    output_dir.mkdir(parents=True, exist_ok=True)

    for gff_path in gff_dir.glob("*.gff3"):
        protein_id = gff_path.stem
        residue_map = parse_gff(gff_path)

        # Save as JSON
        output_path = output_dir / f"{protein_id}_domain_map.json"
        with output_path.open("w") as f:
            json.dump(residue_map, f)

        print(f"Generated: {output_path}")

def make_residue_domain_map(protein_id: str) -> dict[int, str]:
    """
    Load residue-to-domain mapping from precomputed JSON file
    Returns empty dict if not found
    """
    map_path = DOMAIN_MAP_DIR / f"{protein_id}_domain_map.json"

    if not map_path.exists():
        return {}

    with map_path.open() as f:
        # JSON keys are strings, convert to integers
        str_map = json.load(f)
        return {int(k): v for k, v in str_map.items()}

_DOMAIN_MAP_CACHE = {}

def cached_domain_map(protein_id: str) -> dict[int, str]:
    """Get domain map with caching"""
    if protein_id not in _DOMAIN_MAP_CACHE:
        _DOMAIN_MAP_CACHE[protein_id] = make_residue_domain_map(protein_id)
    return _DOMAIN_MAP_CACHE[protein_id]

if __name__ == "__main__":
    gff_directory = Path("/home/zhenli/git/valleyfevermutation/mutation_pipeline/dataset/proteins/gff/")
    output_directory = Path("/home/zhenli/git/valleyfevermutation/mutation_pipeline/dataset/proteins/gff")
    generate_all_domain_maps(gff_directory, output_directory)
