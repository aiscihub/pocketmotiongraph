import pandas as pd
from pathlib import Path
import torch

class ContactLookup:
    def __init__(self, csv_path: Path):
        df = pd.read_csv(csv_path)
        # residue_id ↔ contact bool per protein
        self._tbl = (df.assign(is_contact=~df['protein'].isna())
                     .groupby(['protein', 'residue_id'])['is_contact']
                     .any())
        self.contacts = set(
            (row['protein'], int(row['residue_id']))
            for _, row in df.iterrows()
        )

    def label(self, prot_id: str, residue_ids: torch.Tensor | list[int]) -> int:
        """
        Returns 1 iff *any* residue in `residue_ids` is annotated as contact for
        `prot_id`, else 0.
        """
        # self.mask(...) must already return a Bool / 0-1 tensor of same length
        return int(self.mask(prot_id, residue_ids).any())

    def mask(self, protein_id: str, residue_ids: torch.Tensor) -> torch.Tensor:
        """
        Return a 0/1 tensor (same length as `residue_ids`) indicating whether
        each residue participates in a PLIP contact for the given protein.
        """
        contacts = self._tbl[protein_id]            #

        return torch.tensor([int(r.item() in contacts) for r in residue_ids],
                            device=residue_ids.device)

    def __call__(self, prot: str, rid_tensor):
        # rid_tensor: (N,) int32 on CPU
        key = pd.MultiIndex.from_product([[prot], rid_tensor.cpu().numpy()],
                                         names=['protein','residue_id'])
        return torch.tensor(self._tbl.reindex(key, fill_value=False).values,
                            dtype=torch.float32, device=rid_tensor.device)

    def is_contact(self, prot_id, residue_id):
        return int((prot_id, int(residue_id)) in self.contacts)
