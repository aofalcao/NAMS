"""Double-bond E/Z assignment for the Open Babel 3 NAMS reference backend.

The original NAMS module generated 2D coordinates and inferred cis/trans from
geometry.  Open Babel 3 already stores cis/trans stereochemistry as structured
stereo data, so this rewrite uses ``OBStereoFacade`` and retains the historical
NAMS result convention:

* ``+1``: Z (higher-priority substituents cis)
* ``-1``: E (higher-priority substituents trans)
* ``0``: unspecified or unresolved

Ligand priority uses the legacy NAMS breadth-first signature rather than claiming
full modern CIP coverage.  This is deterministic and substantially less fragile
than regenerating 2D coordinates.
"""

from __future__ import annotations

from typing import Iterable

try:
    from openbabel import openbabel, pybel
except ImportError as exc:  # pragma: no cover
    openbabel = None  # type: ignore[assignment]
    pybel = None  # type: ignore[assignment]
    _OPENBABEL_IMPORT_ERROR: ImportError | None = exc
else:
    _OPENBABEL_IMPORT_ERROR = None

from chirality import ligand_signature


def require_openbabel() -> None:
    if _OPENBABEL_IMPORT_ERROR is not None:
        raise RuntimeError(
            "Open Babel Python bindings are required. Install Open Babel 3 and "
            "verify: 'from openbabel import openbabel, pybel'."
        ) from _OPENBABEL_IMPORT_ERROR


def _normalise_input_format(input_type: str) -> str:
    return "smi" if input_type.lower() == "can" else input_type.lower()


def _implicit_ref() -> int:
    require_openbabel()
    try:
        return int(openbabel.OBStereo.ImplicitRef)
    except (AttributeError, TypeError):
        # Value documented by Open Babel's Python stereo guide.
        return 0xFFFFFFFE


def _eligible_double_bond(bond: object) -> bool:
    """Return whether *bond* is a NAMS-compatible C=C stereo bond.

    Open Babel 3's ``OBBond`` API exposes ``GetBondOrder()`` but does not
    define ``IsDouble()``.  Similarly, carbon identity is obtained from
    ``OBAtom.GetAtomicNum()`` rather than the non-portable ``IsCarbon()``
    convenience method found in some older bindings.
    """

    begin = bond.GetBeginAtom()
    end = bond.GetEndAtom()

    if begin is None or end is None:
        return False

    return bool(
        int(bond.GetBondOrder()) == 2
        and not bool(bond.IsAromatic())
        and int(begin.GetAtomicNum()) == 6
        and int(end.GetAtomicNum()) == 6
        and int(begin.CountBondsOfOrder(2)) == 1
        and int(end.CountBondsOfOrder(2)) == 1
    )


def _side_ligands(atom: object, other: object) -> list[int]:
    result = [
        int(neighbour.GetId())
        for neighbour in openbabel.OBAtomAtomIter(atom)
        if int(neighbour.GetIdx()) != int(other.GetIdx())
    ]
    while len(result) < 2:
        result.append(_implicit_ref())
    return result[:2]


def _signature_for_ref(mol: object, center: object, ref: int) -> tuple:
    if ref == _implicit_ref():
        return ((1.0,),)
    ligand = mol.GetAtomById(ref)
    if ligand is None:
        return tuple()
    return ligand_signature(
        mol,
        int(center.GetIdx()),
        int(ligand.GetIdx()),
        isotopes=False,
    )


def _highest_priority_ref(mol: object, center: object, refs: list[int]) -> int | None:
    ranked = [(_signature_for_ref(mol, center, ref), ref) for ref in refs]
    if ranked[0][0] == ranked[1][0]:
        return None
    return max(ranked, key=lambda item: item[0])[1]


class Stereodoubleb:
    """Expose NAMS E/Z values by Open Babel bond index."""

    def __init__(self, line_notation: str, input_type: str = "smi") -> None:
        require_openbabel()
        source = pybel.readstring(_normalise_input_format(input_type), line_notation)
        self.can_smi = source.write("can").strip()
        canonical = pybel.readstring("smi", self.can_smi)
        self.obmol = canonical.OBMol
        self.n_bonds = int(self.obmol.NumBonds())
        self.n_atoms = int(self.obmol.NumAtoms())
        self.e_z: list[int] = [0] * self.n_bonds

        facade = openbabel.OBStereoFacade(self.obmol)
        for bond in openbabel.OBMolBondIter(self.obmol):
            if not _eligible_double_bond(bond):
                continue
            bond_id = int(bond.GetId())
            if not facade.HasCisTransStereo(bond_id):
                continue
            stereo = facade.GetCisTransStereo(bond_id)
            if stereo is None or not stereo.IsSpecified():
                continue

            begin = bond.GetBeginAtom()
            end = bond.GetEndAtom()
            begin_ref = _highest_priority_ref(
                self.obmol, begin, _side_ligands(begin, end)
            )
            end_ref = _highest_priority_ref(
                self.obmol, end, _side_ligands(end, begin)
            )
            if begin_ref is None or end_ref is None:
                continue

            try:
                if stereo.IsCis(begin_ref, end_ref):
                    self.e_z[int(bond.GetIdx())] = 1
                elif stereo.IsTrans(begin_ref, end_ref):
                    self.e_z[int(bond.GetIdx())] = -1
            except (AttributeError, TypeError):
                # Some SWIG builds expose GetCisRef/GetTransRef but not IsCis.
                try:
                    if int(stereo.GetCisRef(begin_ref)) == end_ref:
                        self.e_z[int(bond.GetIdx())] = 1
                    elif int(stereo.GetTransRef(begin_ref)) == end_ref:
                        self.e_z[int(bond.GetIdx())] = -1
                except (AttributeError, TypeError, ValueError):
                    self.e_z[int(bond.GetIdx())] = 0

    def get_e_z_at(self, at1: object, at2: object) -> int:
        bond = at2.GetBond(at1)
        if bond is None:
            raise ValueError("the supplied atoms are not bonded")
        return self.get_e_z_bond(int(bond.GetIdx()))

    def get_e_z_bond(self, bond_idx: int) -> int:
        if bond_idx < 0 or bond_idx >= len(self.e_z):
            raise IndexError(f"bond index out of range: {bond_idx}")
        return self.e_z[bond_idx]
