"""Legacy-compatible tetrahedral chirality support for NAMS 2.

This module is a Python 3 / Open Babel 3 rewrite of the original NAMS
``chirality.py``.  It keeps the public ``Chirality`` class and its zero-based
``get_chirality()`` method, but replaces Python 2 idioms and recursive mutable
container handling with deterministic Python 3 code.

The returned NAMS convention is:

* ``+1``: R-like orientation under the legacy NAMS priority algorithm
* ``-1``: S-like orientation under the legacy NAMS priority algorithm
* ``0``: unspecified, non-stereogenic, or unresolved due to tied ligands

This is intentionally a *reference compatibility algorithm*, not a complete
implementation of all modern CIP sequence rules.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Iterable, Sequence

try:
    from openbabel import openbabel, pybel
except ImportError as exc:  # pragma: no cover - exercised on installations without OB
    openbabel = None  # type: ignore[assignment]
    pybel = None  # type: ignore[assignment]
    _OPENBABEL_IMPORT_ERROR: ImportError | None = exc
else:
    _OPENBABEL_IMPORT_ERROR = None


IMPLICIT_REF_FALLBACK = 0xFFFFFFFE


def require_openbabel() -> None:
    """Raise a useful error when the Open Babel Python bindings are absent."""

    if _OPENBABEL_IMPORT_ERROR is not None:
        raise RuntimeError(
            "Open Babel Python bindings are required. Install Open Babel 3 and "
            "verify: 'from openbabel import openbabel, pybel'."
        ) from _OPENBABEL_IMPORT_ERROR


def _normalise_input_format(input_type: str) -> str:
    # ``can`` is an output format name; its content is still SMILES.
    return "smi" if input_type.lower() == "can" else input_type.lower()


def scan_smiles_rotations(canonical_smiles: str) -> list[int]:
    """Return the legacy NAMS ``@``/``@@`` orientation sequence.

    The historical implementation scans the canonical SMILES text and maps
    ``@`` to -1 and ``@@`` to +1.  We retain that convention because it is the
    reference point for old NAMS databases.
    """

    rotations: list[int] = []
    index = 0
    while index < len(canonical_smiles):
        if canonical_smiles[index] != "@":
            index += 1
            continue
        if index + 1 < len(canonical_smiles) and canonical_smiles[index + 1] == "@":
            rotations.append(1)
            index += 2
        else:
            rotations.append(-1)
            index += 1
    return rotations


def permutation_parity_sign(permutation: Sequence[int]) -> int:
    """Return +1 for an even permutation and -1 for an odd permutation."""

    inversions = 0
    for left in range(len(permutation)):
        for right in range(left + 1, len(permutation)):
            if permutation[left] > permutation[right]:
                inversions += 1
    return -1 if inversions % 2 else 1


def _bond_multiplicity(bond: object) -> int:
    """Return the legacy integer multiplicity used in ligand signatures."""

    # Open Babel 2 exposed GetBO(); Open Babel 3 normally exposes
    # GetBondOrder(). Keep both for compatibility with binary distributions.
    if hasattr(bond, "GetBondOrder"):
        value = int(getattr(bond, "GetBondOrder")())
    elif hasattr(bond, "GetBO"):
        value = int(getattr(bond, "GetBO")())
    else:  # defensive; a normal OBBond always provides one of these
        value = 1
    return max(1, value)


def _atom_value(atom: object, isotopes: bool) -> float:
    if isotopes:
        return float(getattr(atom, "GetExactMass")())
    return float(getattr(atom, "GetAtomicNum")())


def ligand_signature(
    mol: object,
    center_idx: int,
    ligand_idx: int,
    *,
    isotopes: bool = False,
) -> tuple[tuple[float, ...], ...]:
    """Build the historical NAMS breadth-first ligand signature.

    Atom indices are Open Babel's one-based ``Idx`` values.  The bond from the
    stereocentre to the ligand is excluded. Multiple bonds duplicate the
    destination atomic number according to bond order, matching the original
    algorithm's intent.
    """

    require_openbabel()
    ligand = mol.GetAtom(ligand_idx)
    if ligand is None:
        return tuple()

    levels: list[tuple[float, ...]] = [(_atom_value(ligand, isotopes),)]
    blocked = frozenset((center_idx, ligand_idx))
    visited_edges: set[frozenset[int]] = {blocked}
    frontier: list[int] = [ligand_idx]

    while frontier:
        next_frontier: list[int] = []
        next_seen: set[int] = set()
        level_values: list[float] = []

        for atom_idx in frontier:
            atom = mol.GetAtom(atom_idx)
            if atom is None:
                continue
            for neighbour in openbabel.OBAtomAtomIter(atom):
                neighbour_idx = int(neighbour.GetIdx())
                edge = frozenset((atom_idx, neighbour_idx))
                if edge in visited_edges:
                    continue
                visited_edges.add(edge)

                bond = atom.GetBond(neighbour)
                multiplicity = _bond_multiplicity(bond)
                level_values.extend(
                    [_atom_value(neighbour, isotopes)] * multiplicity
                )
                if neighbour_idx not in next_seen:
                    next_seen.add(neighbour_idx)
                    next_frontier.append(neighbour_idx)

        if not level_values:
            break
        level_values.sort(reverse=True)
        levels.append(tuple(level_values))
        frontier = next_frontier

    return tuple(levels)


class Chirality:
    """Determine legacy NAMS tetrahedral atom chirality using Open Babel 3."""

    def __init__(
        self,
        line_notation: str,
        input_type: str = "smi",
        isotopes: bool = False,
    ) -> None:
        require_openbabel()
        self.obmol = None
        self.n_atoms = 0
        self.n_bonds = 0
        self.can_smi = ""
        self.chiralities: list[int] = []
        self._isotopes = isotopes

        source = pybel.readstring(_normalise_input_format(input_type), line_notation)
        self.can_smi = source.write("can").strip()
        canonical = pybel.readstring("smi", self.can_smi)
        self.obmol = canonical.OBMol
        self.obmol.AddHydrogens()
        self.n_atoms = int(self.obmol.NumAtoms())
        self.n_bonds = int(self.obmol.NumBonds())

        rotations = scan_smiles_rotations(self.can_smi)
        chiral_atoms = [
            atom
            for atom in openbabel.OBMolAtomIter(self.obmol)
            if bool(atom.IsChiral())
        ]

        # The historical implementation deliberately discarded all chirality
        # when the SMILES did not specify every perceived stereocentre.
        if len(rotations) != len(chiral_atoms):
            self.chiralities = [0] * self.n_atoms
            return

        rotation_by_idx = {
            int(atom.GetIdx()): rotations[position]
            for position, atom in enumerate(chiral_atoms)
        }

        for atom in openbabel.OBMolAtomIter(self.obmol):
            atom_idx = int(atom.GetIdx())
            rotation = rotation_by_idx.get(atom_idx)
            if rotation is None:
                self.chiralities.append(0)
                continue
            self.chiralities.append(self._calculate_atom(atom, rotation))

    def _calculate_atom(self, atom: object, rotation: int) -> int:
        neighbours = list(openbabel.OBAtomAtomIter(atom))
        if len(neighbours) != 4:
            return 0

        center_idx = int(atom.GetIdx())
        signatures: list[tuple[tuple[tuple[float, ...], ...], int]] = []
        for original_position, neighbour in enumerate(neighbours):
            signature = ligand_signature(
                self.obmol,
                center_idx,
                int(neighbour.GetIdx()),
                isotopes=self._isotopes,
            )
            signatures.append((signature, original_position))

        only_signatures = [entry[0] for entry in signatures]
        if len(set(only_signatures)) != len(only_signatures):
            return 0

        # Highest-priority ligand first, as in the original reverse-sort.
        ordered = sorted(signatures, reverse=True)
        permutation = [original_position for _, original_position in ordered]
        return rotation * permutation_parity_sign(permutation)

    def get_chirality(self, atom_id: int) -> int:
        """Return chirality for a zero-based atom index."""

        if atom_id < 0 or atom_id >= len(self.chiralities):
            raise IndexError(f"atom index out of range: {atom_id}")
        return self.chiralities[atom_id]
