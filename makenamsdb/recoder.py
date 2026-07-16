"""Python 3 / Open Babel 3 reference encoder for NAMS databases.

The encoder preserves the original NAMS text database layout and its core ABA
bond descriptors while fixing Python 2 incompatibilities, the historical
fragment-test bug, nondeterministic dictionary traversal, recursive graph
walking, and unsafe binary string handling.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import struct
from typing import BinaryIO, Iterable, Literal, TextIO

try:
    from openbabel import openbabel, pybel
except ImportError as exc:  # pragma: no cover
    openbabel = None  # type: ignore[assignment]
    pybel = None  # type: ignore[assignment]
    _OPENBABEL_IMPORT_ERROR: ImportError | None = exc
else:
    _OPENBABEL_IMPORT_ERROR = None

from chirality import Chirality
from doubleb_e_z import Stereodoubleb


FragmentPolicy = Literal["legacy-first", "largest"]
AbaBond = tuple[int, int, int, int, int, int, int, int, float, int]
MolInfo = dict[int, list[list[AbaBond]]]


class NamsEncodingError(ValueError):
    """Raised when a molecule cannot be encoded into the NAMS format."""


def require_openbabel() -> None:
    if _OPENBABEL_IMPORT_ERROR is not None:
        raise RuntimeError(
            "Open Babel Python bindings are required. Install Open Babel 3 and "
            "verify: 'from openbabel import openbabel, pybel'."
        ) from _OPENBABEL_IMPORT_ERROR


def openbabel_version() -> str:
    require_openbabel()
    try:
        return str(openbabel.OBReleaseVersion())
    except AttributeError:
        return "unknown"


def _normalise_input_format(input_type: str) -> str:
    return "smi" if input_type.lower() == "can" else input_type.lower()


def _canonical_smiles(molecule: object) -> str:
    value = molecule.write("can").strip()
    if not value:
        raise NamsEncodingError("Open Babel produced an empty canonical SMILES")
    # Some writers append a title separated by whitespace.
    return value.split()[0]


def _fragment_heavy_atom_count(smiles: str) -> int:
    molecule = pybel.readstring("smi", smiles)
    return sum(1 for atom in molecule.atoms if atom.atomicnum > 1)


def select_fragment(canonical_smiles: str, policy: FragmentPolicy) -> str:
    """Select a connected component using an explicit, documented policy."""

    fragments = canonical_smiles.split(".")
    if len(fragments) == 1:
        return canonical_smiles
    if policy == "legacy-first":
        return fragments[0]
    if policy == "largest":
        # Stable tie-breaker: heavy atoms, total atoms, then canonical text.
        return max(
            fragments,
            key=lambda fragment: (
                _fragment_heavy_atom_count(fragment),
                len(pybel.readstring("smi", fragment).atoms),
                fragment,
            ),
        )
    raise NamsEncodingError(f"unsupported fragment policy: {policy}")


@dataclass(frozen=True, slots=True)
class EncodedMatrices:
    aba_types: list[AbaBond]
    atom_aba_indices: list[list[int]]
    atom_levels: list[list[int]]
    atom_count: int
    bond_count: int


class Recoder:
    """Convert molecules to the original NAMS ABA/level representation."""

    def __init__(self, fragment_policy: FragmentPolicy = "legacy-first") -> None:
        require_openbabel()
        self.fragment_policy = fragment_policy
        self.last_canonical_smiles: str | None = None

    @staticmethod
    def get_num_rings(mol: object, atom: object) -> int:
        """Count SSSR rings containing ``atom``, matching original NAMS."""

        count = 0
        for ring in mol.GetSSSR():
            if ring.IsMember(atom):
                count += 1
        return count

    @staticmethod
    def _adjacency(mol: object) -> dict[int, list[int]]:
        adjacency: dict[int, list[int]] = {
            int(atom.GetIdx()): [] for atom in openbabel.OBMolAtomIter(mol)
        }
        for bond in openbabel.OBMolBondIter(mol):
            begin = int(bond.GetBeginAtomIdx())
            end = int(bond.GetEndAtomIdx())
            adjacency[begin].append(end)
            adjacency[end].append(begin)
        return adjacency

    @staticmethod
    def _bond_levels(
        root: int,
        adjacency: dict[int, list[int]],
    ) -> list[list[tuple[int, int]]]:
        """Return each molecular bond exactly once, layered from ``root``."""

        frontier: list[int] = [root]
        visited_edges: set[tuple[int, int]] = set()
        levels: list[list[tuple[int, int]]] = []

        while frontier:
            current: list[tuple[int, int]] = []
            next_frontier: list[int] = []
            next_seen: set[int] = set()

            for at1 in frontier:
                for at2 in adjacency.get(at1, []):
                    edge = (at1, at2) if at1 < at2 else (at2, at1)
                    if edge in visited_edges:
                        continue
                    visited_edges.add(edge)
                    current.append((at1, at2))
                    if at2 not in next_seen:
                        next_seen.add(at2)
                        next_frontier.append(at2)

            if not current:
                break
            levels.append(current)
            frontier = next_frontier

        return levels

    @staticmethod
    def _get_bond(atom1: object, atom2: object) -> object:
        bond = atom1.GetBond(atom2)
        if bond is None:
            raise NamsEncodingError("internal graph traversal produced a missing bond")
        return bond

    def canonicalize(self, input_type: str, molecular_text: str) -> tuple[object, str]:
        try:
            molecule = pybel.readstring(
                _normalise_input_format(input_type), molecular_text
            )
        except (OSError, IOError, RuntimeError, ValueError) as exc:
            raise NamsEncodingError(f"cannot parse molecule: {exc}") from exc

        canonical = select_fragment(
            _canonical_smiles(molecule), self.fragment_policy
        )
        try:
            selected = pybel.readstring("smi", canonical)
        except (OSError, IOError, RuntimeError, ValueError) as exc:
            raise NamsEncodingError(
                f"cannot parse selected canonical fragment {canonical!r}: {exc}"
            ) from exc
        self.last_canonical_smiles = canonical
        return selected, canonical

    def get_mol_info(
        self,
        typ: str,
        mol_str: str,
        hydros: bool = False,
        DoIsomerism: bool = False,
    ) -> tuple[object | bool, MolInfo | bool]:
        """Return ``(pybel_molecule, mol_info)`` using the legacy public API.

        ``mol_info`` keys remain one-based heavy-atom indices.  Added explicit
        hydrogens participate in bond environments but do not become assignment
        atoms, matching the historical implementation.
        """

        molecule, canonical = self.canonicalize(typ, mol_str)
        obmol = molecule.OBMol
        n_big_atoms = len(molecule.atoms)
        if n_big_atoms < 2:
            return False, False

        chirality = Chirality(canonical, "smi") if DoIsomerism else None
        double_stereo = Stereodoubleb(canonical, "smi") if DoIsomerism else None

        if hydros:
            obmol.AddHydrogens()
        n_total_atoms = int(obmol.NumAtoms())

        atom_features: dict[int, tuple[object, int, int, int]] = {}
        for atom in openbabel.OBMolAtomIter(obmol):
            atom_idx = int(atom.GetIdx())
            chirality_value = (
                chirality.get_chirality(atom_idx - 1)
                if chirality is not None and atom_idx <= n_big_atoms
                else 0
            )
            atom_features[atom_idx] = (
                atom,
                int(atom.GetAtomicNum()),
                self.get_num_rings(obmol, atom),
                int(chirality_value),
            )

        adjacency = self._adjacency(obmol)
        mol_info: MolInfo = {}

        for root in range(1, n_big_atoms + 1):
            encoded_levels: list[list[AbaBond]] = []
            for level in self._bond_levels(root, adjacency):
                encoded_level: list[AbaBond] = []
                for at1_idx, at2_idx in level:
                    at1, num1, rings1, chir1 = atom_features[at1_idx]
                    at2, num2, rings2, chir2 = atom_features[at2_idx]
                    bond = self._get_bond(at1, at2)
                    db_stereo = (
                        int(double_stereo.get_e_z_at(at1, at2))
                        if double_stereo is not None
                        else 0
                    )
                    aromatic = bool(bond.IsAromatic())
                    order = 1.5 if aromatic else float(bond.GetBondOrder())
                    encoded_level.append(
                        (
                            num1,
                            rings1,
                            chir1,
                            num2,
                            rings2,
                            chir2,
                            int(bool(bond.IsInRing())),
                            int(aromatic),
                            order,
                            db_stereo,
                        )
                    )
                encoded_levels.append(encoded_level)
            mol_info[root] = encoded_levels

        return molecule, mol_info

    @staticmethod
    def calc_btypes(mol_info: MolInfo) -> tuple[dict[AbaBond, int], dict[int, list[tuple[int, int]]]]:
        bond_types: dict[AbaBond, int] = {}
        atom_bonds: dict[int, list[tuple[int, int]]] = {}
        for atom_id in sorted(mol_info):
            atom_bonds[atom_id] = []
            for level_index, level in enumerate(mol_info[atom_id]):
                for bond in level:
                    if bond not in bond_types:
                        bond_types[bond] = len(bond_types)
                    atom_bonds[atom_id].append((bond_types[bond], level_index))
        return bond_types, atom_bonds

    @staticmethod
    def _matrices(mol_info: MolInfo) -> EncodedMatrices:
        aba_to_index: dict[AbaBond, int] = {}
        atom_aba_indices: list[list[int]] = []
        atom_levels: list[list[int]] = []
        expected_bonds: int | None = None

        for atom_id in sorted(mol_info):
            aba_row: list[int] = []
            level_row: list[int] = []
            for level_index, level in enumerate(mol_info[atom_id]):
                if level_index >= 128:
                    raise NamsEncodingError(
                        "molecule requires a bond level >= 128, exceeding NAMS MAX_LEVELS"
                    )
                for aba in level:
                    aba_index = aba_to_index.setdefault(aba, len(aba_to_index))
                    aba_row.append(aba_index)
                    level_row.append(level_index)
            if expected_bonds is None:
                expected_bonds = len(aba_row)
            elif len(aba_row) != expected_bonds:
                raise NamsEncodingError(
                    "inconsistent bond counts between atom environments"
                )
            atom_aba_indices.append(aba_row)
            atom_levels.append(level_row)

        ordered_aba = [None] * len(aba_to_index)
        for aba, index in aba_to_index.items():
            ordered_aba[index] = aba

        return EncodedMatrices(
            aba_types=ordered_aba,  # type: ignore[arg-type]
            atom_aba_indices=atom_aba_indices,
            atom_levels=atom_levels,
            atom_count=len(atom_aba_indices),
            bond_count=expected_bonds or 0,
        )

    @staticmethod
    def _aba_as_ints(aba: AbaBond) -> tuple[int, ...]:
        return (
            int(aba[0]),
            int(aba[1]),
            int(aba[2]),
            int(aba[3]),
            int(aba[4]),
            int(aba[5]),
            int(aba[6]),
            int(aba[7]),
            int(round(aba[8] * 10.0)),
            int(aba[9]),
        )

    def export_mol_info(
        self,
        mol_info: MolInfo,
        cid: int,
        cname: str,
        fil: TextIO,
        molwt: float = -1.0,
    ) -> None:
        matrices = self._matrices(mol_info)
        fil.write(f"{int(cid)} {int(molwt * 10)} {cname}\n")
        fil.write(
            f"{matrices.atom_count} {matrices.bond_count} {len(matrices.aba_types)}\n"
        )
        for aba in matrices.aba_types:
            fil.write(" ".join(str(value) for value in self._aba_as_ints(aba)) + "\n")
        for row in matrices.atom_aba_indices:
            fil.write(" ".join(map(str, row)) + "\n")
        for row in matrices.atom_levels:
            fil.write(" ".join(map(str, row)) + "\n")

    def export_mol_info_bin(
        self,
        mol_info: MolInfo,
        cid: int,
        cname: str,
        fil: BinaryIO,
        molwt: float = -1.0,
    ) -> None:
        """Write the legacy little-endian compact binary representation.

        The original binary reader stores each ABA field in one unsigned byte.
        Consequently, signed stereo values are not representable. Use text
        output for stereo-enabled reference databases until the NAMS2 binary
        format is versioned in the C++ reader.
        """

        matrices = self._matrices(mol_info)
        if matrices.atom_count > 32767 or matrices.bond_count > 32767:
            raise NamsEncodingError("legacy binary counts exceed signed 16-bit range")
        if len(matrices.aba_types) > 255:
            raise NamsEncodingError("legacy binary ABA type count exceeds one byte")
        if any(any(value > 255 for value in row) for row in matrices.atom_aba_indices):
            raise NamsEncodingError("legacy binary ABA index exceeds one byte")

        encoded_name = cname.encode("utf-8")[:31].ljust(32, b"\0")
        fil.write(struct.pack("<II32s", int(cid), int(molwt * 10), encoded_name))
        fil.write(
            struct.pack(
                "<hhh",
                matrices.atom_count,
                matrices.bond_count,
                len(matrices.aba_types),
            )
        )
        for aba in matrices.aba_types:
            values = self._aba_as_ints(aba)
            if any(value < 0 or value > 255 for value in values):
                raise NamsEncodingError(
                    "legacy binary format cannot encode negative stereo values; "
                    "use text output"
                )
            fil.write(struct.pack("<10B", *values))
        for row in matrices.atom_aba_indices:
            fil.write(bytes(row))
        for row in matrices.atom_levels:
            if any(level < 0 or level > 255 for level in row):
                raise NamsEncodingError("legacy binary level exceeds one byte")
            fil.write(bytes(row))
