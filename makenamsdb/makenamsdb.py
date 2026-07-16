"""Build a NAMS database using Python 3 and Open Babel 3/Pybel."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Iterator, TextIO

try:
    from openbabel import openbabel, pybel
except ImportError as exc:  # pragma: no cover
    openbabel = None  # type: ignore[assignment]
    pybel = None  # type: ignore[assignment]
    _OPENBABEL_IMPORT_ERROR: ImportError | None = exc
else:
    _OPENBABEL_IMPORT_ERROR = None

from recoder import NamsEncodingError, Recoder, openbabel_version


@dataclass(slots=True)
class BuildStats:
    read: int = 0
    written: int = 0
    skipped_blank: int = 0
    skipped_invalid: int = 0


def require_openbabel() -> None:
    if _OPENBABEL_IMPORT_ERROR is not None:
        raise RuntimeError(
            "Open Babel Python bindings are required. Install Open Babel 3 and "
            "verify: 'from openbabel import openbabel, pybel'."
        ) from _OPENBABEL_IMPORT_ERROR


def mol2svg(output_path: Path, molecule: object) -> None:
    """Write an indexed monochrome SVG using Open Babel."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    conversion = openbabel.OBConversion()
    if not conversion.SetOutFormat("svg"):
        raise NamsEncodingError("Open Babel SVG output plugin is unavailable")
    for option, value in (("d", ""), ("r", "1"), ("c", "1"), ("i", "1"), ("u", "1")):
        conversion.AddOption(option, conversion.OUTOPTIONS, value)
    if not conversion.WriteFile(molecule.OBMol, str(output_path)):
        raise NamsEncodingError(f"cannot write SVG {output_path}")
    conversion.CloseOutFile()


def _iter_records(
    input_file: TextIO,
    *,
    delimiter: str,
    smiles_column: int,
    id_column: int,
    skip_header: bool,
) -> Iterator[tuple[int, str | None, int | None, str | None]]:
    """Yield parsed rows without letting one malformed row stop streaming."""

    for line_number, raw_line in enumerate(input_file, start=1):
        if line_number == 1 and skip_header:
            continue
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            yield line_number, None, None, None
            continue
        columns = stripped.split(delimiter)
        try:
            molecular_text = columns[smiles_column].strip()
            cid = int(columns[id_column].strip())
        except (IndexError, ValueError):
            yield (
                line_number,
                None,
                None,
                f"expected integer ID and molecular text in columns "
                f"{id_column + 1} and {smiles_column + 1}",
            )
            continue
        if not molecular_text:
            yield line_number, None, None, "molecular text is empty"
            continue
        yield line_number, molecular_text, cid, None


def recode_file(
    fname_in: str | Path,
    fname_out: str | Path,
    *,
    input_format: str = "smi",
    binary: bool = False,
    make_svg: bool = False,
    svg_directory: str | Path | None = None,
    stereo: bool = False,
    explicit_hydrogens: bool = False,
    fragment_policy: str = "legacy-first",
    delimiter: str = "\t",
    smiles_column: int = 0,
    id_column: int = 1,
    skip_header: bool = False,
    strict: bool = False,
    quiet: bool = False,
) -> BuildStats:
    """Convert a delimited molecular input file into a NAMS database."""

    require_openbabel()
    recoder = Recoder(fragment_policy=fragment_policy)  # type: ignore[arg-type]
    input_path = Path(fname_in)
    output_path = Path(fname_out)
    svg_path = Path(svg_directory) if svg_directory else output_path.parent / "svg"
    stats = BuildStats()

    output_mode = "wb" if binary else "w"
    output_kwargs = {} if binary else {"encoding": "utf-8", "newline": "\n"}

    with input_path.open("r", encoding="utf-8-sig", newline=None) as input_file, output_path.open(
        output_mode, **output_kwargs
    ) as output_file:
        for line_number, molecular_text, cid, row_error in _iter_records(
            input_file,
            delimiter=delimiter,
            smiles_column=smiles_column,
            id_column=id_column,
            skip_header=skip_header,
        ):
            if molecular_text is None and row_error is None:
                stats.skipped_blank += 1
                continue
            if row_error is not None:
                stats.skipped_invalid += 1
                message = f"line {line_number}: {row_error}"
                if strict:
                    raise NamsEncodingError(message)
                print(f"WARNING: {message}", file=sys.stderr)
                continue
            assert molecular_text is not None and cid is not None
            stats.read += 1
            try:
                molecule, mol_info = recoder.get_mol_info(
                    input_format,
                    molecular_text,
                    hydros=explicit_hydrogens,
                    DoIsomerism=stereo,
                )
                if molecule is False or mol_info is False:
                    raise NamsEncodingError("molecule has fewer than two assignment atoms")
                canonical = recoder.last_canonical_smiles or molecule.write("can").strip()
                if make_svg:
                    mol2svg(svg_path / f"{cid}.svg", molecule)
                if binary:
                    recoder.export_mol_info_bin(
                        mol_info, cid, canonical, output_file, molecule.molwt
                    )
                else:
                    recoder.export_mol_info(
                        mol_info, cid, canonical, output_file, molecule.molwt
                    )
            except (NamsEncodingError, OSError, IOError, RuntimeError, ValueError) as exc:
                stats.skipped_invalid += 1
                message = f"line {line_number}, molecule {cid}: {exc}"
                if strict:
                    raise NamsEncodingError(message) from exc
                print(f"WARNING: {message}", file=sys.stderr)
                continue

            stats.written += 1
            if not quiet:
                print(
                    f"{cid}: {canonical}  MW={molecule.molwt:.3f}  written",
                    flush=True,
                )

    return stats


def _decode_delimiter(value: str) -> str:
    decoded = bytes(value, "utf-8").decode("unicode_escape")
    if not decoded:
        raise argparse.ArgumentTypeError("delimiter cannot be empty")
    return decoded


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a NAMS text or legacy binary database using Open Babel 3.",
    )
    parser.add_argument("input", type=Path, help="Delimited input molecular file")
    parser.add_argument("output", type=Path, help="Output NAMS database")
    parser.add_argument(
        "--input-format",
        default="smi",
        help="Open Babel input format for the molecular column (default: smi)",
    )
    parser.add_argument("--binary", action="store_true", help="Write legacy compact binary")
    parser.add_argument(
        "--stereo",
        action="store_true",
        help="Encode atom chirality and double-bond E/Z (text output recommended)",
    )
    parser.add_argument(
        "--explicit-hydrogens",
        action="store_true",
        help="Add explicit hydrogens to bond environments",
    )
    parser.add_argument(
        "--fragment-policy",
        choices=("legacy-first", "largest"),
        default="legacy-first",
        help="Disconnected-component policy (default: legacy-first)",
    )
    parser.add_argument(
        "--delimiter",
        type=_decode_delimiter,
        default="\t",
        help=r"Column delimiter; escape sequences accepted, e.g. '\t'",
    )
    parser.add_argument(
        "--smiles-column",
        type=int,
        default=1,
        help="One-based molecular-text column (default: 1)",
    )
    parser.add_argument(
        "--id-column",
        type=int,
        default=2,
        help="One-based integer-ID column (default: 2)",
    )
    parser.add_argument("--skip-header", action="store_true")
    parser.add_argument("--write-images", action="store_true", help="Write indexed SVG files")
    parser.add_argument("--svg-directory", type=Path)
    parser.add_argument("--strict", action="store_true", help="Stop on the first invalid record")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Manifest path (default: OUTPUT.manifest.json)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.smiles_column < 1 or args.id_column < 1:
        parser.error("column numbers are one-based and must be >= 1")
    if args.binary and args.stereo:
        parser.error(
            "the legacy binary format cannot represent -1 stereo values safely; "
            "use text output for --stereo"
        )

    try:
        stats = recode_file(
            args.input,
            args.output,
            input_format=args.input_format,
            binary=args.binary,
            make_svg=args.write_images,
            svg_directory=args.svg_directory,
            stereo=args.stereo,
            explicit_hydrogens=args.explicit_hydrogens,
            fragment_policy=args.fragment_policy,
            delimiter=args.delimiter,
            smiles_column=args.smiles_column - 1,
            id_column=args.id_column - 1,
            skip_header=args.skip_header,
            strict=args.strict,
            quiet=args.quiet,
        )
    except (NamsEncodingError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    manifest_path = args.manifest or Path(str(args.output) + ".manifest.json")
    manifest = {
        "nams_database_version": "1-compatible-reference",
        "python_reference_version": "2.0.0-ob3-ref1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "chemistry_backend": "openbabel-pybel",
        "openbabel_version": openbabel_version(),
        "python_version": sys.version.split()[0],
        "input": str(args.input.resolve()),
        "output": str(args.output.resolve()),
        "output_format": "legacy-binary" if args.binary else "text",
        "input_format": args.input_format,
        "stereo": bool(args.stereo),
        "explicit_hydrogens": bool(args.explicit_hydrogens),
        "fragment_policy": args.fragment_policy,
        "columns": {
            "molecule": args.smiles_column,
            "id": args.id_column,
            "delimiter_repr": repr(args.delimiter),
        },
        "stats": asdict(stats),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if not args.quiet:
        print(
            f"Completed: {stats.written} written, {stats.skipped_invalid} invalid, "
            f"{stats.skipped_blank} blank/comment lines."
        )
        print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
