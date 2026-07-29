"""
Scan existing .nwb files under NWB_DIR and report whether the SpikeGLX
gate (g0, g1, ...) can be recovered from anything stored inside them.

Approach
--------
The converter (convKilosort2NWB_laser.convert_session) does not write the
gate as a dedicated NWB field. But several string fields written by
neuroconv / the KiloSortSortingInterface tend to contain the SpikeGLX
folder or filenames, which look like:

    .../{monkey}/{date}/{session}_{gate}/{session}_{gate}_imec0/...
    {session}_{gate}_t0.imec0.ap.bin
    {session}_{gate}_t0.nidq.bin

so tokens like `_g0_`, `_g0/`, `_g1_`, `_g1/`, or a trailing `_g0` / `_g1`
often survive somewhere in the file (e.g. `source_script`, electrode-group
descriptions, `session_id`).

This script opens each .nwb file read-only with h5py (no pynwb needed),
walks the entire HDF5 tree, and searches every string dataset and every
string attribute for a regex matching `_g<N>` where N is one or more
digits. It then prints:

  - the set of distinct gate values found in that file,
  - up to a few example locations (dataset/attr path + snippet) so you can
    see where the information came from.

Nothing is modified. Safe to run against production NWBs.
"""

from __future__ import annotations

from pathlib import Path
import re

import h5py


# ---- edit these if needed ----
NWB_DIR = Path(r"F:/ProjectData/BSD-variability/NWBs/Moderate+")
MAX_EXAMPLES_PER_FILE = 5     # how many hit locations to print per file
SNIPPET_CONTEXT = 40          # chars of context around each match
# ------------------------------


# Match `_g0`, `_g12`, ... when followed by `_`, `/`, `\`, `.`, or end-of-string.
# Also allow the same pattern at the very start of a string (rare but cheap).
GATE_RE = re.compile(r"(?:^|[^A-Za-z0-9])g(\d+)(?=[_./\\]|$)")


def _iter_strings_from_value(value) -> list[str]:
    """Return a list of Python strs extracted from an h5py dataset value.

    Handles scalar bytes, scalar str, 0-d/1-d ndarrays of bytes or str.
    Non-string data returns an empty list.
    """
    out: list[str] = []

    def _add(v):
        if isinstance(v, bytes):
            try:
                out.append(v.decode("utf-8", errors="replace"))
            except Exception:
                pass
        elif isinstance(v, str):
            out.append(v)

    if isinstance(value, (bytes, str)):
        _add(value)
        return out

    # numpy arrays / h5py results
    try:
        flat = value.flatten() if hasattr(value, "flatten") else [value]
    except Exception:
        return out
    for item in flat:
        _add(item)
    return out


def _scan_string(text: str, path: str, hits: list[tuple[str, str, str]]) -> None:
    """Append (path, gate, snippet) for each gate match found in `text`."""
    for m in GATE_RE.finditer(text):
        gate = f"g{m.group(1)}"
        start = max(0, m.start() - SNIPPET_CONTEXT)
        end = min(len(text), m.end() + SNIPPET_CONTEXT)
        snippet = text[start:end].replace("\n", " ").strip()
        hits.append((path, gate, snippet))


def scan_file(nwb_path: Path) -> list[tuple[str, str, str]]:
    """Return list of (hdf5_path, gate, snippet) for every match in one file."""
    hits: list[tuple[str, str, str]] = []

    with h5py.File(nwb_path, "r") as f:
        # Root-level attributes
        for aname, aval in f.attrs.items():
            for s in _iter_strings_from_value(aval):
                _scan_string(s, f"/@{aname}", hits)

        def visit(name: str, obj) -> None:
            # Attributes on this group/dataset
            for aname, aval in obj.attrs.items():
                for s in _iter_strings_from_value(aval):
                    _scan_string(s, f"/{name}@{aname}", hits)
            # Dataset contents (only if plausibly stringy and not huge)
            if isinstance(obj, h5py.Dataset):
                dtype = obj.dtype
                # Skip datasets that are clearly numeric or too large.
                is_stringy = (
                    dtype.kind in ("O", "S", "U")
                    or (dtype.kind == "V" and dtype.names is None)
                )
                if not is_stringy:
                    return
                if obj.size > 10_000:
                    return
                try:
                    val = obj[()]
                except Exception:
                    return
                for s in _iter_strings_from_value(val):
                    _scan_string(s, f"/{name}", hits)

        f.visititems(visit)

    return hits


def _summarize_gates(hits: list[tuple[str, str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for _, gate, _ in hits:
        counts[gate] = counts.get(gate, 0) + 1
    return counts


def main() -> None:
    if not NWB_DIR.is_dir():
        raise FileNotFoundError(f"NWB_DIR does not exist: {NWB_DIR}")

    nwb_files = sorted(NWB_DIR.glob("*.nwb"))
    if not nwb_files:
        print(f"No .nwb files found in {NWB_DIR}")
        return

    print(f"Scanning {len(nwb_files)} file(s) in {NWB_DIR}\n")

    per_file: list[tuple[Path, dict[str, int], list[tuple[str, str, str]]]] = []

    for nwb_path in nwb_files:
        print("=" * 72)
        print(f"{nwb_path.name}")
        try:
            hits = scan_file(nwb_path)
        except Exception as exc:
            print(f"  ERROR opening file: {type(exc).__name__}: {exc}")
            continue

        counts = _summarize_gates(hits)
        per_file.append((nwb_path, counts, hits))

        if not counts:
            print("  gate: NOT FOUND -- no `g<N>` token seen anywhere")
        else:
            gate_summary = ", ".join(
                f"{g} ({counts[g]} hit{'s' if counts[g] != 1 else ''})"
                for g in sorted(counts)
            )
            print(f"  gate candidates: {gate_summary}")
            print("  examples:")
            for path, gate, snippet in hits[:MAX_EXAMPLES_PER_FILE]:
                print(f"    [{gate}] {path}")
                print(f"        ...{snippet}...")

    # Final compact summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    for nwb_path, counts, _ in per_file:
        if not counts:
            verdict = "unknown"
        elif len(counts) == 1:
            verdict = next(iter(counts))
        else:
            verdict = "AMBIGUOUS: " + ",".join(sorted(counts))
        print(f"  {nwb_path.name}: {verdict}")


if __name__ == "__main__":
    main()
