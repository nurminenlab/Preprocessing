"""
Batch driver: re-run the Kilosort -> NWB conversion for every .nwb file
already sitting in NWB_DIR, so that each file gets the fixes shipped in the
current version of convKilosort2NWB_laser.py (notably the `expt_info_json`
scratch entry, which was silently dropped by the old version).

Assumptions
-----------
- Each old .nwb file name has the shape:
      {YYYY-MM-DD}-{MONKEY_NAME}-{PROJECT}.nwb
  where PROJECT is fixed (see PROJECT below) and MONKEY_NAME may contain
  hyphens (e.g. "MM005-Desta").
- The raw source data (Kilosort output, SpikeGLX .ap/.nidq .bin/.meta, and
  the session .mat file) still live under DATA_ROOT following the layout
  hard-coded in convKilosort2NWB_laser.py.
- SESSION_NAME and GATE are the same for every session (defaults 'BSD' /
  'g0'); adjust below if that ever stops being true.
- MONKEY_INFO is a per-monkey lookup for birthday (BD) and sex (SEX). The
  starter dict below only covers MM005-Desta; add every other monkey you
  have data for before running.

The rewritten .nwb files are written back into NWB_DIR, overwriting the
originals. Consider backing that folder up first.
"""

from __future__ import annotations

from pathlib import Path
import traceback

from convKilosort2NWB_laser import convert_session


# ============== CONFIGURATION ==============
NWB_DIR = Path(r"F:/ProjectData/BSD-variability/NWBs/Moderate+")
DATA_ROOT = "F:/localDATA/Electrophysiology"

PROJECT = "Inhibition-neural-information-processing"
SESSION_NAME = "BSD"
GATE = "g0"
SPECIES = "Callithrix jacchus"

# Per-monkey static info: birthday (YYYY-MM-DD) and sex. Fill in every
# monkey that appears in NWB_DIR before running.
MONKEY_INFO: dict[str, dict[str, str]] = {
    "MM005-Desta": {"BD": "2021-04-23", "SEX": "F"},
    # "MMxxx-Name": {"BD": "YYYY-MM-DD", "SEX": "M"},
}

# If True, overwrite existing .nwb files in NWB_DIR. If False, sessions
# whose output already exists will be skipped.
OVERWRITE = True
# ============== END CONFIGURATION ==============


def parse_nwb_filename(nwb_path: Path, project: str) -> tuple[str, str]:
    """
    Parse '{YYYY-MM-DD}-{MONKEY_NAME}-{PROJECT}.nwb' into (date_str, monkey_name).

    Raises ValueError if the filename does not match the expected pattern.
    """
    stem = nwb_path.stem  # strip .nwb
    project_suffix = "-" + project
    if not stem.endswith(project_suffix):
        raise ValueError(
            f"Filename does not end with '-{project}': {nwb_path.name}"
        )
    stem_no_project = stem[: -len(project_suffix)]

    # Date is the first 10 chars: YYYY-MM-DD
    if len(stem_no_project) < 12 or stem_no_project[10] != "-":
        raise ValueError(
            f"Cannot find YYYY-MM-DD prefix in filename: {nwb_path.name}"
        )
    date_str = stem_no_project[:10]
    monkey_name = stem_no_project[11:]

    if not date_str or not monkey_name:
        raise ValueError(
            f"Empty date or monkey_name parsed from filename: {nwb_path.name}"
        )
    return date_str, monkey_name


def main() -> None:
    if not NWB_DIR.is_dir():
        raise FileNotFoundError(f"NWB_DIR does not exist: {NWB_DIR}")

    nwb_files = sorted(NWB_DIR.glob("*.nwb"))
    if not nwb_files:
        print(f"No .nwb files found in {NWB_DIR}")
        return

    print(f"Found {len(nwb_files)} .nwb file(s) in {NWB_DIR}\n")

    succeeded: list[Path] = []
    skipped: list[tuple[Path, str]] = []
    failed: list[tuple[Path, str]] = []

    for nwb_path in nwb_files:
        print("=" * 72)
        print(f"Session file: {nwb_path.name}")

        try:
            date_str, monkey_name = parse_nwb_filename(nwb_path, PROJECT)
        except ValueError as exc:
            reason = f"filename parse error: {exc}"
            print(f"  SKIP -- {reason}")
            skipped.append((nwb_path, reason))
            continue

        info = MONKEY_INFO.get(monkey_name)
        if info is None:
            reason = (
                f"no MONKEY_INFO entry for {monkey_name!r}; "
                "add it to MONKEY_INFO and rerun"
            )
            print(f"  SKIP -- {reason}")
            skipped.append((nwb_path, reason))
            continue

        print(
            f"  monkey={monkey_name}  date={date_str}  "
            f"BD={info['BD']}  SEX={info['SEX']}"
        )

        try:
            out_path = convert_session(
                monkey_name=monkey_name,
                date_str=date_str,
                birthday_str=info["BD"],
                sex=info["SEX"],
                session_name=SESSION_NAME,
                gate=GATE,
                project=PROJECT,
                data_root=DATA_ROOT,
                output_dir=str(NWB_DIR),
                species=SPECIES,
                overwrite=OVERWRITE,
            )
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            print(f"  FAIL -- {reason}")
            traceback.print_exc()
            failed.append((nwb_path, reason))
            continue

        print(f"  OK   -- wrote {out_path}")
        succeeded.append(out_path)

    # ---- Summary ----
    print("\n" + "=" * 72)
    print("BATCH SUMMARY")
    print(f"  succeeded : {len(succeeded)}")
    print(f"  skipped   : {len(skipped)}")
    print(f"  failed    : {len(failed)}")

    if skipped:
        print("\nSkipped:")
        for p, r in skipped:
            print(f"  - {p.name}: {r}")
    if failed:
        print("\nFailed:")
        for p, r in failed:
            print(f"  - {p.name}: {r}")


if __name__ == "__main__":
    main()
