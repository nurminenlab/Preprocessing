"""
Batch driver: re-run the Kilosort -> NWB conversion for a hard-coded list
of sessions (see SESSIONS below), so that each output .nwb file gets the
fixes shipped in the current version of convKilosort2NWB_laser.py (notably
the `expt_info_json` scratch entry, which was silently dropped by the old
version).

Assumptions
-----------
- Output .nwb file names have the shape:
      {YYYY-MM-DD}-{MONKEY_NAME}-{PROJECT}.nwb
  where PROJECT is fixed (see PROJECT below) and MONKEY_NAME may contain
  hyphens (e.g. "MM005-Desta").
- The raw source data (Kilosort output, SpikeGLX .ap/.nidq .bin/.meta, and
  the session .mat file) still live under DATA_ROOT following the layout
  hard-coded in convKilosort2NWB_laser.py.
- SESSION_NAME is the same for every session (default 'BSD'); the gate may
  differ per session ('g0' or 'g1') and is specified per-entry in SESSIONS.
- MONKEY_INFO is a per-monkey lookup for birthday (BD) and sex (SEX). Add
  an entry for every monkey referenced in SESSIONS before running.

Output .nwb files are written into NWB_DIR. If OVERWRITE is True, existing
files with the same name are replaced. Consider backing that folder up
first.

How to edit the session list
----------------------------
Each entry in SESSIONS is a dict with three keys:
    "date"   : "YYYY-MM-DD"     session date
    "monkey" : "MMxxx-Name"     must match a key in MONKEY_INFO
    "gate"   : "g0" or "g1"     SpikeGLX gate index for this session

Add / remove / edit lines freely; blank lines and `# ...` comments between
entries are fine.
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
SPECIES = "Callithrix jacchus"

# Per-monkey static info: birthday (YYYY-MM-DD) and sex. Fill in every
# monkey referenced in SESSIONS before running.
MONKEY_INFO: dict[str, dict[str, str]] = {
    "MM005-Desta": {"BD": "2021-04-22", "SEX": "F"},
    "MM001-Sansa": {"BD": "2019-08-27", "SEX": "F"},
}

# ---- Sessions to re-analyze ----
# One dict per session. `gate` is per-session because some recordings use
# g0 and others g1. Order does not matter.
SESSIONS: list[dict[str, str]] = [
    {"date": "2024-07-03", "monkey": "MM001-Sansa", "gate": "g1"},
    {"date": "2024-10-29", "monkey": "MM001-Sansa", "gate": "g1"},
    {"date": "2024-10-31", "monkey": "MM001-Sansa", "gate": "g0"},
    {"date": "2025-11-13", "monkey": "MM005-Desta", "gate": "g1"},
    {"date": "2025-11-17", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-11-19", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-11-24", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-11-25", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-01", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-02", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-08", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-09", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-12", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-12-17", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-03-17", "monkey": "MM005-Desta", "gate": "g0"},
    {"date": "2025-03-24", "monkey": "MM005-Desta", "gate": "g0"},
]

# If True, overwrite existing .nwb files in NWB_DIR. If False, sessions
# whose output already exists will be skipped by convert_session().
OVERWRITE = True
# ============== END CONFIGURATION ==============


def _describe(session: dict[str, str]) -> str:
    """Short human-readable label for log lines."""
    return (
        f"{session.get('date', '?')}-{session.get('monkey', '?')} "
        f"[{session.get('gate', '?')}]"
    )


def main() -> None:
    NWB_DIR.mkdir(parents=True, exist_ok=True)

    if not SESSIONS:
        print("SESSIONS is empty; nothing to do.")
        return

    print(f"Processing {len(SESSIONS)} session(s); output -> {NWB_DIR}\n")

    succeeded: list[Path] = []
    skipped: list[tuple[str, str]] = []
    failed: list[tuple[str, str]] = []

    for session in SESSIONS:
        label = _describe(session)        

        date_str = session.get("date", "").strip()
        monkey_name = session.get("monkey", "").strip()
        gate = session.get("gate", "").strip()

        missing = [k for k, v in (
            ("date", date_str), ("monkey", monkey_name), ("gate", gate)
        ) if not v]
        if missing:
            reason = f"missing required field(s): {', '.join(missing)}"
            print(f"  SKIP -- {reason}")
            skipped.append((label, reason))
            continue

        info = MONKEY_INFO.get(monkey_name)
        if info is None:
            reason = (
                f"no MONKEY_INFO entry for {monkey_name!r}; "
                "add it to MONKEY_INFO and rerun"
            )
            print(f"  SKIP -- {reason}")
            skipped.append((label, reason))
            continue

        print(
            f"  monkey={monkey_name}  date={date_str}  gate={gate}  "
            f"BD={info['BD']}  SEX={info['SEX']}"
        )

        try:
            out_path = convert_session(
                monkey_name=monkey_name,
                date_str=date_str,
                birthday_str=info["BD"],
                sex=info["SEX"],
                session_name=SESSION_NAME,
                gate=gate,
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
            failed.append((label, reason))
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
        for lbl, r in skipped:
            print(f"  - {lbl}: {r}")
    if failed:
        print("\nFailed:")
        for lbl, r in failed:
            print(f"  - {lbl}: {r}")


if __name__ == "__main__":
    main()
