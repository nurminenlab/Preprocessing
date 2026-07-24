"""
Inspect an existing NWB file to check whether expt_metadata was actually
persisted into experiment_description (or anywhere else).
"""
from pathlib import Path
import h5py
from pynwb import NWBHDF5IO

# --- EDIT THIS: point to an NWB file your script has already produced ---
NWB_PATH = r"F:/ProjectData/BSD-variability/NWBs/Moderate+/2025-11-17-MM005-Desta-Inhibition-neural-information-processing.nwb"

nwb_path = Path(NWB_PATH)
assert nwb_path.exists(), f"File not found: {nwb_path}"
print(f"Inspecting: {nwb_path}\n")

# --- 1. What pynwb sees ---
with NWBHDF5IO(str(nwb_path), "r") as io:
    nwb = io.read()
    desc = nwb.experiment_description or ""
    print("=== via pynwb ===")
    print(f"experiment_description length: {len(desc)}")
    print(f"contains 'Experiment Info:'   : {'Experiment Info:' in desc}")
    print("--- first 500 chars ---")
    print(desc[:500])
    print("--- last 500 chars ---")
    print(desc[-500:] if len(desc) > 500 else "(same as above)")
    print(f"\nscratch keys                  : {list(nwb.scratch.keys()) if nwb.scratch else '(none)'}")

# --- 2. What's actually on disk (bypass pynwb) ---
print("\n=== via raw h5py ===")
with h5py.File(str(nwb_path), "r") as f:
    # In NWB 2.x, experiment_description is stored as a dataset under /general
    if "general/experiment_description" in f:
        raw = f["general/experiment_description"][()]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        print(f"/general/experiment_description length: {len(raw)}")
        print(f"contains 'Experiment Info:'           : {'Experiment Info:' in raw}")
        print("--- first 500 chars ---")
        print(raw[:500])
        print("--- last 500 chars ---")
        print(raw[-500:] if len(raw) > 500 else "(same as above)")
    else:
        print("/general/experiment_description not present on disk")

    # Also look for any scratch or unexpected place expt_metadata might have landed
    print("\nTop-level groups:", list(f.keys()))
    if "general" in f:
        print("/general contents:", list(f["general"].keys()))
    if "scratch" in f:
        print("/scratch contents:", list(f["scratch"].keys()))