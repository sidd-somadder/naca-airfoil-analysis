"""
Drives XFOIL 6.99 via subprocess and parses its output, providing an independent
reference for validating the VPM and TAT solvers.

XFOIL is NOT bundled with this repo (GPL, platform-specific). Install from
http://web.mit.edu/drela/Public/web/xfoil/ and either place the executable on
PATH or set the XFOIL_PATH environment variable.
"""

import os
import shutil
import subprocess
import tempfile
import numpy as np
import pandas as pd

def run_xfoil_solver(dat_path, angle_param, n_panels, viscous=False):
    """
    Runs XFOIL over an angle of attack sweep and returns lift and moment coefficients.

    Inputs:
    Absolute path to the airfoil .dat file (dat_path).
    Angle of Attack sweep (angle_param) in DEGREES (size N).
    Panel count (n_panels) passed to XFOIL's PPAR repaneling.
    Viscous mode flag (viscous); when True, runs at the Reynolds number set in run_xfoil_polar.

    Outputs:
    Coefficient of Lift (cl), Coefficient of Moment about the Leading Edge (cmLE) and
    about the Quarter Chord (cmqc), and the angle array XFOIL actually returned.

    Note: XFOIL reports CM about the quarter chord; cmLE is derived as cmqc - 0.25*cl.

    The returned angle array is provided separately because viscous runs omit angles where
    the boundary layer solution fails to converge, giving fewer points than requested.
    """

    xf = run_xfoil_polar(dat_path, angle_param, n_panels, viscous=viscous)

    if viscous:
        alpha_xfVisc = xf["alpha"].to_numpy()
        cl_xfVisc    = xf["CL"].to_numpy()
        cmqc_xfVisc    = xf["CM"].to_numpy()
        cmLE_xfVisc = cmqc_xfVisc - 0.25*cl_xfVisc
        
        return cl_xfVisc, cmLE_xfVisc, cmqc_xfVisc, alpha_xfVisc   
    else:     
        alpha_xf = xf["alpha"].to_numpy()
        cl_xf    = xf["CL"].to_numpy()
        cmqc_xf    = xf["CM"].to_numpy()
        cmLE_xf = cmqc_xf - 0.25*cl_xf

        return cl_xf, cmLE_xf, cmqc_xf, alpha_xf

def find_xfoil():
    """
    Locates the XFOIL executable on the host system.

    Inputs:
    None. Checks the XFOIL_PATH environment variable, then PATH, then common install locations.

    Outputs:
    Absolute path to the XFOIL executable.

    Raises FileNotFoundError with install instructions if XFOIL cannot be located.
    """

    explicit = os.environ.get("XFOIL_PATH")
    if explicit and os.path.isfile(explicit):
        return explicit

    found = shutil.which("xfoil") or shutil.which("xfoil.exe")
    if found:
        return found
    # Common local install locations, checked last
    for candidate in [
        os.path.expanduser(r"~\XFOIL6.99\xfoil.exe"),
    ]:
        if os.path.isfile(candidate):
            return candidate
        
    raise FileNotFoundError(
        "XFOIL executable not found. Install from "
        "http://web.mit.edu/drela/Public/web/xfoil/ and either place it on PATH "
        "or set the XFOIL_PATH environment variable."
    )   

def run_xfoil_polar(dat_path, alphas, n_panels=None, viscous=False, Re=1e6, timeout=60):
    """
    Drives XFOIL via subprocess to produce a polar over an angle of attack sweep.

    Inputs:
    Absolute path to the airfoil .dat file (dat_path).
    Angle of Attack sweep (alphas) in DEGREES (size N).
    Optional panel count (n_panels) for XFOIL's PPAR repaneling; XFOIL's default is used if None.
    Viscous mode flag (viscous) and Reynolds number (Re), which has no effect in inviscid mode.
    Subprocess timeout in seconds (timeout).

    Outputs:
    DataFrame with columns alpha, CL, CD, and CM, containing only the angles XFOIL converged on.

    Note: the coordinate file is copied into a temporary working directory and referenced by bare
    filename, because XFOIL 6.99 truncates paths containing a Windows drive letter colon.
    """
    
    xfoil = find_xfoil()
    dat_path = os.path.abspath(dat_path)

    with tempfile.TemporaryDirectory() as workdir:
        # XFOIL 6.99 truncates paths containing ':' (Windows drive letters),
        # so copy the coordinate file in and reference it by bare filename.
        local_dat = "airfoil.dat"
        shutil.copy(dat_path, os.path.join(workdir, local_dat))

        # bare name for XFOIL
        polar_name = "polar.txt"
        # full path for Python
        polar_file = os.path.join(workdir, polar_name) 

        cmds = []
        # bare filename, no drive letter
        cmds.append(f"LOAD {local_dat}")     
        cmds.append("")

        if n_panels is not None:
            cmds.append("PPAR")
            cmds.append(f"N {n_panels}")
            cmds.append("")
            cmds.append("")

        cmds.append("ITER 100")

        cmds.append("OPER")
        if viscous:
            cmds.append(f"VISC {Re}")
        cmds.append("PACC")
        # bare filename here too
        cmds.append(polar_name)             
        cmds.append("")

        for a in alphas:
            cmds.append(f"ALFA {a:.4f}")

        cmds.append("PACC")
        cmds.append("")
        cmds.append("QUIT")

        stdin = "\n".join(cmds) + "\n"

        try:
            proc = subprocess.run(
                [xfoil], input=stdin, cwd=workdir,
                capture_output=True, text=True, timeout=timeout
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"XFOIL timed out after {timeout}s. This usually means a "
                f"command was malformed and XFOIL is waiting at a prompt."
            )

        if not os.path.isfile(polar_file):
            raise RuntimeError(
                "XFOIL produced no polar file. Last output:\n"
                + proc.stdout[-2000:]
            )

        return parse_polar(polar_file)

def parse_polar(polar_file):
    """Parse an XFOIL PACC polar file into a DataFrame."""
    rows = []
    with open(polar_file, "r") as f:
        lines = f.readlines()

    # The numeric block starts after a line of dashes.
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("---"):
            start = i + 1
            break

    if start is None:
        raise RuntimeError(f"Could not find data block in {polar_file}")

    for line in lines[start:]:
        tokens = line.split()
        if len(tokens) < 5:
            continue
        try:
            rows.append([float(t) for t in tokens[:5]])
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric rows parsed from {polar_file}")

    arr = np.array(rows)
    # XFOIL polar columns: alpha  CL  CD  CDp  CM
    return pd.DataFrame({
        "alpha": arr[:, 0],
        "CL":    arr[:, 1],
        "CD":    arr[:, 2],
        "CM":    arr[:, 4],
    })

def run_xfoil_cp(dat_path, alpha, n_panels=None, timeout=60):
    """
    Drives XFOIL via subprocess to produce a pressure distribution at a single angle of attack.

    Inputs:
    Absolute path to the airfoil .dat file (dat_path).
    Angle of Attack (alpha) in DEGREES.
    Optional panel count (n_panels) for XFOIL's PPAR repaneling.
    Subprocess timeout in seconds (timeout).

    Outputs:
    DataFrame with columns x, y, and Cp in standard selig order (TE -> upper -> LE -> lower -> TE).
    Older XFOIL versions omit the y column, in which case only x and Cp are returned.
    Runs inviscid only.
    """

    xfoil = find_xfoil()
    dat_path = os.path.abspath(dat_path)

    with tempfile.TemporaryDirectory() as workdir:
        local_dat = "airfoil.dat"
        shutil.copy(dat_path, os.path.join(workdir, local_dat))

        cp_name = "cp.txt"                         # bare name for XFOIL
        cp_file = os.path.join(workdir, cp_name)      # full path for Python

        cmds = [f"LOAD {local_dat}", ""]
        if n_panels is not None:
            cmds += ["PPAR", f"N {n_panels}", "", ""]
        cmds += ["OPER", f"ALFA {alpha:.4f}", f"CPWR {cp_name}", "", "QUIT"]

        proc = subprocess.run([xfoil], input="\n".join(cmds) + "\n", cwd=workdir,
                              capture_output=True, text=True, timeout=timeout)

        if not os.path.isfile(cp_file):
            raise RuntimeError(
                "XFOIL produced no Cp file.\n"
                f"--- stdout ---\n{proc.stdout}\n"
                f"--- stderr ---\n{proc.stderr}\n"
                f"--- commands sent ---\n" + "\n".join(cmds)
                )

        # XFOIL's CPWR header varies by version (often two lines, one containing
        # "Alfa"). Skip anything non-numeric rather than assuming a fixed count.
        rows = []
        with open(cp_file, "r") as f:
            for line in f:
                tokens = line.split()
                if len(tokens) < 2:
                    continue
                try:
                    rows.append([float(t) for t in tokens])
                except ValueError:
                    continue      # header line

        if not rows:
            raise RuntimeError(f"No numeric rows parsed from {cp_file}")

        data = np.array(rows)

        if data.shape[1] >= 3:
            return pd.DataFrame({"x": data[:,0], "y": data[:,1], "Cp": data[:,2]})
        else:
            return pd.DataFrame({"x": data[:,0], "Cp": data[:,1]})