# This script will run a subprocess method on the XFOIL software to validate results from the VPM and TAT scripts
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# Calling this method is expected to yield the most accurate results since it combines VPM with integral boundary layer equations

# xfoil_wrapper.py
# Drives XFOIL 6.99 via subprocess and parses its output.
# XFOIL is NOT bundled with this repo (GPL, platform-specific).
# Install from http://web.mit.edu/drela/Public/web/xfoil/ and either put
# xfoil.exe on PATH or set the XFOIL_PATH environment variable.

import os;
import shutil;
import subprocess;
import tempfile;
import numpy as np;
import pandas as pd;
from post_processing import plot_coeffs, export_VPM_pressure;

def run_xfoil_solver(dat_path, angle_param, n_panels, viscous=False):
    xf = run_xfoil_polar(dat_path, angle_param, n_panels, viscous=viscous);

    if viscous:
        alpha_xfVisc = xf["alpha"].to_numpy();
        cl_xfVisc    = xf["CL"].to_numpy();
        cmqc_xfVisc    = xf["CM"].to_numpy();
        cmLE_xfVisc = cmqc_xfVisc - 0.25*cl_xfVisc;
        
        return cl_xfVisc, cmLE_xfVisc, cmqc_xfVisc, alpha_xfVisc;   
    else:      
        alpha_xf = xf["alpha"].to_numpy();
        cl_xf    = xf["CL"].to_numpy();
        cmqc_xf    = xf["CM"].to_numpy();
        cmLE_xf = cmqc_xf - 0.25*cl_xf;

        return cl_xf, cmLE_xf, cmqc_xf, alpha_xf;

def find_xfoil():
    explicit = os.environ.get("XFOIL_PATH");
    if explicit and os.path.isfile(explicit):
        return explicit;

    found = shutil.which("xfoil") or shutil.which("xfoil.exe");
    if found:
        return found;

    # Common local install locations, checked last
    for candidate in [
        os.path.expanduser(r"~\XFOIL6.99\xfoil.exe"),
    ]:
        if os.path.isfile(candidate):
            return candidate;

    raise FileNotFoundError(...);


def run_xfoil_polar(dat_path, alphas, n_panels=None, viscous=False, Re=1e6,
                    timeout=60):
    xfoil = find_xfoil();
    dat_path = os.path.abspath(dat_path);

    with tempfile.TemporaryDirectory() as workdir:
        # XFOIL 6.99 truncates paths containing ':' (Windows drive letters),
        # so copy the coordinate file in and reference it by bare filename.
        local_dat = "airfoil.dat";
        shutil.copy(dat_path, os.path.join(workdir, local_dat));

        polar_name = "polar.txt";                       # bare name for XFOIL
        polar_file = os.path.join(workdir, polar_name); # full path for Python

        cmds = [];
        cmds.append(f"LOAD {local_dat}");     # <-- bare filename, no drive letter
        cmds.append("");

        if n_panels is not None:
            cmds.append("PPAR");
            cmds.append(f"N {n_panels}");
            cmds.append("");
            cmds.append("");

        cmds.append("ITER 100");

        cmds.append("OPER");
        if viscous:
            cmds.append(f"VISC {Re}");
        cmds.append("PACC");
        cmds.append(polar_name);              # <-- bare filename here too
        cmds.append("");

        for a in alphas:
            cmds.append(f"ALFA {a:.4f}");

        cmds.append("PACC");
        cmds.append("");
        cmds.append("QUIT");

        stdin = "\n".join(cmds) + "\n";

        try:
            proc = subprocess.run(
                [xfoil], input=stdin, cwd=workdir,
                capture_output=True, text=True, timeout=timeout
            );
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"XFOIL timed out after {timeout}s. This usually means a "
                f"command was malformed and XFOIL is waiting at a prompt."
            );

        if not os.path.isfile(polar_file):
            raise RuntimeError(
                "XFOIL produced no polar file. Last output:\n"
                + proc.stdout[-2000:]
            );

        return parse_polar(polar_file);


def parse_polar(polar_file):
    """Parse an XFOIL PACC polar file into a DataFrame."""
    rows = [];
    with open(polar_file, "r") as f:
        lines = f.readlines();

    # The numeric block starts after a line of dashes.
    start = None;
    for i, line in enumerate(lines):
        if line.strip().startswith("---"):
            start = i + 1;
            break;

    if start is None:
        raise RuntimeError(f"Could not find data block in {polar_file}");

    for line in lines[start:]:
        tokens = line.split();
        if len(tokens) < 5:
            continue;
        try:
            rows.append([float(t) for t in tokens[:5]]);
        except ValueError:
            continue;

    if not rows:
        raise RuntimeError(f"No numeric rows parsed from {polar_file}");

    arr = np.array(rows);
    # XFOIL polar columns: alpha  CL  CD  CDp  CM
    return pd.DataFrame({
        "alpha": arr[:, 0],
        "CL":    arr[:, 1],
        "CD":    arr[:, 2],
        "CM":    arr[:, 4],
    });


def run_xfoil_cp(dat_path, alpha, n_panels=None, timeout=60):
    xfoil = find_xfoil();
    dat_path = os.path.abspath(dat_path);

    with tempfile.TemporaryDirectory() as workdir:
        local_dat = "airfoil.dat";
        shutil.copy(dat_path, os.path.join(workdir, local_dat));

        cp_name = "cp.txt";                            # bare name for XFOIL
        cp_file = os.path.join(workdir, cp_name);      # full path for Python

        cmds = [f"LOAD {local_dat}", ""];
        if n_panels is not None:
            cmds += ["PPAR", f"N {n_panels}", "", ""];
        cmds += ["OPER", f"ALFA {alpha:.4f}", f"CPWR {cp_name}", "", "QUIT"];

        proc = subprocess.run([xfoil], input="\n".join(cmds) + "\n", cwd=workdir,
                              capture_output=True, text=True, timeout=timeout);

        if not os.path.isfile(cp_file):
            raise RuntimeError(
                "XFOIL produced no Cp file.\n"
                f"--- stdout ---\n{proc.stdout}\n"
                f"--- stderr ---\n{proc.stderr}\n"
                f"--- commands sent ---\n" + "\n".join(cmds)
                );

        # XFOIL's CPWR header varies by version (often two lines, one containing
        # "Alfa"). Skip anything non-numeric rather than assuming a fixed count.
        rows = [];
        with open(cp_file, "r") as f:
            for line in f:
                tokens = line.split();
                if len(tokens) < 2:
                    continue;
                try:
                    rows.append([float(t) for t in tokens]);
                except ValueError:
                    continue;      # header line

        if not rows:
            raise RuntimeError(f"No numeric rows parsed from {cp_file}");

        data = np.array(rows);

        if data.shape[1] >= 3:
            return pd.DataFrame({"x": data[:,0], "y": data[:,1], "Cp": data[:,2]});
        else:
            return pd.DataFrame({"x": data[:,0], "Cp": data[:,1]});