from pathlib import Path
import os

# Find repo root
def find_repo_root(path=None):
    current_dir = path or Path.cwd()
    
    for path in [current_dir, *current_dir.parents]:
        if (path / ".git").exists():
            return path
        
    raise RuntimeError(f"{current_dir} is not inside a git repository")

def find_apref_dir(path):
    
    # Find APREF folder
    apref_dirs = list(path.glob("APREF*"))

    # If not exactly one folder, raise error
    if len(apref_dirs) == 0:
        raise FileNotFoundError(f"No APREF folder found in {NGCA_DIR}")

    if len(apref_dirs) > 1:
        raise RuntimeError(f"Expected one APREF directory, found {len(apref_dirs)}: {apref_dirs}")

    return apref_dirs[0]

def get_dynadjust_env():
    """Return an environment containing the DynAdjust paths."""
    env = os.environ.copy()

    dynadjust_dir = str(DYNADJUST_DIR)

    path_dirs = env.get("PATH", "").split(os.pathsep)
    if dynadjust_dir not in path_dirs:
        env["PATH"] = dynadjust_dir + os.pathsep + env.get("PATH", "")

    library_dirs = env.get("LD_LIBRARY_PATH", "").split(os.pathsep)
    if dynadjust_dir not in library_dirs:
        env["LD_LIBRARY_PATH"] = (
            dynadjust_dir + os.pathsep + env.get("LD_LIBRARY_PATH", "")
        )

    return env

REPO_ROOT = find_repo_root()

# ---- Directories ------------------------------------------------------------

NGCA_DIR = REPO_ROOT / "ngca"
NGCA_APREF_DIR = find_apref_dir(NGCA_DIR)
NGCA_TRANS_TABLES_DIR = NGCA_DIR / "transTables"
NGCA_SENT_DIR = NGCA_DIR / "sent"

DYNADJUST_DIR = REPO_ROOT / "dynadjust-linux-static"

# ---- Env --------------------------------------------------------------------

DYNADJUST_ENV = get_dynadjust_env()