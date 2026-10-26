from pathlib import Path

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

REPO_ROOT = find_repo_root()

# ---- Directories ------------------------------------------------------------

NGCA_DIR = REPO_ROOT / "ngca"
NGCA_APREF_DIR = find_apref_dir(NGCA_DIR)
NGCA_TRANS_TABLES_DIR = NGCA_DIR / "transTables"

DYNADJUST_DIR = REPO_ROOT / "dynadjust-linux-mkl"