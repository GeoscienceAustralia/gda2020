from pathlib import Path

# Find repo root
def find_repo_root(path=None):
    current_dir = path or Path.cwd()
    
    for path in [current_dir, *current_dir.parents]:
        if (path / ".git").exists():
            return path
        
    raise RuntimeError(f"{current_dir} is not inside a git repository")

REPO_ROOT = find_repo_root()

# ---- Directories ------------------------------------------------------------

NGCA_DIR = REPO_ROOT / "ngca"
NGCA_