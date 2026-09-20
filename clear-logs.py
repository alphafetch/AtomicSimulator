from pathlib import Path
import shutil

def clear_logs():
    path = Path("logs/")

    for item in path.iterdir():
        try:
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        except Exception as e:
            print(f"Failed to delete {item}. Reason: {e}")

clear_logs()