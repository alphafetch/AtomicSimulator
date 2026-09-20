from pathlib import Path
import shutil

def clear_logs():
    path = Path("./")

    for item in path.glob('atoms-*.csv'):
        if item.is_file():
            item.unlink()

clear_logs()