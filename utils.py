from pathlib import Path
import shutil


def clear_tmp_dir(tmp_dir="data/tmp"):
    """
    Delete all files and subdirectories inside tmp_dir,
    but keep tmp_dir itself.
    """
    tmp_path = Path(tmp_dir)

    if not tmp_path.exists():
        print(f"[clear_tmp_dir] {tmp_dir} does not exist — nothing to delete.")
        return

    for item in tmp_path.iterdir():
        try:
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
        except Exception as e:
            print(f"[clear_tmp_dir] Failed to delete {item}: {e}")

    print(f"[clear_tmp_dir] Cleared contents of {tmp_dir}")
