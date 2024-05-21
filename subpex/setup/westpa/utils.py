import os


def link_file(
    path_original: str | None,
    link_name: str,
    write_dir: str = "",
    f_ext: str | None = None,
) -> None:
    if path_original is None:
        raise ValueError(f"Original final path for {link_name} cannot be None")
    if f_ext is None:
        f_ext = os.path.splitext(path_original)[1][1:]
    os.symlink(path_original, os.path.join(write_dir, f"{link_name}.{f_ext}"))
