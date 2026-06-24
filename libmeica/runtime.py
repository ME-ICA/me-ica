# libmeica/runtime.py
import os
import shlex
import shutil
import sys


def _prepend_env_path(varname: str, new_path: str) -> None:
    """Prepend a path to an environment variable without duplicate insertion."""
    if not new_path:
        return

    new_path_abs = os.path.abspath(os.path.expanduser(new_path))
    old_value = os.environ.get(varname, "")
    old_parts = [p for p in old_value.split(":") if p]
    old_parts = [
        p for p in old_parts if os.path.abspath(os.path.expanduser(p)) != new_path_abs
    ]
    os.environ[varname] = ":".join([new_path_abs] + old_parts)


def find_distributed_afni_runtime() -> str:
    """
    Locate the ME-ICA bundled AFNI runtime.

    Search order:
      1. MEICA_AFNI_RUNTIME environment variable
      2. runtime directory installed inside the ME-ICA source tree:
           me-ica/
             meica.py
             libmeica/
             meica-afni-runtime/
    """
    candidates = []

    env_runtime = os.environ.get("MEICA_AFNI_RUNTIME", "")
    if env_runtime:
        candidates.append(env_runtime)

    # libmeica/runtime.py -> libmeica -> me-ica-X.Y.Z
    this_file = os.path.abspath(os.path.realpath(__file__))
    libmeica_dir = os.path.dirname(this_file)
    meica_src_dir = os.path.dirname(libmeica_dir)

    candidates.append(os.path.join(meica_src_dir, "meica-afni-runtime"))

    for runtime in candidates:
        runtime = os.path.abspath(os.path.expanduser(runtime))
        runtime_bin = os.path.join(runtime, "bin")
        activate_script = os.path.join(runtime, "activate_meica_afni")

        if (
            os.path.isdir(runtime)
            and os.path.isdir(runtime_bin)
            and os.path.isfile(activate_script)
        ):
            return runtime

    return ""


def activate_distributed_afni_runtime(
    *,
    require: bool | None = None,
    verbose: bool = True,
) -> str:
    """
    Activate bundled AFNI runtime for this Python process and its children only.

    This intentionally does not alter the user's parent shell. It only changes
    os.environ inside meica.py, so AFNI calls launched by meica.py inherit the
    distributed runtime.
    """
    if require is None:
        require = os.environ.get("MEICA_REQUIRE_DISTRIBUTED_AFNI", "0") == "1"

    runtime = find_distributed_afni_runtime()

    if not runtime:
        if require:
            print(
                "*+ MEICA_REQUIRE_DISTRIBUTED_AFNI=1, but the bundled AFNI "
                "runtime could not be found. Expected MEICA_AFNI_RUNTIME or "
                "a sibling directory named meica-afni-runtime.",
                file=sys.stderr,
            )
            sys.exit(11)
        return ""

    runtime_bin = os.path.join(runtime, "bin")
    runtime_lib = os.path.join(runtime, "lib")
    runtime_dri = os.path.join(runtime_lib, "dri")
    runtime_activate = os.path.join(runtime, "activate_meica_afni")

    _prepend_env_path("PATH", runtime_bin)

    os.environ["MEICA_AFNI_RUNTIME"] = runtime
    os.environ["MEICA_AFNI_RUNTIME_ACTIVATE"] = runtime_activate

    os.environ["AFNI_PLUGINPATH"] = runtime_bin
    os.environ["AFNI_NOMOTD"] = "YES"

    # The runtime binaries/libs should use RPATH, not LD_LIBRARY_PATH.
    # Mesa is scoped here only for this process and child scripts.
    if os.path.isdir(runtime_dri):
        os.environ["LIBGL_DRIVERS_PATH"] = runtime_dri
        os.environ["MESA_LOADER_DRIVER_OVERRIDE"] = "swrast"
        os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"

    if verbose:
        print(f"++ Using ME-ICA distributed AFNI runtime: {runtime}")
        print(f"++ 3dinfo: {shutil.which('3dinfo')}")
        print(f"++ 3dSkullStrip: {shutil.which('3dSkullStrip')}")
        print(f"++ 3dSeg: {shutil.which('3dSeg')}")

    return runtime


def shell_activation_lines() -> list[str]:
    """
    Return shell lines to source the bundled AFNI runtime from generated scripts.

    These lines make _meica_*.sh explicit and reproducible when run later.
    """
    runtime_activate = os.environ.get("MEICA_AFNI_RUNTIME_ACTIVATE", "")

    if runtime_activate:
        return [
            "# Activate ME-ICA distributed AFNI runtime",
            f"source {shlex.quote(runtime_activate)}",
        ]

    runtime = find_distributed_afni_runtime()
    if runtime:
        runtime_activate = os.path.join(runtime, "activate_meica_afni")
        return [
            "# Activate ME-ICA distributed AFNI runtime",
            f"source {shlex.quote(runtime_activate)}",
        ]

    return [
        "# No ME-ICA distributed AFNI runtime found at script-generation time.",
        "# Falling back to AFNI on PATH.",
    ]


def runtime_summary() -> str:
    runtime = (
        os.environ.get("MEICA_AFNI_RUNTIME", "") or find_distributed_afni_runtime()
    )
    return "\n".join(
        [
            f"MEICA_AFNI_RUNTIME={runtime}",
            f"MEICA_AFNI_RUNTIME_ACTIVATE={os.environ.get('MEICA_AFNI_RUNTIME_ACTIVATE', '')}",
            f"3dinfo={shutil.which('3dinfo')}",
            f"3dSkullStrip={shutil.which('3dSkullStrip')}",
            f"3dSeg={shutil.which('3dSeg')}",
        ]
    )
