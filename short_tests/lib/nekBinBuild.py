import os
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path


def build_tools(
    tools_root,
    tools_bin,
    f77=None,
    cc=None,
    bigmem=None,
    targets=("clean", "all"),
    verbose=False,
):

    tools_root = Path(tools_root)

    print("Compiling tools... ")
    print(f'    Using output directory "{tools_bin}"')
    print(f'    Using FC "{f77}"')
    print(f'    Using CC "{cc}"')

    maketools_in = tools_root / "maketools"

    my_env = os.environ.copy()
    if f77:
        my_env["FC"] = f77
    if cc:
        my_env["CC"] = cc
    my_env["bin_nek_tools"] = tools_bin

    if targets[0] == "core":
        with Popen(
            [maketools_in, "core"], env=my_env, cwd=tools_root, stderr=STDOUT
        ) as proc:
            proc.wait()
        if proc.returncode != 0:
            targets = [t for t in os.listdir(tools_root) if "maketools" not in t]
            for t in targets:
                logfile = tools_root / t / "build.log"
                print(logfile, end="")
                try:
                    with open(logfile, "r") as file:
                        text = file.read()
                    print(":")
                    print(text)
                except FileNotFoundError:
                    print("  (File does not exist)")
            exit(-1)
        return

    if targets[0] == "all":
        targets = [t for t in os.listdir(tools_root) if "maketools" not in t]
        print("Targets:", targets)

    for t in targets:
        with Popen(
            [maketools_in, t], env=my_env, cwd=tools_root, stderr=STDOUT
        ) as proc:
            proc.wait()
        logfile = tools_root / t / "build.log"
        if proc.returncode != 0:
            with open(logfile, "r") as file:
                text = file.read()
            print(text)
            exit(-1)


def build_nek(source_root, usr_file, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(NEK_SOURCE_ROOT=source_root)

    print("Compiling nek5000...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    for key, val in list(_opts.items()):
        print(f'    Using {key}="{val}"')

    my_env = os.environ.copy()
    if source_root:
        my_env["NEK_SOURCE_ROOT"] = source_root
    if _opts.get("F77"):
        my_env["FC"] = _opts.get("F77")
    if _opts.get("CC"):
        my_env["CC"] = _opts.get("CC")
    if _opts.get("PPLIST"):
        my_env["PPLIST"] = _opts.get("PPLIST")
    if _opts.get("FFLAGS"):
        my_env["FFLAGS"] = _opts.get("FFLAGS")

    makenek_in = Path(source_root) / "bin" / "makenek"
    logfile = Path(cwd) / "build.log"

    with Popen(
        [makenek_in, "clean"], cwd=cwd, env=my_env, stdin=PIPE, text=True
    ) as proc:
        proc.communicate(input="N\n")

    with Popen(
        [makenek_in, usr_file], cwd=cwd, env=my_env, stderr=STDOUT
    ) as proc:
        proc.wait()

    if proc.returncode != 0:
        with open(logfile, "r") as file:
            text = file.read()
        print(text)
        exit(-1)


def build_libtest(source_root, driver_path, cwd=None, f77=None, verbose=False):
    """Link a caller-supplied Fortran driver against libnek5000.a only
    (no separate *_mod.o), exercising the library/embedding path:
      make libtest LIBTEST_DRV=<driver.f>
    Produces ./nek_libtest in cwd. Exits nonzero on link failure so the
    calling test fails.
    """
    print("Linking library test (make libtest)...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using driver "{driver_path}"')

    my_env = os.environ.copy()
    if source_root:
        my_env["NEK_SOURCE_ROOT"] = source_root
    if f77:
        my_env["FC"] = f77

    logfile = Path(cwd) / "libtest_build.log"
    with open(logfile, "w") as f:
        with Popen(
            ["make", "libtest", f"LIBTEST_DRV={driver_path}"],
            cwd=cwd,
            env=my_env,
            stdout=f,
            stderr=STDOUT,
        ) as proc:
            proc.wait()

    if proc.returncode != 0:
        with open(logfile, "r") as file:
            print(file.read())
        exit(-1)
