import atexit
import functools
import traceback
import inspect
import sys
import os
import gc
import uuid
import numpy as np
from types import ModuleType, FunctionType, BuiltinFunctionType, SimpleNamespace
from collections.abc import Mapping, Container
from contextlib import contextmanager
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None


def format_size(obj):
    try:
        if np and isinstance(obj, np.ndarray):
            return f"{obj.nbytes / 1024**2:.2f} MB (ndarray shape={obj.shape})"
        elif hasattr(obj, "__len__") and hasattr(obj, "__getitem__"):
            return f"{sys.getsizeof(obj) / 1024**2:.2f} MB (len={len(obj)})"
        else:
            return f"{sys.getsizeof(obj) / 1024**2:.2f} MB"
    except Exception as e:
        return f"Could not measure size: {e}"


def trace_and_profile_inputs(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print("\n=== [Memory Profile + Call Trace] ===")
        print("→ Large object inventory:")
        # list_big_live_objects(min_size_mb=1000)
        print("→ Call stack (last 4 calls):")
        for line in traceback.format_stack()[-5:-1]:
            print(line.strip())

        sig = inspect.signature(func)
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()

        print("\n→ Input sizes:")
        for name, val in bound_args.arguments.items():
            print(f"  {name}: {format_size(val)}")

        print("=====================================\n")
        return func(*args, **kwargs)

    return wrapper


def deep_getsizeof(obj, seen):
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    seen.add(obj_id)

    size = 0

    if isinstance(obj, np.ndarray):
        size += obj.nbytes
    else:
        try:
            size += sys.getsizeof(obj)
        except Exception:
            return 0  # skip if size can't be obtained

    # Don't recurse into non-containers or irrelevant types
    if isinstance(
        obj,
        (str, bytes, bytearray, FunctionType, BuiltinFunctionType, ModuleType, type),
    ):
        return size

    # Traverse mappings (dict-like)
    if isinstance(obj, Mapping):
        for k, v in obj.items():
            size += deep_getsizeof(k, seen)
            size += deep_getsizeof(v, seen)

    # Traverse containers and SimpleNamespace
    elif isinstance(obj, (Container, SimpleNamespace)) and not isinstance(
        obj, (str, bytes, bytearray)
    ):
        try:
            items = (
                list(obj)
                if not isinstance(obj, SimpleNamespace)
                else vars(obj).values()
            )
            for item in items:
                size += deep_getsizeof(item, seen)
        except Exception:
            pass

    return size


def list_big_live_objects(limit=15, min_size_mb=1):
    seen = set()
    results = []
    # print("Scanning all live objects...")
    for obj in gc.get_objects():
        try:
            size = deep_getsizeof(obj, seen)
            if size >= min_size_mb * 1024 * 1024:
                results.append((type(obj).__name__, size, repr(obj)[:100]))
        except Exception:
            continue
    # Sort by size descending
    results.sort(key=lambda x: x[1], reverse=True)
    print(f"\nTop {limit} largest live objects (> {min_size_mb} MB):\n")
    for typename, size, preview in results[:limit]:
        print(f"{typename:20s} : {size / (1024**2):6.2f} MB  | {preview}")
    return results


def _unlink(p):
    p = Path(p)
    p.unlink()


def create_memmap(name, shape, dtype=np.float32, *, loc=None):
    try:
        memmap_dir = os.environ["MEMMAP_DIR"]
    except:
        if loc is not None:
            memmap_dir = str(loc)
        else:
            memmap_dir = os.path.join(os.getcwd(), "memmaps")
        os.makedirs(memmap_dir, exist_ok=True)
    mmpath = os.path.join(memmap_dir, str(uuid.uuid4()) + f"{name}.dat")
    atexit.register(_unlink, mmpath)
    return np.memmap(mmpath, dtype=dtype, mode="w+", shape=shape)


@contextmanager
def flush_after(*memmaps):
    try:
        yield
    finally:
        for m in memmaps:
            m.flush()
