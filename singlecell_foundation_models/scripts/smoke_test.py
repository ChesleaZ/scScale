#!/usr/bin/env python
from __future__ import annotations

import importlib
import os
import platform
import sys
from pathlib import Path


def check(module_name: str) -> None:
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        print(f"[FAIL] {module_name}: {exc}")
        return

    version = getattr(module, "__version__", "version unknown")
    print(f"[ OK ] {module_name}: {version}")


print(f"Python: {platform.python_version()}")

ROOT = Path(__file__).resolve().parents[1]

for name in [
    "torch",
    "scanpy",
    "anndata",
    "scvi",
    "transformers",
    "datasets",
    "huggingface_hub",
    "scgpt",
    "geneformer",
]:
    check(name)

try:
    from geneformer import Classifier, EmbExtractor, InSilicoPerturber, TranscriptomeTokenizer

    _ = (Classifier, EmbExtractor, InSilicoPerturber, TranscriptomeTokenizer)
    print("[ OK ] geneformer classes")
except Exception as exc:
    print(f"[FAIL] geneformer classes: {exc}")


def check_path(path: Path) -> None:
    if path.exists():
        size_gb = path.stat().st_size / (1024**3)
        print(f"[ OK ] {path.relative_to(ROOT)}: {size_gb:.2f} GB")
    else:
        print(f"[FAIL] missing {path.relative_to(ROOT)}")


cwd = Path.cwd()
try:
    os.chdir(ROOT / "sources" / "UCE")
    sys.path.insert(0, str(ROOT / "sources" / "UCE"))
    import model as uce_model  # noqa: F401
    from evaluate import AnndataProcessor  # noqa: F401

    print("[ OK ] UCE source imports")
    check_path(ROOT / "models" / "uce" / "model_files" / "4layer_model.torch")
    check_path(ROOT / "models" / "uce" / "model_files" / "all_tokens.torch")
except Exception as exc:
    print(f"[FAIL] UCE source imports: {exc}")
finally:
    os.chdir(cwd)

try:
    os.chdir(ROOT / "sources" / "scFoundation" / "model")
    sys.path.insert(0, str(ROOT / "sources" / "scFoundation" / "model"))
    import load as scfoundation_load  # noqa: F401
    import get_embedding as scfoundation_get_embedding  # noqa: F401

    print("[ OK ] scFoundation source imports")
    check_path(ROOT / "models" / "scfoundation" / "models.ckpt")
except Exception as exc:
    print(f"[FAIL] scFoundation source imports: {exc}")
finally:
    os.chdir(cwd)

try:
    import torch

    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"CUDA device count: {torch.cuda.device_count()}")
        print(f"CUDA device 0: {torch.cuda.get_device_name(0)}")
except Exception as exc:
    print(f"[WARN] CUDA check failed: {exc}")
