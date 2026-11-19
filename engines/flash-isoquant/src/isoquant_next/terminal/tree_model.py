"""Zero-XGBoost L1-cache packed float32 GBDT ensemble loader & inference wrapper."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from isoquant_next.terminal._cy_tree_eval import eval_gbtree_proba_cy

FEATURE_COLUMNS: List[str] = [
    "var",
    "skew",
    "peak_count",
    "peak_width",
    "entropy",
    "mean_height",
    "peak_heights",
    "relative_height",
]


def _parse_base_score(raw_val) -> float:
    if isinstance(raw_val, (int, float)):
        return float(raw_val)
    s = str(raw_val).strip()
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    return float(s)


def _logit(p: float) -> float:
    eps = 1e-12
    p = min(max(p, eps), 1.0 - eps)
    return math.log(p / (1.0 - p))


class L1TreeEnsemble:
    """Compact contiguous float32/int32 decision-tree ensemble loaded directly from XGBoost JSON."""

    __slots__ = (
        "model_path",
        "feature_names",
        "n_features",
        "num_raw_trees",
        "num_active_trees",
        "const_margin",
        "tree_roots",
        "left_child",
        "right_child",
        "feature_idx",
        "default_left",
        "split_val",
    )

    def __init__(self, model_path: Optional[Union[str, Path]] = None) -> None:
        self.model_path: Optional[Path] = None
        self.feature_names: List[str] = list(FEATURE_COLUMNS)
        self.n_features: int = len(FEATURE_COLUMNS)
        self.num_raw_trees: int = 0
        self.num_active_trees: int = 0
        self.const_margin: float = 0.0
        self.tree_roots = np.empty(0, dtype=np.int32)
        self.left_child = np.empty(0, dtype=np.int32)
        self.right_child = np.empty(0, dtype=np.int32)
        self.feature_idx = np.empty(0, dtype=np.int32)
        self.default_left = np.empty(0, dtype=np.uint8)
        self.split_val = np.empty(0, dtype=np.float32)
        if model_path is not None:
            self.load_model(model_path)

    def load_model(self, model_path: Union[str, Path]) -> "L1TreeEnsemble":
        path = Path(model_path)
        self.model_path = path
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)

        learner = data["learner"]
        lmp = learner["learner_model_param"]
        self.n_features = int(lmp.get("num_feature", len(FEATURE_COLUMNS)))
        self.feature_names = list(learner.get("feature_names") or FEATURE_COLUMNS)

        base_p = _parse_base_score(lmp.get("base_score", 0.5))
        margin_bias = _logit(base_p)

        trees = learner["gradient_booster"]["model"]["trees"]
        self.num_raw_trees = len(trees)

        roots: List[int] = []
        lefts: List[int] = []
        rights: List[int] = []
        fidxs: List[int] = []
        dlefts: List[int] = []
        svals: List[float] = []

        offset = 0
        for t in trees:
            lc = t["left_children"]
            rc = t["right_children"]
            fi = t["split_indices"]
            dl = t["default_left"]
            sc = t["split_conditions"]
            n_nodes = len(lc)
            if n_nodes == 1 and lc[0] == -1:
                # Pruned single-leaf tree: fold constant float32 leaf weight into const_margin
                margin_bias += float(np.float32(sc[0]))
                continue
            roots.append(offset)
            for i in range(n_nodes):
                lefts.append(-1 if lc[i] == -1 else offset + int(lc[i]))
                rights.append(-1 if rc[i] == -1 else offset + int(rc[i]))
                fidxs.append(int(fi[i]))
                dlefts.append(1 if dl[i] else 0)
                svals.append(float(sc[i]))
            offset += n_nodes

        self.const_margin = float(margin_bias)
        self.num_active_trees = len(roots)
        self.tree_roots = np.asarray(roots, dtype=np.int32)
        self.left_child = np.asarray(lefts, dtype=np.int32)
        self.right_child = np.asarray(rights, dtype=np.int32)
        self.feature_idx = np.asarray(fidxs, dtype=np.int32)
        self.default_left = np.asarray(dlefts, dtype=np.uint8)
        self.split_val = np.asarray(svals, dtype=np.float32)
        return self

    def _prepare_matrix(self, X) -> np.ndarray:
        if hasattr(X, "to_numpy"):
            arr = X.to_numpy(dtype=np.float32, copy=False)
        else:
            arr = np.asarray(X, dtype=np.float32)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if not arr.flags.c_contiguous:
            arr = np.ascontiguousarray(arr, dtype=np.float32)
        return arr

    def predict_proba_positive(self, X) -> np.ndarray:
        if self.num_raw_trees == 0:
            arr = self._prepare_matrix(X)
            return np.full(arr.shape[0], 0.5, dtype=np.float64)
        arr = self._prepare_matrix(X)
        if arr.shape[0] == 0:
            return np.empty(0, dtype=np.float64)
        return eval_gbtree_proba_cy(
            arr,
            self.tree_roots,
            self.left_child,
            self.right_child,
            self.feature_idx,
            self.default_left,
            self.split_val,
            self.const_margin,
        )

    def predict_proba(self, X) -> np.ndarray:
        p1 = self.predict_proba_positive(X)
        return np.column_stack((1.0 - p1, p1))

    def predict(self, X) -> np.ndarray:
        p1 = self.predict_proba_positive(X)
        return (p1 > 0.5).astype(np.int32)


def reference_eval_json_tree_py(json_path: Union[str, Path], X) -> np.ndarray:
    """Pure-Python XGBoost JSON reference traversal for <1e-6 mathematical verification."""
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    learner = data["learner"]
    base_p = _parse_base_score(learner["learner_model_param"].get("base_score", 0.5))
    base_margin = _logit(base_p)
    trees = learner["gradient_booster"]["model"]["trees"]

    arr = np.asarray(X, dtype=np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    probs = np.empty(arr.shape[0], dtype=np.float64)
    for r in range(arr.shape[0]):
        margin = base_margin
        for t in trees:
            lc = t["left_children"]
            rc = t["right_children"]
            fi = t["split_indices"]
            dl = t["default_left"]
            sc = t["split_conditions"]
            node = 0
            while lc[node] != -1:
                fval = float(arr[r, fi[node]])
                if math.isnan(fval):
                    node = lc[node] if dl[node] else rc[node]
                else:
                    thresh = float(np.float32(sc[node]))
                    node = lc[node] if fval < thresh else rc[node]
            margin += float(np.float32(sc[node]))
        probs[r] = 1.0 / (1.0 + math.exp(-margin))
    return probs
