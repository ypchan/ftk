from __future__ import annotations

from enum import Enum
from typing import List, Optional


class GcPeaksMethod(str, Enum):
    kde = "kde"
    gmm = "gmm"


def count_gc_peaks(
    gc_values: List[float],
    method: GcPeaksMethod,
    *,
    kde_grid: int = 1024,
    kde_min_prominence: float = 0.01,
    gmm_max_components: int = 6,
) -> Optional[int]:
    if not gc_values:
        return None

    if method == GcPeaksMethod.kde:
        return _count_kde(gc_values, grid=kde_grid, min_prominence=kde_min_prominence)
    if method == GcPeaksMethod.gmm:
        return _count_gmm(gc_values, max_components=gmm_max_components)

    raise ValueError(f"Unknown method: {method}")


def _count_kde(gc_values: List[float], *, grid: int, min_prominence: float) -> int:
    import numpy as np
    from scipy.stats import gaussian_kde

    x = np.array(gc_values, dtype=float)
    if x.size < 3:
        return 1

    xs = np.linspace(0.0, 1.0, grid)
    ys = gaussian_kde(x)(xs)

    max_y = float(ys.max()) if ys.size else 0.0
    if max_y <= 0:
        return 1

    peaks = 0
    for i in range(1, len(xs) - 1):
        if ys[i] > ys[i - 1] and ys[i] > ys[i + 1]:
            if float(ys[i] / max_y) >= min_prominence:
                peaks += 1
    return max(peaks, 1)


def _count_gmm(gc_values: List[float], *, max_components: int) -> int:
    import numpy as np
    from sklearn.mixture import GaussianMixture

    x = np.array(gc_values, dtype=float).reshape(-1, 1)
    if x.shape[0] < 3:
        return 1

    best_k = 1
    best_bic = None
    for k in range(1, max_components + 1):
        model = GaussianMixture(n_components=k, random_state=0)
        model.fit(x)
        bic = model.bic(x)
        if best_bic is None or bic < best_bic:
            best_bic = bic
            best_k = k
    return best_k
