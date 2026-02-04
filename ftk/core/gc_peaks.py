from typing import List, Optional


def gc_peaks_kde(gc_values: List[float], grid: int = 1024, bw_method: Optional[float] = None,
                 min_prominence: float = 0.01) -> int:
    if len(gc_values) < 5:
        return 0

    try:
        from scipy.stats import gaussian_kde
    except Exception as e:
        raise RuntimeError("scipy is required for KDE-based gc_peaks.") from e

    xs = [i / (grid - 1) for i in range(grid)]
    kde = gaussian_kde(gc_values, bw_method=bw_method)
    ys = kde(xs)

    y_max = max(ys) if ys is not None else 0.0
    if y_max <= 0:
        return 0
    ys = [y / y_max for y in ys]

    peaks = 0
    for i in range(1, grid - 1):
        if ys[i] > ys[i - 1] and ys[i] > ys[i + 1] and ys[i] >= min_prominence:
            peaks += 1
    return peaks


def gc_peaks_gmm(gc_values: List[float], max_components: int = 6) -> int:
    if len(gc_values) < 10:
        return 0

    try:
        from sklearn.mixture import GaussianMixture
    except Exception as e:
        raise RuntimeError("scikit-learn is required for GMM-based gc_peaks.") from e

    import numpy as np
    X = np.array(gc_values, dtype=float).reshape(-1, 1)

    best_k = 1
    best_bic = None
    for k in range(1, max_components + 1):
        gmm = GaussianMixture(n_components=k, covariance_type="full", random_state=0)
        gmm.fit(X)
        bic = gmm.bic(X)
        if best_bic is None or bic < best_bic:
            best_bic = bic
            best_k = k
    return best_k
