"""Reusable flux-based importance sampling for dissertation dataset figures."""
from ._common import *


def sample_uniform_vessel(n_points, rng=None):
    """Uniformly sample points inside the TCV first-wall polygon."""
    if n_points < 1: raise ValueError(f'n_points must be positive, got {n_points}.')
    rng = np.random.default_rng() if rng is None else rng
    points = np.empty((n_points, 2), dtype=DTYPE); n_sampled = 0
    while n_sampled < n_points:
        n_remaining = n_points-n_sampled
        candidates = np.column_stack((rng.uniform(RRD[0,0], RRD[0,-1], 2*n_remaining), rng.uniform(ZZD[0,0], ZZD[-1,0], 2*n_remaining)))
        selected = candidates[inpoly(candidates, VESS_CORNERS)][:n_remaining]
        points[n_sampled:n_sampled+len(selected)] = selected; n_sampled += len(selected)
    return points


def sample_flux_importance(n_points, flux, rq, zq, sigma=0.3, rng=None):
    """Sample vessel points using the flux-distance acceptance rule from ``prepare_dataset.ipynb``."""
    if n_points < 1: raise ValueError(f'n_points must be positive, got {n_points}.')
    if sigma <= 0: raise ValueError(f'sigma must be positive, got {sigma}.')
    rng = np.random.default_rng() if rng is None else rng
    rq = np.clip(rq, RRD[0,0], RRD[0,-1]); zq = np.clip(zq, ZZD[0,0], ZZD[-1,0])
    boundary_flux = np.mean(interp_pts(flux, np.stack((rq, zq), axis=-1), method='linear'))
    flux_range = np.max(flux)-np.min(flux)
    if not np.isfinite(flux_range) or flux_range <= 0: raise ValueError('Flux map must have a positive finite range.')
    points = np.empty((n_points, 2), dtype=DTYPE); n_sampled = 0
    for _ in range(200):
        n_remaining = n_points-n_sampled
        candidates = np.column_stack((rng.uniform(RRD[0,0], RRD[0,-1], 3*n_remaining), rng.uniform(ZZD[0,0], ZZD[-1,0], 3*n_remaining)))
        flux_distance = np.abs(interp_pts(flux, candidates, method='linear')-boundary_flux)/flux_range
        accepted = (rng.normal(0, sigma, len(candidates)) > flux_distance) & inpoly(candidates, VESS_CORNERS)
        selected = candidates[accepted][:n_remaining]; points[n_sampled:n_sampled+len(selected)] = selected; n_sampled += len(selected)
        if n_sampled == n_points: return points
    raise RuntimeError(f'Could not draw {n_points} importance samples within 200 iterations; increase sigma.')
