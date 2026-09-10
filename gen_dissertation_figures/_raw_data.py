"""Reusable loading and preprocessing of raw LIUQE ``.mat`` samples."""
from dataclasses import dataclass
from pathlib import Path

from ._common import *

MIN_NUMBER_NONZERO_IY = 240
MAX_DELTA_FX, MIN_FX, MAX_FX = 0.8, -2.0, 2.0
RAW_NAMES = INPUT_NAMES + [IP, FX, IY]


@dataclass(frozen=True)
class RawSample:
    """One raw time slice converted to the inputs and full-grid targets used by the network figure."""
    shot: str
    time_index: int
    time: float
    x: np.ndarray
    fg: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def _raw_array(data, name, nt):
    expected = (*DS_SIZES[name], nt)
    if name not in data: raise KeyError(f'Missing raw variable {name!r}.')
    try: value = np.asarray(data[name]).reshape(expected)
    except ValueError as error: raise ValueError(f'Raw variable {name!r} cannot be reshaped to {expected}: {np.asarray(data[name]).shape}.') from error
    return value


def raw_valid_mask(data):
    """Return the per-time validity mask used by ``prepare_dataset.ipynb``."""
    nt = np.asarray(data['t']).size
    iy, fx = _raw_array(data, IY, nt), _raw_array(data, FX, nt)
    iy_ok = np.sum(np.abs(iy) > 1e-2, axis=(0, 1)) > MIN_NUMBER_NONZERO_IY
    ip_ok = np.abs(np.asarray(data[IP]).reshape(nt)) > 30_000
    flat_fx = fx.reshape(-1, nt); delta_fx = np.max(flat_fx, axis=0) - np.min(flat_fx, axis=0)
    fx_ok = (delta_fx < MAX_DELTA_FX) & (np.max(flat_fx, axis=0) < MAX_FX) & (np.min(flat_fx, axis=0) > MIN_FX)
    return iy_ok & ip_ok & fx_ok & np.isfinite(flat_fx).all(axis=0)


def flux_to_full_grid_fields(fx, k_finer=7, interp_method=INTERP_METHOD):
    """Derive $B_r$ and $B_z$ on a configurable intermediate grid; defaults reproduce the notebook."""
    if fx.shape != RRD.shape: raise ValueError(f'Flux shape must be {RRD.shape}, got {fx.shape}.')
    if k_finer < 1: raise ValueError(f'k_finer must be at least 1, got {k_finer}.')
    if k_finer == 1: return tuple(value.reshape(RRD.shape).astype(DTYPE) for value in meqBrBz(fx))
    rf = np.linspace(RRD[0, 0], RRD[0, -1], (RRD.shape[1]-1)*k_finer+1)
    zf = np.linspace(ZZD[0, 0], ZZD[-1, 0], (ZZD.shape[0]-1)*k_finer+1)
    rrf, zzf = np.meshgrid(rf, zf); fine_pts = np.stack((rrf, zzf), axis=-1)
    fine_fx = interp_pts(fx, fine_pts.reshape(-1, 2), method=interp_method).reshape(rrf.shape)
    fine_br, fine_bz = [value.reshape(rrf.shape) for value in meqBrBz(fine_fx, rr=rrf, zz=zzf)]
    grid_pts = np.stack((RRD, ZZD), axis=-1).reshape(-1, 2)
    br = interp_pts(fine_br, grid_pts, gr=rf, gz=zf, method='linear').reshape(RRD.shape)
    bz = interp_pts(fine_bz, grid_pts, gr=rf, gz=zf, method='linear').reshape(RRD.shape)
    return br.astype(DTYPE), bz.astype(DTYPE)


def preprocess_raw_sample(data, shot, time_index, require_valid=True, field_k_finer=7, field_interp_method=INTERP_METHOD):
    """Convert one loaded raw time slice to a network input and notebook-equivalent full-grid maps."""
    nt = np.asarray(data['t']).size
    if not -nt <= time_index < nt: raise IndexError(f'Time index {time_index} is outside a shot with {nt} samples.')
    time_index %= nt
    if require_valid and not raw_valid_mask(data)[time_index]: raise ValueError(f'Shot {shot}, time index {time_index} is rejected by the dataset filters.')
    values = {name:_raw_array(data, name, nt)[...,time_index].astype(DTYPE) for name in INPUT_NAMES + [FX, IY]}
    for name, value in values.items():
        if not np.isfinite(value).all(): raise ValueError(f'Shot {shot}, time index {time_index}: {name!r} contains non-finite values.')
    x = np.concatenate([values[name].reshape(-1) for name in INPUT_NAMES]).astype(DTYPE)
    iy = np.zeros_like(values[FX]); iy[1:-1,1:-1] = values[IY]
    br, bz = flux_to_full_grid_fields(values[FX], k_finer=field_k_finer, interp_method=field_interp_method)
    return RawSample(str(shot), time_index, float(np.asarray(data['t']).reshape(-1)[time_index]), x, (values[FX], iy, br, bz))


def load_raw_sample(path, time_index, require_valid=True, field_k_finer=7, field_interp_method=INTERP_METHOD):
    """Load and preprocess one time slice directly from a raw MATLAB shot file."""
    path = Path(path)
    if not path.is_file(): raise FileNotFoundError(f'Raw shot file not found: {path}')
    data = loadmat(path, variable_names=['t', *RAW_NAMES])
    return preprocess_raw_sample(data, path.stem, time_index, require_valid=require_valid, field_k_finer=field_k_finer, field_interp_method=field_interp_method)


def load_raw_samples(raw_dir, selections, require_valid=True, field_k_finer=7, field_interp_method=INTERP_METHOD):
    """Load ``(shot, time_index)`` selections from ``raw_dir`` in the requested order."""
    raw_dir = Path(raw_dir)
    return [load_raw_sample(raw_dir/f'{shot}.mat', time_index, require_valid=require_valid, field_k_finer=field_k_finer, field_interp_method=field_interp_method) for shot, time_index in selections]
