# utils functions
DARK_THEME = False # dark theme for plotting

import numpy as np
np.set_printoptions(precision=2)
from numpy.random import uniform

## just for plotting settings
import matplotlib.pyplot as plt
plt.style.use('dark_background' if DARK_THEME else 'default')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3 if DARK_THEME else 0.7
plt.rcParams['font.family'] = 'monospace'
# plt.rcParams['image.cmap'] = 'inferno'
plt.rcParams['image.cmap'] = 'plasma'
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.colormaps['plasma'](np.array([0.2,0.4,0.6,0.8,1.0])))
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.colormaps['plasma'](np.array([.1,.5,.2,.6,.3,.7,.4,.8,.9,1.0])))

# print("Current savefig.dpi:", plt.rcParams['savefig.dpi'])
# plt.rcParams['savefig.dpi'] = 300

from scipy.io import loadmat, savemat

import torch
from torch.nn import functional as F, Module, Linear, Conv2d, MaxPool2d, BatchNorm2d, ReLU, Sequential, ConvTranspose2d, Tanh
from torch.utils.data import Dataset

import builtins
import os

CUDA = torch.device('cuda')
CPU = torch.device('cpu')
MPS = torch.device('mps')

DEV = CPU
try: 
    JOBID = os.environ["SLURM_JOB_ID"] # get job id from slurm, when training on cluster
    NCORES = int(os.environ.get("SLURM_CPUS_PER_TASK", 1)) # number of cores per task
    if torch.cuda.is_available():
        DEV = CUDA
        GPU_MEM, _ = torch.cuda.mem_get_info()
    LOCAL = False # for plotting or saving images
except:
    if torch.backends.mps.is_available():
        DEV = MPS
        GPU_MEM = 16 * 1024**3  # 16 GB in bytes
    JOBID = "local"
    LOCAL = True
    NCORES = os.cpu_count() or 1
print(f"Running JOBID: {JOBID}, on {DEV}, GPU_MEM: {GPU_MEM/1e9:.2f} GB" if DEV != CPU else f"Running JOBID: {JOBID}, on cpu. " + f"NCORES: {NCORES}")
SAVE_DIR = f"data/{JOBID}"

N_CTRL_PTS = 25 # number of control points to estimate

NDOM = 4 # number of plasma domains, should be matched with create_ds.m
##################################################################################################################################
## Dataset inputs/outputs (see create_ds.m)
# input names
BM = 'Bm' # magnetic probe measurements
FF = 'Ff' # flux loops poloidal flux 
IA = 'Ia' # poloidal field coil currents
IP = 'Ip' # plasma current % not used
UF = 'Uf' # 
RBT = 'rBt' #

# output names
FX = 'Fx' # Plasma poloidal flux map
IY = 'Iy' # Plasma current density map
BR = 'Br' # recalculated in python, not in matlab
BZ = 'Bz' # recalculated in python, not in matlab
#more names
PHYS = 'phys' # physics inputs
PTS = 'pts' # points 
FG = 'fg' # full grid
RT = 'rt' # real time output
B = 'B' # magnetic field, combined Br, Bz 
BD = 'Bd' # combined Br, Bz from derivative
IYD = 'Iyd' # current density from derivative

INPUT_NAMES = [BM, FF, IA, UF, RBT] # input names
OUTPUT_NAMES = [FX, IY, BR, BZ] # output names
DS_NAMES = INPUT_NAMES + OUTPUT_NAMES # dataset names
DS_SIZES = { 
    BM:(38,), FF:(38,), IA:(19,), UF:(38,), RBT:(1,),  # input sizes
    FX:(65,28), IY:(63,26), BR:(65,28), BZ:(65,28),  # output sizes
}
INPUT_MASKS = {
    name: np.arange(sum(DS_SIZES[n][0] for n in INPUT_NAMES[:i]), 
                    sum(DS_SIZES[n][0] for n in INPUT_NAMES[:i+1]))
    for i, name in enumerate(INPUT_NAMES)
}
DTYPE = 'float32'

π = torch.pi # constant pi
μ0 = 4e-7 * π  # Permeability of free space

from scipy.interpolate import RegularGridInterpolator
# INTERP_METHOD = 'linear' # fast, but less accurate
# INTERP_METHOD = 'quintic' # slowest, but most accurate
INTERP_METHOD = 'cubic' # slower, close accuracy to quintic, but little bit faster than quintic
if INTERP_METHOD == 'linear': print('Warning: using linear interpolation, which is fast but less accurate')

# DS_DIR = 'dss/ds' # where the final dataset will be stored
DS_DIR = 'dss' if LOCAL else '/nfsd/automatica/grandinmat/dss' 
os.makedirs(f'{DS_DIR}/imgs', exist_ok=True)
TRAIN_DS_PATH = f'{DS_DIR}/train_ds.npz'
EVAL_DS_PATH = f'{DS_DIR}/eval_ds.npz'

# paths to the best models
LOSS_NAMES = [FX, IY, B, BD, IYD] # loss names 
def model_path(loss_name, save_dir=SAVE_DIR):
    assert loss_name in LOSS_NAMES, f"loss_name should be one of {LOSS_NAMES}, got {loss_name}"
    return f"{save_dir}/best_{loss_name}.pth"

BEST_MODEL_DIR = 'data/best' # where the best models are saved
STRICT_LOAD = True # for loading the weights, should be true, but for testing varying architectures, set to false

TEST_DIR = 'test' if LOCAL else '/nfsd/automatica/grandinmat/test'
os.makedirs(TEST_DIR, exist_ok=True)

NIN = sum(DS_SIZES[name][0] for name in INPUT_NAMES) # 134

##################################################################################################################################


# read the original grid coordinates
d = loadmat('tcv_params/grid.mat')
rd, zd = d['r'].flatten(), d['z'].flatten() # original grid coordinates (DATA)
RRD, ZZD = np.meshgrid(rd, zd)  # meshgrid for the original grid coordinates (from the data)
RRD2, ZZD2 = np.meshgrid(rd[1:-1], zd[1:-1])  
GRA = (rd[1]-rd[0])*(zd[1]-zd[0])  # grid area
del d, rd, zd


# load the vessel perimeter
m = loadmat('tcv_params/vess.mat')
vr, vz = m['vess_r'], m['vess_z']
vri, vzi = m['r_in'], m['z_in']
vro, vzo = m['r_out'], m['z_out']
v = np.hstack([vr, vz])
vi, vo = np.hstack([vri, vzi]), np.hstack([vro, vzo])
VESS = np.vstack([v, v[0]])[::-1]
VESSI = np.vstack([vi, vi[0]])
VESSO = np.vstack([vo, vo[0]])
# _vess_corners_idxs = [92, 101, 140, 165, 201, 211, 220, 338, 372, 411, 420]
# VESS_CORNERS = VESS[_vess_corners_idxs]  # corners of the vessel, used for plotting
_vess_corners_idxs = [39,40,41,42,52,53,54,55,76,77,78,79,80,81,82,83,172,173,174,175,176,177,178,179,200,201,202,203,213,214,215,216]
VESS_CORNERS = VESSI[_vess_corners_idxs]  # corners of the vessel, used for plotting

def inpoly(pts, poly=VESS_CORNERS):
    pts_shape = pts.shape
    assert pts_shape[-1] == 2, f"pts.shape = {pts_shape}, should be (whatever, 2)"
    pts = pts.reshape(-1, 2)  # ensure pts is a 2D array of shape (n, 2)
    # Polygon edges: (xi, yi) to (xj, yj)
    xi, yi = poly[:, 0], poly[:, 1]
    xj, yj = np.roll(xi, -1), np.roll(yi, -1)
    px = pts[:, 0][:, np.newaxis]
    py = pts[:, 1][:, np.newaxis]
    # Check if point is between yi and yj, and if it is to the left of the edge
    intersect = ((yi <= py) & (py < yj)) | ((yj <= py) & (py < yi))
    slope = (xj - xi) / (yj - yi + 1e-12)  # avoid division by zero
    xints = xi + (py - yi) * slope
    crosses = px < xints
    inside = np.sum(intersect & crosses, axis=1) % 2 == 1
    return inside.reshape(pts_shape[:-1])  # reshape back to original shape

IV = inpoly(np.stack([RRD[:,:,np.newaxis], ZZD[:,:,np.newaxis]], axis=-1)).reshape(RRD.shape)  # inside vessel mask
del m, vr, vz, vri, vzi, vro, vzo, v, vi, vo

if not LOCAL: # Redefine the print function to always flush
    def print(*args, **kwargs): builtins.print(*args, **{**kwargs, 'flush': True})

####################################################################################################
## torch stuff
def to_tensor(x, device=torch.device(CPU)): return torch.tensor(x, dtype=torch.float32, device=device)

# custom trainable swish activation function
class ActF(Module): # swish
    def __init__(self, no=1): 
        super(ActF, self).__init__()
        if no == 1: self.beta = torch.nn.Parameter(torch.tensor(1.0), requires_grad=True)
        else: self.beta = torch.nn.Parameter(torch.ones(no), requires_grad=True) # n is the number of outputs
    def forward(self, x): return x*torch.sigmoid(self.beta*x)

# class ActF(Module): # swish beta = 1
#     def __init__(self, no=1): super(ActF, self).__init__()
#     def forward(self, x): return x*torch.sigmoid(x)


# PHYSICS_LS, PTSH, PHYSH, HEADH = 64, 64, 128, 128 # medium (old baseline)
# PHYSICS_LS, PTSH, PHYSH, HEADH = 48, 48, 64, 64 
# PHYSICS_LS, PTSH, PHYSH, HEADH = 48, 48, 48, 48 
PHYSICS_LS, PTSH, PHYSH, HEADH = 32, 32, 48, 48 # small (new baseline) <- 
# PHYSICS_LS, PTSH, PHYSH, HEADH = 32, 32, 32, 32 
# PHYSICS_LS, PTSH, PHYSH, HEADH = 16, 16, 32, 32 # turbosmall 
# PHYSICS_LS, PTSH, PHYSH, HEADH = 3, 16, 32, 32 # for shot mapping
print(f"Using PHYSICS_LS = {PHYSICS_LS}, PTSH = {PTSH}, PHYSH = {PHYSH}, HEADH = {HEADH}")

USE_BIAS = True # whether to use bias in the linear layers
print(f"BIAS = {USE_BIAS}")

class liuqeML(Module):
    def __init__(self, x_mean_std=torch.tensor([[0.0]*NIN, [1.0]*NIN], dtype=torch.float32)):
        super(liuqeML, self).__init__()
        self.pts_enc = Sequential(
            Linear(2, PTSH, bias=USE_BIAS), ActF(PTSH),
            # Linear(PTSH, PTSH, bias=USE_BIAS), ActF(PTSH),
            Linear(PTSH, PHYSICS_LS, bias=USE_BIAS), ActF(PHYSICS_LS),
        )
        self.x_mean_std = x_mean_std
        self.phys_enc = Sequential(
            Linear(NIN, PHYSH, bias=USE_BIAS), ActF(PHYSH),
            # Linear(PHYSH, PHYSH, bias=USE_BIAS), ActF(PHYSH),
            Linear(PHYSH, PHYSICS_LS, bias=USE_BIAS), ActF(PHYSICS_LS),
        )
        self.fx_head = Sequential(
            Linear(PHYSICS_LS, HEADH, bias=USE_BIAS), ActF(HEADH),
            # Linear(HEADH, HEADH, bias=USE_BIAS), ActF(HEADH),
            Linear(HEADH, 1, bias=USE_BIAS), ActF(1),
        )
        self.brz_head = Sequential(
            Linear(PHYSICS_LS, HEADH, bias=USE_BIAS), ActF(HEADH),
            # Linear(HEADH, HEADH, bias=USE_BIAS), ActF(HEADH),
            Linear(HEADH, 2, bias=USE_BIAS), ActF(2),
        )
        self.iy_head = Sequential(
            Linear(PHYSICS_LS, HEADH, bias=USE_BIAS), ActF(HEADH),
            # Linear(HEADH, HEADH, bias=USE_BIAS), ActF(HEADH),
            Linear(HEADH, 1, bias=USE_BIAS), ActF(1),
        )
    def to(self, device):
        super(liuqeML, self).to(device)
        self.x_mean_std = self.x_mean_std.to(device)
        return self
    def forward(self, x, pts):
        assert pts.shape[-1] == 2, f"pts.shape[-1] = {pts.shape[-1]}, should be 2 (r,z)"
        assert pts.dim() == 3, f"pts.dim() = {pts.dim()}, should be 3 (batch, n_points, 2)"
        assert x.dim() == 2, f"x.dim() = {x.dim()}, should be 2 (batch, NIN)"
        assert x.shape[-1] == NIN, f"x.shape[1] = {x.shape[1]}, NIN = {NIN}"
        assert pts.shape[0] == x.shape[0], f"pts.shape[0] = {pts.shape[0]}, x.shape[0] = {x.shape[0]}, should be equal (batch size)"
        npts = pts.shape[1] # number of points
        x = (x - self.x_mean_std[0]) / self.x_mean_std[1] # normalize
        # calculate physics vector ph and repeat it for each point in pts
        ph = self.phys_enc(x) # get physics vector 
        ph = ph.unsqueeze(1).expand(-1, npts, -1) # (BS, NP, PHYSICS_LS)
        # calculate a positional encoding for the points pts
        ps = self.pts_enc(pts) # encode points       
        psh = ps * ph # aggregate points and physics vector
        fx = self.fx_head(psh).view(-1, npts) # get flux
        brz = self.brz_head(psh).view(-1, npts, 2) # get Br, Bz
        iy = self.iy_head(psh).view(-1, npts) # get current density
        return fx, brz, iy

def test_network_io(verbose=True):
    v = verbose
    if v: print('test_network_io')
    n_points = 5 # number of points to sample
    # single sample
    x, pts = (torch.rand(1, NIN), torch.rand(1, n_points, 2)) # x: (1, NIN), pts: (1, n_points, 2)
    lml = liuqeML()
    fx, brz, iy = lml(x, pts)
    r, z = pts[0,:,0], pts[0,:,1] # r, z: (n_points,)
    assert fx.shape == (1, n_points), f"fx.shape = {fx.shape}, should be (1, {n_points})"
    assert brz.shape == (1, n_points, 2), f"brz.shape = {brz.shape}, should be (1, {n_points}, 2)"
    assert iy.shape == (1, n_points), f"iy.shape = {iy.shape}, should be (1, {n_points})"
    assert pts.shape == (1, n_points, 2), f"pts.shape = {pts.shape}, should be (1, {n_points}, 2)"
    if v: print(f"single liuqeML -> in: [{x.shape}, {pts.shape}], \n            out: [{fx.shape}, {brz.shape}, {iy.shape}]")
    # batched
    bs = 7
    x, pts = torch.rand(bs, NIN), torch.rand(bs, n_points, 2) # x: (bs, NIN), pts: (bs, n_points, 2)
    fx, brz, iy = lml(x, pts)
    assert iy.shape == (bs, n_points), f"iy.shape = {iy.shape}, should be ({bs}, {n_points})"
    assert fx.shape == (bs, n_points), f"fx.shape = {fx.shape}, should be ({bs}, {n_points})"
    assert brz.shape == (bs, n_points, 2), f"brz.shape = {brz.shape}, should be ({bs}, {n_points}, 2)"
    if v: print(f"batched liuqeML -> in: [{x.shape}, {pts.shape}], \n            out: [{fx.shape}, {brz.shape}, {iy.shape}]")
    # if v: print(f"No batched versions of LiuqeRTNet and LiuqeIyNet (meant for single sample inference only)")

def save_weights_to_mat(net, title="", save_dir=[SAVE_DIR]):
    weights = {}
    weights['x_mean_std'] = net.x_mean_std.cpu().numpy()  # save normalization parameters
    for name, param in net.named_parameters():
        w = param.data.cpu().numpy()
        n = name.replace('.', '_')
        weights[n] = w
        # print(f"{n}: {w.shape}")
    # print(f'Weights {title}: {weights.keys()}')
    for sd in save_dir:
        savemat(f'{sd}/weights_{title}.mat', weights)
        print(f"Weights saved to {sd}/weights_{title}.mat")

def save_weights_to_json(net, title="", save_dir=[SAVE_DIR]):
    import json
    weights = {}
    weights['x_mean_std'] = net.x_mean_std.cpu().numpy().tolist()  # save normalization parameters
    for name, param in net.named_parameters():
        w = param.data.cpu().numpy()
        n = name.replace('.', '_')
        weights[n] = w.tolist()
        # print(f"{n}: {w.shape}")
    # print(f'Weights {title}: {weights.keys()}')
    for sd in save_dir:
        with open(f'{sd}/weights_{title}.json', 'w') as f:
            json.dump(weights, f)
        print(f"Weights saved to {sd}/weights_{title}.json")

# function to load the dataset
def load_ds(ds_path):
    assert os.path.exists(ds_path), f"Dataset not found: {ds_path}"
    d = np.load(ds_path)
    r = {
        PHYS:d[PHYS], 
        PTS:d[PTS], 
        FX:d[FX], 
        # IY:np.abs(d[IY]), # TEST: make current dens is always positive 
        IY:d[IY]*μ0/GRA, # convert to a better unit 
        # IY:np.abs(d[IY])*μ0/GRA, # convert to a better unit 
        BR:d[BR], 
        BZ:d[BZ], 
        # SEP:d[SEP], 
        "x_mean_std":d["x_mean_std"],
        FG:d[FG] if FG in d else None, # full grid, optional
    }
    return r

####################################################################################################
class LiuqeDataset(Dataset):
    def __init__(self, ds_path, pts_per_sample=None, verbose=True):
        d = load_ds(ds_path)
        x_mean_std = to_tensor(d.pop("x_mean_std")) # mean and std for inputs
        self.fg = d.pop(FG) # pop the full grid 
        print(f"Dataset loaded from {ds_path}, keys: {d.keys()}")
        self.data = {k:to_tensor(v) for k,v in d.items()} # convert to tensors
        # move to DEV (doable bc the dataset is fairly small, check memory usage)
        tot_memory_ds = sum([x.element_size()*x.nelement() for x in self.data.values()]) # total memory in bytes
        gpu_free_mem = torch.cuda.mem_get_info()[0] if DEV == CUDA else np.inf
        self.on_dev = DEV != CPU and tot_memory_ds < ( gpu_free_mem - 3e9 ) # x gb buffer for batches
        if self.on_dev: self.data = {k:v.to(DEV) for k,v in self.data.items()} # move to device
        self.x_mean_std = x_mean_std.to(DEV) if self.on_dev else x_mean_std
        self.tnp = self.data[PTS].shape[1] # total number of points per timestamp
        self.pps = pts_per_sample if pts_per_sample is not None else self.tnp # number of points per sample
        assert self.tnp % self.pps == 0, f'Number of points per sample ({self.tnp}) is not a perfect multiple of pps ({self.pps})'
        if verbose: print(f"Dataset: n_t:{len(self.data[PHYS])}, total_pts_per_t:{self.tnp}, picked:{self.pps} -> N:{len(self)} ({len(self.data[PHYS])}*{self.tnp}//{self.pps}), memory:{tot_memory_ds/1e6}MB, on_dev:{self.on_dev}")
    def __len__(self): return len(self.data[PHYS]) * self.tnp // self.pps # number of samples, all data should have the same length
    def __getitem__(self, idx): 
        # return [x[idx] for x in self.data.values()] # [PHYS, PTS, FX, IY, BR, BZ]
        tnp, pps = self.tnp, self.pps
        phys_idx = idx // (tnp // pps) # physical sample index
        point_start = (idx % (tnp // pps)) * pps # starting point index
        point_end = point_start + pps # ending point index
        return (self.data[PHYS][phys_idx], *[self.data[k][phys_idx, point_start:point_end] for k in [PTS, FX, IY, BR, BZ]]) # return selected points only
def test_dataset(ds:LiuqeDataset, verbose=True):
    if verbose:
        print("test_dataset")
        print(f"Dataset length: {len(ds)}")
        print(f"Inputs: X -> {ds[0][0].shape}, r -> {ds[0][1].shape}, z -> {ds[0][2].shape}")
        print(f"Outputs: Y1 -> {ds[0][3].shape}, Y2 -> {ds[0][4].shape}, Y3 -> {ds[0][5].shape}")
    n_plot = 10
    print(len(ds))
    idxs = np.random.randint(0, len(ds), n_plot)
    plt.figure(figsize=(15, 3*n_plot))
    for i, j in enumerate(idxs):
        x, pts, fx, iy, br, bz = map(lambda x: x.cpu().numpy(), ds[j])
        μ, σ = ds.x_mean_std.cpu().numpy()
        x = (x-μ) / σ  # normalize inputs
        r,z = pts[:,0], pts[:,1]
        ms = max(1, int(1000/len(r)))  # marker size
        plt.subplot(n_plot, 5, i*5+1)
        plt.scatter(r, z, c=fx, s=ms), plt.title('FX')
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+2)
        plt.scatter(r, z, c=iy, s=ms), plt.title('IY')
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+3)
        plt.scatter(r, z, c=br, s=ms), plt.title('Br')
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+4)
        plt.scatter(r, z, c=bz, s=ms), plt.title('Bz')
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+5)
        plt.bar(np.arange(len(x)), x), plt.title('Physical Inputs')
        plt.tight_layout()
    plt.savefig(f"{TEST_DIR}/dataset_examples.png")

    return

####################################################################################################
def sample_random_points(n):
    pts = np.zeros((n, 2), dtype=DTYPE)
    pts[:,0] = uniform(RRD[0,0], RRD[0,-1], n) # r
    pts[:,1] = uniform(ZZD[0,0], ZZD[-1,0], n) # z
    return pts

def sample_random_points_inside_fw(n):
    over_n = int(n * 1.4)  # oversample to ensure we get enough points inside the vessel (Area ratio: 0.8895 ± 0.0081)
    pts = np.zeros((n, 2), dtype=DTYPE)
    cpts = sample_random_points(over_n)  # sample random points
    inside = inpoly(cpts, VESS_CORNERS)  # check if points are inside the vessel
    assert np.sum(inside) >= n, f"Not enough points inside the vessel, got {np.sum(inside)}, need {n}"
    pts = cpts[inside]  # keep only points inside the vessel
    return pts[:n]  # return the first n points

def uniform_sampling_inside_fw(n, gFx=None, rq=None, zq=None, σFx=None):
    return sample_random_points_inside_fw(n)

def sample_points_fx_importance_samp(n, gFx, rq, zq, σFx=0.01):
    tot_pts, iters, pts = 0, 0, np.zeros((n,2), dtype=DTYPE)  # initialize 
    rq = np.clip(rq, RRD[0,0], RRD[0,-1]) # clip to the grid range
    zq = np.clip(zq, ZZD[0,0], ZZD[-1,0])
    μFx = np.mean(interp_pts(gFx, np.stack((rq, zq), axis=-1), method='linear'))  # mean Fx at LCFS
    range_Fx = np.max(gFx) - np.min(gFx)
    while tot_pts < n:  # keep sampling until we have enough points or we reach the maximum number of iterations
        rn = n-tot_pts # remaining points to sample
        _pts = sample_random_points(3*rn)
        _Fx = interp_pts(gFx, _pts, method='linear')
        _d = np.abs(_Fx - μFx)/range_Fx  # distance from mean, normalized to [0, 1]
        # close_LCFS = (1 - np.random.uniform(0, 1, size=_d.shape)) > _d**(1/100) # importance sampling (?)
        # close_LCFS = np.random.uniform(0, 1, size=_d.shape)**2 > _d # importance sampling (?)
        close_LCFS = np.random.normal(0.0, σFx, size=_d.shape) > _d # importance sampling 
        inside_vessel = inpoly(_pts, VESS_CORNERS)
        _pts = _pts[close_LCFS & inside_vessel]
        _pts = _pts[:rn]  # keep only the remaining points
        pts[tot_pts:tot_pts+_pts.shape[0]] = _pts
        tot_pts += _pts.shape[0]
        iters += 1
        assert iters < 200, 'Too many iterations in importance sampling, try increasing σFx'
    return pts

def interp_pts(f:np.array, pts:np.array, gr=RRD[0,:], gz=ZZD[:,0], method=INTERP_METHOD):
    """
    Interpolate the function f at the points pts using the grid gr, gz.
    f should be a 2D array with shape (len(gr), len(gz)).
    pts should be a 2D array with shape (n, 2), where n is the number of points.
    """
    assert f.ndim == 2, f"f.ndim = {f.ndim}, f.shape = {f.shape}"
    assert pts.ndim == 2 and pts.shape[1] == 2, f"pts.ndim = {pts.ndim}, pts.shape = {pts.shape}"
    interp_func = RegularGridInterpolator((gr, gz), f.T, method=method)
    return interp_func(pts).reshape(-1)

def filter_zeros(r,z, thrs=1e-7):
    m = ~((r == 0) & (z == 0)) # exact zeros
    # m = ~((np.abs(r) < thrs) & (np.abs(z) < thrs)) # near zeros
    return r[m], z[m]

def get_in_masks(names:list):
    assert all(name in INPUT_NAMES for name in names), f"All names must be in INPUT_NAMES: {INPUT_NAMES}"
    mask = np.zeros(NIN, dtype=bool)
    if len(names) == 0: return mask
    mask[np.hstack([INPUT_MASKS[name] for name in names])] = True
    return mask

# calculate Br, Bz from the flux map (Fx)
def meqBrBz(f:np.ndarray, rr=RRD, zz=ZZD):
    if f.ndim == 2: f = f.reshape((1, *f.shape)) 
    assert f.ndim == 3, f'Input array must be 2D or 3D, got {f.ndim}D'
    # assert f.shape[1:] == (65, 28), f'Input array must have shape (N, 28, 65), got {f.shape}'

    dr, dz, r = rr[0,1]-rr[0,0], zz[1,0]-zz[0,0], rr[0]
    i4pirdr, i4pirdz = 1/(4*np.pi*dr*r), 1/(4*np.pi*dz*r)
    Br, Bz = np.zeros_like(f, dtype=DTYPE), np.zeros_like(f, dtype=DTYPE)

    # Br = -1/(2*pi*R)* dF/dz
    # Central differences for interior points: F/dz[i] =  F[i-1] - F[i+1]/(2*dz)
    Br[:,1:-1,:]    = -i4pirdz *  (f[:,2:,:] - f[:,:-2,:])
    # At grid boundary i, use: dF/dz[i] = (-F(i+2) + 4*F(i+1) - 3*F(i))/(2*dz)
    Br[:,-1,:]      = -i4pirdz * (+f[:,-3,:] - 4*f[:,-2,:] + 3*f[:,-1,:])
    Br[:,0,:]       = -i4pirdz * (-f[:,2,:]  + 4*f[:,1,:]  - 3*f[:,0,:])

    # Bz = 1/(2*pi*R)* dF/dr
    # Central differences dF/dz[i] =  F[i-1] - F[i+1]/(2*dz)
    Bz[:,:,1:-1]    = i4pirdr[1:-1] * (f[:,:,2:] - f[:,:,:-2])
    # At grid boundary i, use: dF/dz[i] = (-F(i+2) + 4*F(i+1) - 3*F(i))/(2*dz)
    Bz[:,:,-1]      = i4pirdr[-1] * (f[:,:,-3] - 4*f[:,:,-2] + 3*f[:,:,-1])
    Bz[:,:,0]       = i4pirdr[0] *  (-f[:,:,2] + 4*f[:,:,1]  - 3*f[:,:,0])
    return Br, Bz

def find_approx_lcfs(Iy):
    assert Iy.ndim == 2, f"Iy.ndim = {Iy.ndim}, should be 2D"
    assert Iy.shape == RRD.shape, f"Iy.shape = {Iy.shape}, RRD.shape = {RRD.shape}, should be equal"
    abs_Iy = np.abs(Iy)
    inIy = abs_Iy > 0
    Iy_border_mask = (
        (np.abs(np.diff(inIy, axis=0, append=0)) > 0) |
        (np.abs(np.diff(inIy, axis=1, append=0)) > 0)
    )
    return RRD[Iy_border_mask].reshape(-1), ZZD[Iy_border_mask].reshape(-1)

####################################################################################################
## PLOTTING
from matplotlib.patches import PathPatch
from matplotlib.path import Path
def _fill_between_polygons(ax, inner_poly, outer_poly, **kwargs):
    """
    Fill the area between two polygons.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to draw to
    inner_poly : array-like, shape (n, 2)
        Vertices of the inner polygon as (x, y) pairs
    outer_poly : array-like, shape (m, 2)
        Vertices of the outer polygon as (x, y) pairs
    **kwargs : dict
        Additional arguments passed to matplotlib.patches.PathPatch
        
    Returns:
    --------
    patch : matplotlib.patches.PathPatch
        The patch representing the filled area
    """
    # Convert to numpy arrays if they aren't already
    inner_poly = np.asarray(inner_poly)
    outer_poly = np.asarray(outer_poly)
    # Ensure polygons are closed
    if not np.array_equal(inner_poly[0], inner_poly[-1]):
        inner_poly = np.vstack([inner_poly, inner_poly[0]])
    if not np.array_equal(outer_poly[0], outer_poly[-1]):
        outer_poly = np.vstack([outer_poly, outer_poly[0]])
    # Create a path with two polygons
    # The inner polygon needs to be in reverse order to create a hole
    vertices = np.vstack([outer_poly, inner_poly[::-1]])
    # Create path codes
    n_outer = len(outer_poly)
    n_inner = len(inner_poly)
    codes = np.ones(n_outer + n_inner, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO  # Start of outer polygon
    codes[n_outer] = Path.MOVETO  # Start of inner polygon
    # Create the path and patch
    path = Path(vertices, codes)
    patch = PathPatch(path, **kwargs)
    # Add to the axes
    ax.add_patch(patch)
    return patch

def plot_vessel(ax=None, lw=1.5, alpha=0.3):
    if ax is None: ax = plt.gca()
    c = 'white' if DARK_THEME else 'black'
    ax.plot(VESS[:,0], VESS[:,1], color=c, lw=lw, alpha=alpha) # most inner
    ax.plot(VESSI[:,0], VESSI[:,1], color=c, lw=lw, alpha=alpha) # inner
    ax.plot(VESSO[:,0], VESSO[:,1], color=c, lw=lw, alpha=alpha) # outer
    _fill_between_polygons(ax, VESSI, VESSO, color='gray', alpha=alpha*0.8, lw=0)
    _fill_between_polygons(ax, VESS, VESSI, color='gray', alpha=alpha*0.3, lw=0)
    return ax

def test_plot_vessel():
    print("test_plot_vessel")
    plt.figure()
    plot_vessel()
    plt.axis('equal')
    plt.title('Vessel 2')
    # plt.show()
    plt.savefig(f"{TEST_DIR}/vessel.png")
    # plt.close()s

def prep_plot_tcv(ax):
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plot_vessel(ax) # plot vessel

def plot_network_outputs(ds:LiuqeDataset, model:liuqeML, title="test", save_dir=SAVE_DIR, nplt=2 if LOCAL else 50, use_db=True, use_diy=True):
    assert ds.fg is not None, "Dataset does not have full grid data (fg)"
    if use_diy: assert use_db, "use_diy requires use_db to be True"
    _ds_pps = ds.pps # tmp workaround for proper plotting
    ds.pps = ds.tnp # make ds to sample every point for every timestamp
    main_col = 'w' if DARK_THEME else 'k'
    if not use_db: model.eval()
    model.to(CPU)
    os.makedirs(f"{save_dir}/imgs", exist_ok=True)
    for pi, i in enumerate(np.random.randint(0, len(ds), nplt)):  
        on_grid = pi >= (nplt // 2) # on grid or not
        if not on_grid: continue # # TODO remove this
        fs = (20, 8) if on_grid else (15, 8) # figure size
        plt.figure(figsize=fs, facecolor='white' if DARK_THEME else 'black')
        x, p, fx, iy, br, bz = ds[i]
        x, p, fx, iy, br, bz = map(lambda x: x.to(CPU), [x, p, fx, iy, br, bz])
        if on_grid: 
            μ0 = 4e-7 * np.pi  # Permeability of free space
            p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)  # shape (n, 2)
            fx, iy, br, bz = [to_tensor(v.reshape(-1)) for v in ds.fg[i]] # reshape to (n,)
            iy = iy*μ0/GRA # workaround, on grid points dont have the conversion applied
        if use_db: r = p.clone()[:,0]; p.requires_grad = True # need r for PINN
        reverse_cmap = iy.mean() < 0  # reverse cmap if current density is negative
        fxp, brzp, iyp = model(x.unsqueeze(0), p.unsqueeze(0))
        if use_db:
            dfx_drz = torch.autograd.grad(outputs=fxp.sum(), inputs=p, create_graph=True)[0].unsqueeze(0) # compute gradients wrt p
            π = torch.pi
            dfx_dr, dfx_dz = dfx_drz[:,:,0], dfx_drz[:,:,1]  # shape (1, n)
            brp = -dfx_dz / (2*π*r) # Br from derivative
            bzp = dfx_dr / (2*π*r)  # Bz from derivative
            if use_diy: # calculate second derivative
                d2fx_dr2 = torch.autograd.grad(outputs=dfx_dr.sum(), inputs=p, create_graph=True)[0][:,0].unsqueeze(0)
                d2fx_dz2 = torch.autograd.grad(outputs=dfx_dz.sum(), inputs=p, create_graph=True)[0][:,1].unsqueeze(0)
                iyp = - (1/(2*π*r) * d2fx_dz2 + 1/(2*π*r) * d2fx_dr2 - 1/(2*π*r**2) * dfx_dr)

        else: brp, bzp = brzp[:,:,0], brzp[:,:,1]  # split Br, Bz
        if True:
            fxp, iyp, brp, bzp = map(lambda t: t.squeeze(0), [fxp, iyp, brp, bzp]) # remove batch dimension
            p, fx, iy, br, bz, fxp, iyp, brp, bzp = map(lambda t: t.detach().cpu().numpy(), [p, fx, iy, br, bz, fxp, iyp, brp, bzp]) # move to CPU and numpy
            iv = inpoly(p) # inside vessel indices
            # mfx, Mfx = np.min([fx[iv], fxp[iv]]), np.max([fx[iv], fxp[iv]]) # both maps
            # miy, Miy = np.min([iy[iv], iyp[iv]]), np.max([iy[iv], iyp[iv]])
            # mbr, Mbr = np.min([br[iv], brp[iv]]), np.max([br[iv], brp[iv]])
            # mbz, Mbz = np.min([bz[iv], bzp[iv]]), np.max([bz[iv], bzp[iv]])
            mfx, Mfx = np.min(fx[iv]), np.max(fx[iv]) # only actual
            miy, Miy = np.min(iy[iv]), np.max(iy[iv])
            mbr, Mbr = np.min(br[iv]), np.max(br[iv])
            mbz, Mbz = np.min(bz[iv]), np.max(bz[iv])
            # mfx, Mfx = np.min(fxp[iv]), np.max(fxp[iv]) # only predicted
            # miy, Miy = np.min(iyp[iv]), np.max(iyp[iv])
            # mbr, Mbr = np.min(brp[iv]), np.max(brp[iv])
            # mbz, Mbz = np.min(bzp[iv]), np.max(bzp[iv]) 
            lw, col, ms, cm = 1.5, 'gray', 4, plt.rcParams['image.cmap'] + ('_r' if reverse_cmap else '') # line width, color, marker size, colormap
            limits = [(mfx, Mfx), (miy, Miy), (mbr, Mbr), (mbz, Mbz)] # variable limits
            vals = [(fx, fxp), (iy, iyp), (br, brp), (bz, bzp)]
            names = ["Fx", "Iyd" if use_diy else "Iy", "Brd" if use_db else "Br", "Bzd" if use_db else "Bz"]
            nraws = 2
            ncols = 8 if on_grid else 6 # number of columns
            α_sep = 0.5
            for j, ((v, vp), l, n) in enumerate(zip(vals, limits, names)):
                plt.subplot(nraws, ncols, j*(ncols//2)+1), plt.title(f"{n} Actual")
                plt.scatter(p[iv,0], p[iv,1], c=v[iv], s=ms, vmin=l[0], vmax=l[1], cmap=cm)
                prep_plot_tcv(plt.gca()), plt.colorbar(); reverse_cmap and plt.gcf().axes and plt.gcf().axes[-1].invert_yaxis(); 
                plt.subplot(nraws, ncols, j*(ncols//2) + 2); plt.title(f"{n} Predicted")
                plt.scatter(p[iv,0], p[iv,1], c=vp[iv], s=ms, vmin=l[0], vmax=l[1], cmap=cm)
                prep_plot_tcv(plt.gca()), plt.colorbar(); reverse_cmap and plt.gcf().axes and plt.gcf().axes[-1].invert_yaxis();
                plt.subplot(nraws, ncols, j*(ncols//2) + ncols//2), plt.title(f"{n} Error %")
                ep = 100*np.abs(v[iv]-vp[iv])/(l[1]-l[0]) # error in %
                # plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=np.max(ep)) # variable bar for everything
                plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=np.max(ep) if n in [IY, IYD] else 1.0) # fixed bar but not for Iy
                prep_plot_tcv(plt.gca()), plt.colorbar()

                if on_grid: # plot contours if on grid
                    plt.subplot(nraws, ncols, j*(ncols//2) + 3), plt.title(f"{n} Contours")
                    plt.contour(RRD, ZZD, v.reshape(65,28), levels=20, linestyles='solid', cmap=cm)
                    plt.contour(RRD, ZZD, vp.reshape(65,28), levels=20, linestyles='dashed', cmap=cm)
                    prep_plot_tcv(plt.gca()), plt.colorbar()
            fig = plt.gcf() # Get the current figure
            dx, xoff = 0.004, 0.01 #0.005, 0.008
            # dy, yoff = dx*fs[0]/fs[1], -0.012 #dx*fs[0]/fs[1], 0.003
            dy, yoff = dx*fs[0]/fs[1], -0.018 #dx*fs[0]/fs[1], 0.003
            tvs = 0.041 # title verticale space
            from matplotlib.patches import Rectangle
            # for x0, y0, lx, ly in [(0.0+dx, 0.5+dy/2+yoff, (1-3*dx)/2+xoff, (1-3*dy)/2-yoff-tvs),
            #                        (0.5+dx/2+xoff, 0.5+dy/2+yoff, (1-3*dx)/2-xoff, (1-3*dy)/2-yoff-tvs),
            #                        (0.0+dx, 0.0+dy, (1-3*dx)/2+xoff, (1-3*dy)/2+yoff),
            #                        (0.5+dx/2+xoff, 0.0+dy, (1-3*dx)/2-xoff, (1-3*dy)/2+yoff)]:
            for x0, y0, lx, ly in [(0.0+dx,         0.5+dy/2+yoff, (1-3*dx)/2+xoff, (1-3*dy)/2-yoff-tvs),
                                (0.5+dx/2+xoff,  0.5+dy/2+yoff, (1-3*dx)/2-xoff, (1-3*dy)/2-yoff-tvs),
                                (0.0+dx,         0.0+dy, (1-3*dx)/2+xoff, (1-3*dy)/2+yoff),
                                (0.5+dx/2+xoff,  0.0+dy, (1-3*dx)/2-xoff, (1-3*dy)/2+yoff)]:
                fig.patches.append(Rectangle( (x0, y0), lx, ly, transform=fig.transFigure, color='black' if DARK_THEME else 'white', alpha=1, zorder=-1))
            plt.suptitle(f"[{JOBID}] liuqeML: {title} {i}", y=0.985, color='black' if DARK_THEME else 'white')
            plt.tight_layout()
            plt.savefig(f"{save_dir}/imgs/net_example_{'g_' if on_grid else ''}{title}_{i}.png")
            if LOCAL: plt.show()  
            plt.close()

    ds.pps = _ds_pps # restore dataset pps
    return


if __name__ == '__main__':
    test_network_io()
    test_plot_vessel()
    if LOCAL: plt.show()

    