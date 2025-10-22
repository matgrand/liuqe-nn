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
plt.rcParams['grid.alpha'] = 0.3
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

PINN = True # if True, use PINN, otherwise use standard model
NDOM = 4 # number of plasma domains, should be matched with create_ds.m
##################################################################################################################################
## Dataset inputs/outputs (see create_ds.m)
# input names
BM = 'Bm' # magnetic probe measurements
FF = 'Ff' # flux loops poloidal flux 
IA = 'Ia' # poloidal field coil currents
IP = 'Ip' # plasma current
# output names
FX = 'Fx' # Plasma poloidal flux map
IY = 'Iy' # Plasma current density map
BR = 'Br' # recalculated in python, not in matlab
BZ = 'Bz' # recalculated in python, not in matlab
# RQ = 'rq' # LCFS r coordinates
# ZQ = 'zq' # LCFS z coordinates
#more names
# SEP = 'sep' # this is the LCFS, last closed flux surface, or separatrix
PHYS = 'phys' # physics inputs
PTS = 'pts' # points 
FG = 'fg' # full grid
RT = 'rt' # real time output
BRP = 'Brp' # Br from PINN  
BZP = 'Bzp' # Bz from PINN

INPUT_NAMES = [BM, FF, IA, IP] # input names
# OUTPUT_NAMES = [FX, IY, BR, BZ, RQ, ZQ] # output names
OUTPUT_NAMES = [FX, IY, BR, BZ] # output names
DS_NAMES = INPUT_NAMES + OUTPUT_NAMES # dataset names
DS_SIZES = { 
    BM:(38,), FF:(38,), IA:(19,), IP:(1,),  # input sizes
    # FX:(65,28), IY:(63,26), BR:(65,28), BZ:(65,28), RQ:(32*NDOM,), ZQ:(32*NDOM,),  # output sizes
    FX:(65,28), IY:(63,26), BR:(65,28), BZ:(65,28),  # output sizes
    # SEP:(2*129,), # LCFS/SEP size, 2*NLCFS points (r,z) for the last closed flux surface
}
DTYPE = 'float32'

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
# LOSS_NAMES = [FX, IY, BR, BZ, SEP, BRP, BZP] # loss names
LOSS_NAMES = [FX, IY, BR, BZ, BRP, BZP] # loss names
def model_path(loss_name, save_dir=SAVE_DIR):
    assert loss_name in LOSS_NAMES, f"loss_name should be one of {LOSS_NAMES}, got {loss_name}"
    return f"{save_dir}/best_{loss_name}.pth"

BEST_MODEL_DIR = 'data/best' # where the best models are saved
STRICT_LOAD = True # for loading the weights, should be true, but for testing varying architectures, set to false

TEST_DIR = 'test' if LOCAL else '/nfsd/automatica/grandinmat/test'
os.makedirs(TEST_DIR, exist_ok=True)

NIN = sum(DS_SIZES[name][0] for name in INPUT_NAMES) # 96

NLCFS = 32*NDOM # number of LCFS points 

# assert NLCFS == DS_SIZES[RQ][0] == DS_SIZES[ZQ][0], f"NLCFS = {NLCFS}, DS_SIZES[RQ] = {DS_SIZES[RQ]}, DS_SIZES[ZQ] = {DS_SIZES[ZQ]}"

##################################################################################################################################


# read the original grid coordinates
d = loadmat('tcv_params/grid.mat')
rd, zd = d['r'].flatten(), d['z'].flatten() # original grid coordinates (DATA)
RRD, ZZD = np.meshgrid(rd, zd)  # meshgrid for the original grid coordinates (from the data)
RRD2, ZZD2 = np.meshgrid(rd[1:-1], zd[1:-1])  
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
# class ActF(Module): # swish
#     def __init__(self): 
#         super(ActF, self).__init__()
#         self.beta = torch.nn.Parameter(torch.tensor(1.0), requires_grad=True)
#     def forward(self, x): return x*torch.sigmoid(self.beta*x)

class ActF(Module): # swish
    def __init__(self, no=1): 
        super(ActF, self).__init__()
        if no == 1: self.beta = torch.nn.Parameter(torch.tensor(1.0), requires_grad=True)
        else: self.beta = torch.nn.Parameter(torch.ones(no), requires_grad=True) # n is the number of outputs
    def forward(self, x): return x*torch.sigmoid(self.beta*x)

# class ActF(Module): # Beta Relu
#     def __init__(self, no=1): 
#         super(ActF, self).__init__()
#         if no == 1: self.beta = torch.nn.Parameter(torch.tensor(1.0), requires_grad=True)
#         else: self.beta = torch.nn.Parameter(torch.ones(no), requires_grad=True) # n is the number of outputs
#     def forward(self, x): return self.beta * F.relu(x)

def percentage_loss(a, b): # percentage loss (a = true, b = predicted)
    assert a.shape == b.shape, f"a.shape = {a.shape}, b.shape = {b.shape}, should be equal"
    assert a.dim() >= 2, f"a.dim() = {a.dim()}, should be >= 2"
    span = torch.max(a, dim=-1).values - torch.min(a, dim=-1).values 
    assert torch.all(span > 1e-10), f"span should be positive, but got {span.min():.3e}"
    err = torch.abs(a - b) / span.unsqueeze(-1)  # calculate the error as a percentage of the span
    return torch.mean(err)  # return the mean error for each batch 

# PHYSICS_LS = 64 # 128 physics latent size [ph] 64 <- 64 -> 40μ
# PHYSICS_LS = 128 #128 -> 50μ
# PHYSICS_LS = 256 # -> 70μ
# PHYSICS_LS = 16 # turbosmall

# PHYSICS_LS, PTSH, INPUTH, HEADH = 64, 64, 128, 128 # medium (old baseline)
PHYSICS_LS, PTSH, INPUTH, HEADH = 32, 32, 48, 48 # small (new baseline) <- 
# PHYSICS_LS, PTSH, INPUTH, HEADH = 16, 16, 32, 32 # turbosmall 
# PHYSICS_LS, PTSH, INPUTH, HEADH = 3, 16, 32, 32 # for shot mapping


class PtsEncoder(Module): # positional encoding for the input vector
    def __init__(self, ni=2, nh=PTSH, no=PHYSICS_LS):
        super(PtsEncoder, self).__init__()
        self.pts_enc = Sequential(
            Linear(ni, nh), ActF(nh),
            Linear(nh, no), ActF(no),
        )
    def forward(self, pts): return self.pts_enc(pts) # pts: (BS, NP, 2) -> (BS, NP, PHYSICS_LS)

class InputNet(Module): # input -> latent physics vector [x -> ph]
    def __init__(self, x_mean_std=torch.tensor([[0.0]*NIN, [1.0]*NIN], dtype=torch.float32),
                 ni=NIN, nh=INPUTH, no=PHYSICS_LS):
        super(InputNet, self).__init__()
        assert x_mean_std.shape == (2, NIN), f"x_mean_std.shape = {x_mean_std.shape}, NIN = {NIN}"
        self.x_mean_std = x_mean_std
        self.phys_enc = Sequential(
            Linear(ni, nh), ActF(nh), 
            Linear(nh, no), ActF(no), 
        )
    def forward(self, x): 
        assert x.shape[1] == NIN, f"x.shape[1] = {x.shape[1]}, NIN = {NIN}"
        x = (x - self.x_mean_std[0]) / self.x_mean_std[1] # normalize
        ph = self.phys_enc(x)
        return ph
    def to(self, device):
        super(InputNet, self).to(device)
        self.x_mean_std = self.x_mean_std.to(device)
        return self

class FHead(Module): # [pt, ph] -> [1] function (flux/Br/Bz/curr density) 
    def __init__(self, ni=PHYSICS_LS, nh=HEADH, no=1):
        super(FHead, self).__init__()
        self.head = Sequential(
            Linear(ni, nh), ActF(nh), # <-
            Linear(nh, no), ActF(no),
        )
        self.no = no # number of outputs
    def forward(self, v): return self.head(v) if self.no > 1 else self.head(v).squeeze(-1)

class LCFSHead(Module): # physics -> LCFS [ph -> LCFS]
    def __init__(self):
        super(LCFSHead, self).__init__()
        self.lcfs = Sequential(
            # Linear(PHYSICS_LS, 64), ActF(),
            # Linear(64, 32), ActF(),
            # Linear(32, NLCFS*2), ActF(),

            # Linear(PHYSICS_LS, 128), ActF(),
            # Linear(128, 64), ActF(),
            # Linear(64, NLCFS*2), ActF(),  # output is (NLCFS*2) for r and z coordinates

            # Linear(PHYSICS_LS, 64), ActF(), # 1 layer
            # Linear(64, NLCFS*2), ActF(),
    
            Linear(PHYSICS_LS, 16), ActF(), # 1 layer
            Linear(16, NLCFS*2), ActF(),

            # Linear(PHYSICS_LS, 32), ActF(), # 1 layer
            # Linear(32, NLCFS*2), ActF(),

            # Linear(PHYSICS_LS, 2), ActF(2), # unused
            # Linear(2, NLCFS*2), ActF(NLCFS*2), 
        )
    def forward(self, ph): return self.lcfs(ph)

def aggregate_pts_ph(ps, ph): # ps points, ph physics vector
    # simple multiplication
    assert ps.shape == ph.shape, f"ps.shape = {ps.shape}, ph.shape = {ph.shape}, should be equal" 
    psh =  ps * ph
    assert psh.dim() == 3 and psh.shape[2] == PHYSICS_LS, f"psh.dim() = {psh.dim()}, psh.shape[2] = {psh.shape[2]}, should be PHYSICS_LS = {PHYSICS_LS}"
    return psh # multiply the points with the physics vector

# def aggregate_pts_ph(ps, ph): # ps points, ph physics vector
#     # multiply only half of the physics vector with the points
#     assert ps.ndim == ph.ndim, f"ps.ndim = {ps.ndim}, ph.ndim = {ph.ndim}, should be equal"
#     assert ps.shape[:-1] == ph.shape[:-1], f"ps.shape = {ps.shape}, ph.shape = {ph.shape}, should be equal except last dimension"
#     assert ph.shape[-1] == PHYSICS_LS, f"ph.shape = {ph.shape}, PHYSICS_LS = {PHYSICS_LS}"
#     assert ps.shape[-1] == PHYSICS_LS//2, f"ps.shape = {ps.shape}, PHYSICS_LS//2 = {PHYSICS_LS//2}"
#     psh = ps * ph[..., :PHYSICS_LS//2] # multiply only the first half of the physics vector with the points
#     assert psh.dim() == 3 and psh.shape[2] == PHYSICS_LS//2, f"psh.dim() = {psh.dim()}, psh.shape[2] = {psh.shape[2]}, should be PHYSICS_LS = {PHYSICS_LS//2}"
#     return psh # multiply the points with the physics vector, return (BS, NP, PHYSICS_LS//2)

class FullNet(Module): # [pt, ph] -> [1] function (flux/Br/Bz/curr density)
    def __init__(self, input_net:InputNet, pts_enc:PtsEncoder, rt_head:FHead, iy_head:FHead, lcfs_head:LCFSHead):
        super(FullNet, self).__init__()
        self.input_net = input_net
        self.pts_enc = pts_enc
        self.iy_head = iy_head
        self.rt_head = rt_head
        self.lcfs_head = lcfs_head
        assert isinstance(self.input_net, InputNet), f"input_net should be an instance of InputNet, got {type(self.input_net)}"
        assert isinstance(self.pts_enc, PtsEncoder), f"pts_enc should be an instance of PtsEncoder, got {type(self.pts_enc)}"
        assert isinstance(self.rt_head, FHead), f"rt_head should be an instance of FHead, got {type(self.rt_head)}"
        assert isinstance(self.iy_head, FHead), f"iy_head should be an instance of FHead, got {type(self.iy_head)}"
    def to(self, device):
        super(FullNet, self).to(device)
        self.input_net.to(device)
        return self
    def forward(self, x, pts):
        assert pts.shape[-1] == 2, f"pts.shape[-1] = {pts.shape[-1]}, should be 2 (r,z)"
        assert pts.dim() == 3, f"pts.dim() = {pts.dim()}, should be 3 (batch, n_points, 2)"
        assert x.dim() == 2, f"x.dim() = {x.dim()}, should be 2 (batch, NIN)"
        assert x.shape[-1] == NIN, f"x.shape[1] = {x.shape[1]}, NIN = {NIN}"
        assert pts.shape[0] == x.shape[0], f"pts.shape[0] = {pts.shape[0]}, x.shape[0] = {x.shape[0]}, should be equal (batch size)"
        # calculate physics vector ph and repeat it for each point in pts
        ph = self.input_net(x) # get physics vector 
        lcfs = self.lcfs_head(ph) # get LCFS (before expanding ph)
        ph = ph.unsqueeze(1).expand(-1, pts.shape[1], -1) # (BS, NP, PHYSICS_LS)
        # calculate a positional encoding for the points pts
        ps = self.pts_enc(pts) # encode points       
        psh = aggregate_pts_ph(ps, ph) # aggregate points and physics vector
        iy = self.iy_head(psh) # get current density
        rt = self.rt_head(psh) # get flux, Br, Bz (rt quantities)
        assert iy.dim() == 2 and iy.shape[1] == pts.shape[1], f"iy.dim() = {iy.dim()}, iy.shape[1] = {iy.shape[1]}, should be pts.shape[1] = {pts.shape[1]}"
        assert rt.dim() == 3 and rt.shape[1] == pts.shape[1] and rt.shape[2] == 3, f"rt.dim() = {rt.dim()}, rt.shape = {rt.shape}, should be (BS, NP, 4)"
        assert lcfs.dim() == 2 and lcfs.shape[1] == NLCFS*2, f"lcfs.dim() = {lcfs.dim()}, lcfs.shape[1] = {lcfs.shape[1]}, should be 2*NLCFS = {NLCFS*2}"
        return rt, iy, lcfs # return rt (flux, Br, Bz), current density iy, and LCFS

class LiuqeRTNet(Module): # Liuqe Real Time surrogate net
    def __init__(self, input_net:InputNet, pts_enc:PtsEncoder, rt_head:FHead):
        super(LiuqeRTNet, self).__init__()
        self.input_net = input_net
        self.pts_enc = pts_enc
        self.rt_head = rt_head
        assert isinstance(self.input_net, InputNet), f"input_net should be an instance of InputNet, got {type(self.input_net)}"
        assert isinstance(self.pts_enc, PtsEncoder), f"pts_enc should be an instance of PtsEncoder, got {type(self.pts_enc)}"
        assert isinstance(self.rt_head, FHead), f"rt_head should be an instance of FHead, got {type(self.rt_head)}"
    def to(self, device):
        super(LiuqeRTNet, self).to(device)
        self.input_net.to(device)
        return self
    def forward(self, x, r, z):
        assert x.dim() == 1, f"x.dim() = {x.dim()}, should be 1(NIN,)"
        assert x.shape == (NIN,), f"x.shape = {x.shape}, should be (NIN,)"
        assert r.dim() == 1 and z.dim() == 1, f"r.dim() = {r.dim()}, z.dim() = {z.dim()}, should be 1"
        assert r.shape == z.shape, f"r.shape = {r.shape}, z.shape = {z.shape}, should be equal"
        n = r.shape[0] # number of points
        pts = torch.stack([r, z], dim=1).view(1,n,2) # create points tensor (n, 2)
        ph = self.input_net(x.view(1, -1)) # get physics vector
        ph = ph.unsqueeze(1).expand(-1, n, -1)
        ps = self.pts_enc(pts) # encode points
        psh = aggregate_pts_ph(ps, ph) # aggregate points and physics vector
        rt = self.rt_head(psh)
        return rt.view(n, 3) # return flux, Br, Bz as (n, 3) tensor

class LCFSNet(Module):
    def __init__(self, input_net:InputNet, lcfs_head:LCFSHead):
        super(LCFSNet, self).__init__()
        self.input_net = input_net
        self.lcfs_head = lcfs_head
        assert isinstance(self.input_net, InputNet), f"input_net should be an instance of InputNet, got {type(self.input_net)}"
        assert isinstance(self.lcfs_head, LCFSHead), f"lcfs_head should be an instance of LCFSHead, got {type(self.lcfs_head)}"
    def to(self, device):
        super(LCFSNet, self).to(device)
        self.input_net.to(device)
        return self
    def forward(self, x):
        ph = self.input_net(x) # get physics vector
        assert ph.dim() == 2 and ph.shape[1] == PHYSICS_LS, f"ph.dim() = {ph.dim()}, ph.shape[1] = {ph.shape[1]}, PHYSICS_LS = {PHYSICS_LS}"
        lcfs = self.lcfs_head(ph) # get LCFS
        assert lcfs.dim() == 2 and lcfs.shape[1] == NLCFS*2, f"lcfs.dim() = {lcfs.dim()}, lcfs.shape[1] = {lcfs.shape[1]}, should be 2*NLCFS = {NLCFS*2}"
        return lcfs
    
def test_network_io(verbose=True):
    v = verbose
    if v: print('test_network_io')
    n_points = 5 # number of points to sample
    # single sample
    x, pts = (torch.rand(1, NIN), torch.rand(1, n_points, 2)) # x: (1, NIN), pts: (1, n_points, 2)
    fn = FullNet(InputNet(),PtsEncoder(),FHead(no=3),FHead(no=1),LCFSHead())
    rt_net = LiuqeRTNet(fn.input_net,fn.pts_enc,fn.rt_head)
    rt, iy, lcfs = fn(x, pts)
    r, z = pts[0,:,0], pts[0,:,1] # r, z: (n_points,)
    rt2 = rt_net(x.view(-1), r, z) # fx2, br2, bz2: (1, n_points, 3)
    assert rt.shape == (1, n_points, 3), f"rt.shape = {rt.shape}, should be (1, {n_points}, 3)"
    assert rt2.shape == (n_points, 3), f"rt2.shape = {rt2.shape}, should be ({n_points}, 3)"
    assert pts.shape == (1, n_points, 2), f"pts.shape = {pts.shape}, should be (1, {n_points}, 2)"
    assert lcfs.shape == (1, NLCFS*2), f"lcfs.shape = {lcfs.shape}, should be (1, {NLCFS*2})"
    assert torch.allclose(rt, rt2), f"rt and rt2 should be close, but got {torch.max(torch.abs(rt-rt2))}"
    if v: print(f"single FullNet -> in: [{x.shape}, {pts.shape}], \n            out: [{rt.shape}, {iy.shape}, {lcfs.shape}]")
    if v: print(f"single LiuqeRTNet  -> in: [{x.shape}, {pts.shape}], \n            out: [{rt2.shape}]") 
    # batched
    bs = 7
    x, pts = torch.rand(bs, NIN), torch.rand(bs, n_points, 2) # x: (bs, NIN), pts: (bs, n_points, 2)
    rt, iy, lcfs = fn(x, pts)
    assert rt.shape == (bs, n_points, 3), f"rt.shape = {rt.shape}, should be ({bs}, {n_points}, 3)"
    assert lcfs.shape == (bs, NLCFS*2), f"lcfs.shape = {lcfs.shape}, should be ({bs}, {NLCFS*2})"
    if v: print(f"batched FullNet -> in: [{x.shape}, {pts.shape}], \n            out: [{rt.shape}, {iy.shape}, {lcfs.shape}]")

def convert_to_onnx_dyn(net:LiuqeRTNet, save_dir=[SAVE_DIR]):
    net.to(CPU)  
    net.eval()  # Set the network to evaluation mode
    dummy_phys = torch.randn(NIN, device=CPU)
    n = 24 # number of inference points
    dummy_r = torch.randn(n, device=CPU)  # Dummy points for inference
    dummy_z = torch.randn(n, device=CPU)  # Dummy points for inference
    for d in save_dir:
        onnx_net_path1 = f'{d}/net_dyn.onnx'
        onnx_net_path2 = f'{d}/net_dyn_ops17.onnx'
        # Export to ONNX
        torch.onnx.export(
            net,
            args=(dummy_phys, dummy_r, dummy_z),  # Dummy inputs for the network
            f=onnx_net_path1,
            input_names=[PHYS, 'r', 'z'],  # Input names
            output_names=[RT],
            dynamic_axes={"r": {0: "n"}, "z": {0: "n"}, RT: {0: "n"}},
            opset_version=15,  # 15 for onnx 1.11, onnx 1.11 bc GNU 5.3.1
        )
        torch.onnx.export(
            net,
            args=(dummy_phys, dummy_r, dummy_z),  # Dummy inputs for the network
            f=onnx_net_path2,
            input_names=[PHYS, 'r', 'z'],  # Input names
            output_names=[RT],
            dynamic_axes={"r": {0: "n"}, "z": {0: "n"}, RT: {0: "n"}},
            opset_version=17, 
        )
        assert os.path.exists(onnx_net_path1), f"ONNX model not saved to {onnx_net_path1}"
        assert os.path.exists(onnx_net_path2), f"ONNX model not saved to {onnx_net_path2}"
        print(f'ONNX models saved to {onnx_net_path1} and {onnx_net_path2}')
    return

def convert_to_onnx_static(net:LiuqeRTNet, npts=N_CTRL_PTS, save_dir=[SAVE_DIR]):
    net.to(CPU)  
    net.eval()  # Set the network to evaluation mode
    dummy_phys = torch.randn(NIN, device=CPU)
    dummy_r = torch.randn(npts, device=CPU)  # Dummy points for inference
    dummy_z = torch.randn(npts, device=CPU)  # Dummy points for inference
    for d in save_dir:
        onnx_net_path = f'{d}/net_static.onnx'
        # Export to ONNX
        onnx_model = torch.onnx.export(
            net, 
            (dummy_phys, dummy_r, dummy_z), 
            dynamo=True,
            input_names=[PHYS, 'r', 'z'],  # Input names
            output_names=[RT],
            opset_version=15,  # 15 for onnx 1.11, onnx 1.11 bc GNU 5.3.1
        )
        # onnx_model.optimize()
        onnx_model.save(onnx_net_path)
        assert os.path.exists(onnx_net_path), f"ONNX model not saved to {onnx_net_path}"
        print(f'ONNX model saved to {onnx_net_path}')
    return

def save_weights_to_mat(net:LiuqeRTNet, save_dir=[SAVE_DIR]):
    weights = {}
    weights['x_mean_std'] = net.input_net.x_mean_std.cpu().numpy()  # save normalization parameters
    for net in [net.input_net, net.pts_enc, net.rt_head]:
        for name, param in net.named_parameters():
            w = param.data.cpu().numpy()
            n = name.replace('.', '_')
            # print(f"Parameter: {n}: {w.shape}, {w.dtype}")
            weights[n] = w
    for sd in save_dir:
        savemat(f'{sd}/weights.mat', weights)
        print(f"Weights saved to {sd}/weights.mat")

# function to load the dataset
def load_ds(ds_path):
    assert os.path.exists(ds_path), f"Dataset not found: {ds_path}"
    d = np.load(ds_path)
    r = {
        PHYS:d[PHYS], 
        PTS:d[PTS], 
        FX:d[FX], 
        IY:d[IY], 
        BR:d[BR], 
        BZ:d[BZ], 
        # SEP:d[SEP], 
        "x_mean_std":d["x_mean_std"],
        FG:d[FG] if FG in d else None, # full grid, optional
    }
    return r

####################################################################################################
class LiuqeDataset(Dataset):
    def __init__(self, ds_path, verbose=True):
        d = load_ds(ds_path)
        x_mean_std = to_tensor(d.pop("x_mean_std")) # mean and std for inputs
        self.fg = d.pop(FG) # pop the full grid 
        print(f"Dataset loaded from {ds_path}, keys: {d.keys()}")
        self.data = {k:to_tensor(v) for k,v in d.items()} # convert to tensors
        # move to DEV (doable bc the dataset is fairly small, check memory usage)
        tot_memory_ds = sum([x.element_size()*x.nelement() for x in self.data.values()]) # total memory in bytes
        gpu_free_mem = torch.cuda.mem_get_info()[0] if DEV == CUDA else np.inf
        self.on_dev = DEV != CPU and tot_memory_ds < ( gpu_free_mem - 2e9 ) # 2 gb buffer for batches
        if self.on_dev: self.data = {k:v.to(DEV) for k,v in self.data.items()} # move to device
        self.x_mean_std = x_mean_std.to(DEV) if self.on_dev else x_mean_std
        if verbose: print(f"Dataset: N:{len(self)}, memory:{tot_memory_ds/1e6}MB, on_dev:{self.on_dev}")
    def __len__(self): return len(self.data[PHYS]) # number of samples, all data should have the same length
    def __getitem__(self, idx): return [x[idx] for x in self.data.values()] # [PHYS, PTS, FX, IY, BR, BZ]

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
        # x, pts, fx, iy, br, bz, sep = map(lambda x: x.cpu().numpy(), ds[j])
        x, pts, fx, iy, br, bz = map(lambda x: x.cpu().numpy(), ds[j])
        μ, σ = ds.x_mean_std.cpu().numpy()
        x = (x-μ) / σ  # normalize inputs
        r,z = pts[:,0], pts[:,1]
        # rq, zq = filter_zeros(sep[:NLCFS], sep[NLCFS:]) # filter zeros from LCFS
        ms = 1
        plt.subplot(n_plot, 5, i*5+1)
        plt.scatter(r, z, c=fx, s=ms), plt.title('FX')
        # plt.plot(rq, zq, 'gray', lw=2)
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+2)
        plt.scatter(r, z, c=iy, s=ms), plt.title('IY')
        # plt.plot(rq, zq, 'gray', lw=2)
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+3)
        plt.scatter(r, z, c=br, s=ms), plt.title('Br')
        # plt.plot(rq, zq, 'gray', lw=2)
        plot_vessel(), plt.axis('equal'), plt.axis('off'), plt.colorbar()
        plt.subplot(n_plot, 5, i*5+4)
        plt.scatter(r, z, c=bz, s=ms), plt.title('Bz')
        # plt.plot(rq, zq, 'gray', lw=2)
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

def plot_network_outputs(ds:LiuqeDataset, model:FullNet, title="test", save_dir=SAVE_DIR, nplt=2 if LOCAL else 50, show_lcfs=False):
    assert ds.fg is not None, "Dataset does not have full grid data (fg)"
    main_col = 'w' if DARK_THEME else 'k'
    model.eval()
    model.to(CPU)
    os.makedirs(f"{save_dir}/imgs", exist_ok=True)
    for pi, i in enumerate(np.random.randint(0, len(ds), nplt)):  
        on_grid = pi >= (nplt // 2) # on grid or not
        fs = (20, 8) if on_grid else (15, 8) # figure size
        plt.figure(figsize=fs, facecolor='white' if DARK_THEME else 'black')
        x, p, fx, iy, br, bz = ds[i]
        x, p, fx, iy, br, bz = map(lambda x: x.to(CPU), [x, p, fx, iy, br, bz])
        if on_grid: 
            p = torch.stack([torch.tensor(RRD.reshape(-1), dtype=torch.float32), torch.tensor(ZZD.reshape(-1), dtype=torch.float32)], dim=1)  # shape (n, 2)
            fx, iy, br, bz = [to_tensor(v.reshape(-1)) for v in ds.fg[i]] # reshape to (n,)
        rtp, iyp, _ = model(x.unsqueeze(0), p.unsqueeze(0))
        fxp, brp, bzp = rtp[:,:,0], rtp[:,:,1], rtp[:,:,2]
        fxp, iyp, brp, bzp = map(lambda t: t.squeeze(0), [fxp, iyp, brp, bzp]) # remove batch dimension
        p, fx, iy, br, bz, fxp, iyp, brp, bzp = map(lambda t: t.detach().cpu().numpy(), [p, fx, iy, br, bz, fxp, iyp, brp, bzp]) # move to CPU and numpy
        iv = inpoly(p) # inside vessel indices
        mfx, Mfx = np.min([fx[iv], fxp[iv]]), np.max([fx[iv], fxp[iv]])
        miy, Miy = np.min([iy[iv], iyp[iv]]), np.max([iy[iv], iyp[iv]])
        mbr, Mbr = np.min([br[iv], brp[iv]]), np.max([br[iv], brp[iv]])
        mbz, Mbz = np.min([bz[iv], bzp[iv]]), np.max([bz[iv], bzp[iv]])
        lw, col, ms = 1.5, 'gray', 4
        limits = [(mfx, Mfx), (miy, Miy), (mbr, Mbr), (mbz, Mbz)] # variable limits
        vals = [(fx, fxp), (iy, iyp), (br, brp), (bz, bzp)]
        names = ["Fx", "Iy", "Br", "Bz"]
        def prep_plot(ax):
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            plot_vessel(ax) # plot vessel
        nraws = 2
        ncols = 8 if on_grid else 6 # number of columns
        α_sep = 0.5
        for j, ((v, vp), l, n) in enumerate(zip(vals, limits, names)):
            plt.subplot(nraws, ncols, j*(ncols//2)+1), plt.title(f"{n} Actual")
            plt.scatter(p[iv,0], p[iv,1], c=v[iv], s=ms, vmin=l[0], vmax=l[1])
            # if show_lcfs: plt.plot(sep[:NLCFS], sep[NLCFS:], f'--{main_col}', lw=1, alpha=α_sep) # plot LCFS
            prep_plot(plt.gca()), plt.colorbar()
            plt.subplot(nraws, ncols, j*(ncols//2) + 2), plt.title(f"{n} Predicted")
            plt.scatter(p[iv,0], p[iv,1], c=vp[iv], s=ms, vmin=l[0], vmax=l[1])
            # if show_lcfs: plt.plot(sepp[:NLCFS], sepp[NLCFS:], f'--{main_col}', lw=1, alpha=α_sep)
            prep_plot(plt.gca()), plt.colorbar()
            plt.subplot(nraws, ncols, j*(ncols//2) + ncols//2), plt.title(f"{n} Error %")
            ep = 100*np.abs(v[iv]-vp[iv])/(l[1]-l[0]) # error in %
            # plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=np.max(ep)) # variable bar
            plt.scatter(p[iv,0], p[iv,1], c=ep, s=ms, vmin=0, vmax=np.max(ep) if n == IY else 1.0) # fixed bar
            # if show_lcfs: plt.plot(sep[:NLCFS], sep[NLCFS:], f'--{main_col}', lw=1, alpha=α_sep)
            prep_plot(plt.gca()), plt.colorbar()
            if on_grid: # plot contours if on grid
                plt.subplot(nraws, ncols, j*(ncols//2) + 3), plt.title(f"{n} Contours")
                plt.contour(RRD, ZZD, v.reshape(65,28), levels=20, linestyles='solid')
                plt.contour(RRD, ZZD, vp.reshape(65,28), levels=20, linestyles='dashed')
                prep_plot(plt.gca()), plt.colorbar()
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
        plt.suptitle(f"[{JOBID}] FullNet: {title} {i}", y=0.985, color='black' if DARK_THEME else 'white')
        plt.tight_layout()
        plt.savefig(f"{save_dir}/imgs/net_example_{'g_' if on_grid else ''}{title}_{i}.png")
        if LOCAL: plt.show()  
        plt.close()
    return
    
def plot_lcfs_net_out(ds:LiuqeDataset, model:LiuqeRTNet, title='test', save_dir=SAVE_DIR, nplt=5 if LOCAL else 50):
    main_col = 'w' if DARK_THEME else 'k'
    model.eval()
    model.to(CPU)
    os.makedirs(f"{save_dir}/imgs", exist_ok=True)
    lw3 = 1.5
    for i in np.random.randint(0, len(ds), nplt):  
        plt.figure(figsize=(20, 8))
        x, y3 = ds[i][0].to(CPU), ds[i][6].to(CPU)  # x: (NIN,), pts: (n_points, 2), y3: (2*NLCFS,)
        x = x.reshape(1, -1)
        yp3 = model(x)
        yp3 = yp3.detach().numpy().reshape(2 * NLCFS)
        y3 = y3.detach().numpy().reshape(2 * NLCFS)
        # convert to 2d points
        y3 = np.array([y3[:NLCFS], y3[NLCFS:]]).T
        yp3 = np.array([yp3[:NLCFS], yp3[NLCFS:]]).T
        err = np.linalg.norm(y3 - yp3, axis=-1)
        norm_err = err / np.max(err)
        # err_colors = np.array([plt.cm.viridis(norm_err[j]) for j in range(NLCFS)]) #viridis
        err_colors = np.array([plt.cm.plasma(norm_err[j]) for j in range(NLCFS)]) # plasma

        plt.subplot(1, 3, 1)
        plt.plot(y3[:,0], y3[:,1], lw=lw3, label='actual')
        plt.plot(yp3[:,0], yp3[:,1], lw=lw3, label='predicted')
        plot_vessel()
        plt.title("LCFS")
        plt.axis('equal')
        plt.xlabel("R [m]")
        plt.ylabel("Z [m]")
        plt.legend(loc='upper right')

        plt.subplot(1, 3, 2)
        α = 0.8
        plt.plot(y3[:,0], y3[:,1], lw=1, color=main_col, alpha=α)
        plt.plot(yp3[:,0], yp3[:,1], lw=1, color=main_col, alpha=α)
        # sc = plt.scatter(yp3[:,0], yp3[:,1], c=err, s=20, cmap='viridis')
        sc = plt.scatter(yp3[:,0], yp3[:,1], c=err, s=20, cmap='plasma')
        plt.colorbar(sc)
        plt.title("LCFS Absolute Error [m]")
        plt.axis('equal')
        plt.xlabel("R [m]")
        plt.ylabel("Z [m]")

        plt.subplot(1, 3, 3)
        plt.plot(y3[:,0], y3[:,1], lw=1, color=main_col, alpha=α)
        plt.plot(yp3[:,0], yp3[:,1], lw=1, color=main_col, alpha=α)
        for j in range(NLCFS):
            plt.plot([y3[j,0], yp3[j,0]], [y3[j,1], yp3[j,1]], color=err_colors[j], lw=2)
        plt.title("LCFS Absolute Error [m]")
        plt.axis('equal')
        plt.xlabel("R [m]")
        plt.ylabel("Z [m]")

        plt.suptitle(f"[{JOBID}] LCFSNet: {title} {i}")
        plt.tight_layout()
        plt.savefig(f"{save_dir}/imgs/lcfs_example_{title}_{i}.png")
        if LOCAL: plt.show()
        plt.close()
    return

if __name__ == '__main__':
    test_network_io()
    test_plot_vessel()
    if LOCAL: plt.show()

    