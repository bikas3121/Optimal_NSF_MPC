
# %%
import numpy as np
from scipy  import signal
import scipy
import csv
import matplotlib.pyplot as plt
import statistics
import math 
import gurobipy as gp
from gurobipy import GRB
import tqdm


from DirectQantization import quantise_signal, generate_code, generate_dac_output
from  quantiser_configurations import quantiser_configurations, get_measured_levels

from MHOQ import MHOQ
from MHOQ_BIN import MHOQ_BIN
from nsdcal import nsdcal, noise_shaping
from dsm import dsm


from signalProcessing import signalProcessing as SP
from plot_variance import plot_variance
from process_sim_output import process_sim_output
from balreal import balreal
from welch_psd import welch_psd
# Generate test signal
def test_signal(SCALE, MAXAMP, FREQ, Rng,  OFFSET, t):
    """
    Generate a test signal (carrier)

    Arguments
        SCALE - percentage of maximum amplitude
        MAXAMP - maximum amplitude
        FREQ - signal frequency in hertz
        OFFSET - signal offset
        t - time vector
    
    Returns
        x - sinusoidal test signal
    """
    return (SCALE/100)*MAXAMP*np.cos(2*np.pi*FREQ*t) + OFFSET  +Rng/2 
    # return (SCALE/100)*MAXAMP*np.cos(2*np.pi*FREQ*t) + OFFSET  

HEADROOM  = 0
# %% Chose how to compute SINAD
class sinad_comp:
    CFIT = 1        # curve fitting
    FFT = 2         # fft based 
SINAD_COMP_SEL = sinad_comp.CFIT

# %% Quantiser configurations 
Qconfig = 2
Nb, Mq, Vmin, Vmax, Rng, Qstep, YQ, Qtype = quantiser_configurations(Qconfig)
# %% Sampling frequency and rate
Fs = 1e6
Ts = 1/Fs

# %% Output low-pass filter cutoff order and cutoff frequency
N_lp = 3
Fc_lp = 1e5 # cutoff frequency
# %% Carrier signal
Xcs_SCALE = 100
Xcs_FREQ = 999

# %% Generate time vector

match 1:
    case 1:  # specify duration as number of samples and find number of periods
        Nts = 1e5  # no. of time samples
        Np = np.ceil(Xcs_FREQ*Ts*Nts).astype(int) # no. of periods for carrier

    case 2:  # specify duration as number of periods of carrier
        Np = 3 # no. of periods for carrier
        
Npt = 1  # no. of carrier periods to use to account for transients
Np = Np + Npt

t_end = Np/Xcs_FREQ  # time vector duration
t = np.arange(0, t_end, Ts)  # time vector

# %% Generate carrier/test signal
SIGNAL_MAXAMP = Rng/2 - Qstep  # make headroom for noise dither (see below)
SIGNAL_OFFSET = -Qstep/2  # try to center given quantiser type
Xcs = test_signal(Xcs_SCALE, SIGNAL_MAXAMP, Xcs_FREQ, Rng,  SIGNAL_OFFSET, t)


# %% Quantiser levels : Ideal and Measured
YQns = YQ     
MLns = get_measured_levels(Qconfig)
QMODEL = 1

# %% Scaled 
# Xcs = Xcs/Qstep
# YQns = YQ/Qstep
# MLns = MLns/Qstep
# # YQns = (2**(Nb-1))*YQns # Ideal quantiser levels
# # MLns = (2**(Nb-1))*MLns  # Measured quantiser levels
# Qstep = 1

fig,ax = plt.subplots()
ax.plot(t, Xcs)
# %% Reconstruction filter parameters 
match 1:
    case 1: # LPF derived from optimal NTF Optimal NTF
        # Optimal NTF
        nsf_num = scipy.io.loadmat('generate_optimal_NSF/NSF_num_100kHz_1MHz_3|1Mueta.mat')
        nsf_den = scipy.io.loadmat('generate_optimal_NSF/NSF_den_100kHz_1MHz_3|1Mueta.mat')
        bn = nsf_num['br']
        an = nsf_den['ar']
        b = an.squeeze()
        a = bn.squeeze()

A1, B1, C1, D1 = signal.tf2ss(b,a) # Transfer function to StateSpace
A, B, C, D = balreal(A1, B1, C1, D1)
e, v = np.linalg.eig(A)

# %% MHOQ Binary formulation
N_PRED = 1
SWITCH_MIN = True

match QMODEL:
    case 1:
        QL = YQns.squeeze()
    case 2:
        QL = MLns.squeeze()
# Storage container for code
C_Store = []

# Loop length
len_MPC = Xcs.size - N_PRED

# len_MPC = 1
L = 1e6
# State dimension
x_dim =  int(A.shape[0]) 

# Initial state
init_state = np.zeros(x_dim).reshape(-1,1)

u_kminus1 = np.zeros((2**Nb, N_PRED))

# u_kminus1 = np.round(Xcs[:N_PRED]).reshape(-1,1)
# MPC loop
for j in tqdm.tqdm(range(len_MPC)):

    m = gp.Model("MPC- INL")
    u = m.addMVar((2**Nb, N_PRED), vtype=GRB.BINARY, name= "u") # control variable
    x = m.addMVar((x_dim*(N_PRED+1),1), vtype= GRB.CONTINUOUS, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "x")  # State varible 
    Ns = m.addMVar((1,1), lb = 0,  name = "ns")
    w = m.addMVar((2**Nb,1), vtype= GRB.BINARY,  name="w")
    v = m.addMVar((2**Nb,1), vtype= GRB.BINARY,  name="v")

    # Add objective function
    Obj = 0

    # Set initial constraint
    m.addConstr(x[0:x_dim,:] == init_state)
    for i in range(N_PRED):
        k = x_dim * i
        st = x[k:k+x_dim]
        bin_con =  QL.reshape(1,-1) @ u[:,i].reshape(-1,1) 
        con = bin_con - Xcs[j+i]

        # Switching set 
        # Objective update
        e_t = C @ x[k:k+x_dim] + D * con
        Obj = Obj + e_t * e_t 

        # Constraints update
        f_value = A @ st + B * con
        st_next = x[k+x_dim:k+2*x_dim]
        m.addConstr(st_next == f_value)

        # Binary varialble constraint
        consi = gp.quicksum(u[:,i]) 
        m.addConstr(consi == 1)
        # m.update
        
        y_t = bin_con
        y_t_minus1 = QL.reshape(1,-1) @ u_kminus1

        m.addConstr((y_t - y_t_minus1)/Ts <= L)
        m.addConstr( ( y_t_minus1 - y_t)/Ts <= L)

    # if SWITCH_MIN: 
        Sw2 =  u[:,i].reshape(-1,1) - u_kminus1[:,i].reshape(-1,1)
        Ns = 0
        for i in range(0, Sw2.size):
            Ns = Ns + Sw2[i,0]* Sw2[i,0]
            # Ns = Ns + gp.abs_(Sw2[i,0])
        Obj = Obj + 0*Ns




    m.update
    # Set Gurobi objective
    m.setObjective(Obj, GRB.MINIMIZE)

    # 0 - Supress log output, 1- Print log outputs
    m.Params.LogToConsole = 0

    # Gurobi setting for precision  
    m.Params.IntFeasTol = 1e-9
    m.Params.IntegralityFocus = 1

    # Optimization 
    m.optimize()

    # Extract variable values 
    allvars = m.getVars()
    values = m.getAttr("X",allvars)
    values = np.array(values)

    # Variable dimension
    nr, nc = u.shape
    u_val = values[0:nr*nc]
    u_val = np.reshape(u_val, (2**Nb, N_PRED))

    # Extract Code
    C_MPC = []
    for i in range(N_PRED):
        c1 = np.nonzero(u_val[:,i])[0][0]
        c1 = int(c1)
        C_MPC.append(c1)
    C_MPC = np.array(C_MPC)
    C_Store.append(C_MPC[0])

    U_opt = QL[C_MPC[0]] 
    # State prediction 
    con = U_opt - Xcs[j]
    x0_new =  A @ init_state + B * con
    # State update for subsequent prediction horizon 
    init_state = x0_new
    u_kminus1 = u_val
C_MHOQ  = np.array(C_Store).reshape(1,-1)

match QMODEL:
    case 1:
        Xcs_MHOQ = generate_dac_output(C_MHOQ, YQns)
    case 2:
        Xcs_MHOQ = generate_dac_output(C_MHOQ, MLns) 

# %% switch counter
S_counter  = 0
for i in range(1, Xcs_MHOQ.size):
    if Xcs_MHOQ[0,i-1] != Xcs_MHOQ[0,i]:
        S_counter += 1
print("Total number of switches:",S_counter)
# # %%
# sl = C_MHOQ.size
sl = 1000
fig, ax = plt.subplots()
ax.plot(t[0:sl], Xcs.squeeze()[0:sl])
ax.plot(t[0:sl], Xcs_MHOQ.squeeze()[0:sl])
ax.grid()
# %%
tm = t[:Xcs.size - N_PRED]

TRANSOFF = np.floor(Npt*Fs/Xcs_FREQ).astype(int)  # remove transient effects from output
sp = SP(Xcs, b, a, TRANSOFF)
ym = Xcs_MHOQ.squeeze()
ym_avg, ENOB_M = process_sim_output(tm, ym, Fc_lp, Fs, N_lp, TRANSOFF, SINAD_COMP_SEL, True, 'linear')
# %%
