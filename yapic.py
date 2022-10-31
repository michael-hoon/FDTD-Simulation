#!/usr/bin/env python
'''
Usage:
    ./yapic.py [options] [<output>]
'''
from docopt import docopt;
import numpy as np;
from time import sleep;

opts = docopt(__doc__,help=True);
outname = opts['<output>']
if not outname:
    outname='dat.npz';

#
#physical constants
#
ut = 1e9;
c = 2.99792458e10/ut;
mu0 = 4*np.pi*1e-7;
e0  = 1/(mu0*2.998e8**2)

dim  = 2

#
# problem setup
#
lm      =  0.8e-4; #wavelength

start_t =    0.0;
#end_t   =   80.0e-6;
end_t   =   40.0e-6;
dt      =   20.0e-9

#laser amplitude
E0      =  1.0
#"width" of laser in space
width   = 20.0*lm;
t0      =  0.0; #not used
phi     =  0.0; #not used

#dimensions in cm...yes, cm
xlim = [ -25.0e-4,  25.0e-4,
         -25.0e-4,  25.0e-4];
#cell subdivision in space
xspacing = [500,500];
dxs = np.array([
    (xlim[2*i+1] - xlim[2*i])/xspacing[i] for i in range(dim)
]);
dus = [ c*dt/dx for dx in dxs ]
dV = dxs[0]*dxs[1];

posx = np.mgrid[
    xlim[0]:xlim[1]:(xspacing[0]+1)*1j,
    xlim[2]:xlim[3]:(xspacing[1]+1)*1j];
posxh=np.array([
    ix + idx/2.0 for ix,idx in zip(posx,dxs) ]);

#shapes
sh = posx.shape[1:];

E=dict();
E['x'] = np.zeros(sh);
E['y'] = np.zeros(sh);
E['z'] = np.zeros(sh);
B=dict();
B['x'] = np.zeros(sh);
B['y'] = np.zeros(sh);
B['z'] = np.zeros(sh);

def destr(d, *l):
    '''destructure a dict'''
    if type(l[0]) is not str:
        l=l[0];
    return [d[i] for i in l];


def xarange(st,en,dx):
    N = int(np.round((en-st)/dx + 1));
    for i in range(N):
        yield dx*i + st;

dF = { dim:np.zeros(xspacing) for dim in 'xy' };

def sminmax(F):
    return f"{np.min(F):.2e}<->{np.max(F):.2e}";

def mkgauss_spatial2D(pos,lm=lm,w0=1.86e-4,phi=0.0,fp=[0,0,0],dir=[1,0,0]):
    x,y = pos;
    x-=fp[0];
    y-=fp[1];
    #z-=fp[2];
    xf = x*dir[0] + y*dir[1];# + z*dir[2];
    xr = np.pi*w0**2/lm;
    wz=np.sqrt(1+(xf/xr)**2)*w0;
    invR = xf/(xf**2+xr**2);
    dsq = x**2 + y**2; # + z**2;
    rs = dsq - xf**2;
    k = 2*np.pi/lm;
    out=dict(
        sp= np.sqrt(w0/wz)*np.exp(-rs/(wz*wz)),
        ph= k*xf + k*invR*rs/2.0 - np.arctan(xf/xr) + phi*np.pi,
        ts= xf/c - invR*rs/2.0/c,
        w = k*c,);
    return out;

def mkgauss_col2D(pos,lm=lm,w0=1.86e-4,fp=[0,0,0],dir=[1,0,0]):
    x,y = pos;
    x-=fp[0];
    y-=fp[1];
    #z-=fp[2];
    xf = x*dir[0] + y*dir[1];# + z*dir[2];
    xr = np.pi*w0**2/lm;
    wz=np.sqrt(1+(xf/xr)**2)*w0;
    invR = xf/(xf**2+xr**2);
    dsq = x**2 + y**2; # + z**2;
    rs = dsq - xf**2;
    k = 2*np.pi/lm;
    out=dict(
        sp= np.exp(-rs/(wz*wz)),
        ph= 0.0*xf, # spatial phase
        ts= 0.0*xf, # temporal shift
        w = k*c,);
    return out;

def mksincol2d(pos,lm=lm,width=width,dir=[1,0,0]):
    x,y = pos;
    xf = x*dir[0] + y*dir[1];
    dsq= x**2 + y**2;
    rs = dsq - xf**2;
    r  = np.sqrt(rs);
    k  = 2*np.pi/lm;
    space_phase = xf*k
    space_env   = np.cos(r*np.pi/width); # spacial envelope
    space_env[r > width/2.0] = 0.0;
    return dict(sph = space_phase, env=space_env, om = k*c);

def sincol2d(t,width,dat):
    sph,env,om  = destr(dat,'sph','env','om');
    tenv = np.sin(t*c*np.pi/width);
    if c*t > width: tenv = 0.0;
    return np.cos(sph-om*t)*env*tenv;

def gauss(spatial,f,t):
    sp,ph,ts,w = destr(spatial,'sp','ph','ts','w');
    return sp*np.cos(ph - w*t)*f(t+ts);

# the basic trend is for Ei, the j and k are x, the i is xh.
# indexes of xmin on E are E[*][0,:,:]
# just do E['y'], lin pol

Ii = 0;
xbnd = posx[1][:,Ii],posx[0][:,Ii];
sindat = mksincol2d(xbnd, lm=lm, width=width,dir=[1,0,0])
print("space allocated");
import matplotlib.pyplot as plt;
from lspplot.pc import pc;

def mksin(E0,t0,T):
    def _f(t):
        t+=t0;
        out = np.zeros_like(t);
        cut = np.logical_or(t < 0, t > T)
        out[cut] = 0;
        out[~cut]= np.sin(t/T*np.pi)[~cut];
        return out*E0;
    return _f;
def mkTgauss(E0,tau,t0,tcut):
    def _f(t):
        v   = np.exp(-((t-t0)/tau)**2)*E0;
        cut = np.abs(t-t0)>tcut;
        v[ cut] = 0;
        return v;
    return _f;

#pulse = mkTgauss(E0,tau,t0,tcut);

def abc(F,nF,axis,i0,i1,dx):
    I0 = [slice(None),slice(None)];
    I0[axis] = i0;
    I1 = [slice(None),slice(None)];
    I1[axis] = i1;
    for k in 'xy':
        F[k][I0] = F[k][I1] + (c*dt - dx)/(c*dt + dx)*(nF[k][I1] - nF[k][I0]);
    return;

doplot = True;
outfmt = '{:04}EM';
acc = [];
def output(t):
    print(f"saving time = {t}");
    cur = dict(t=t);
    for dim in 'xyz':
        cur[f'E{dim}'] = np.copy(E[dim]);
        cur[f'B{dim}'] = np.copy(B[dim]);
    acc.append(cur);

def plot():
    plt.clf();
    plt.semilogy(posx[0]*1e4,E['y']);
    #plt.show();        
outputs = np.arange(000,1001);
print("starting simulation");
for i,t in enumerate(xarange(start_t,end_t,dt)):
    print(f"...{i}");
    #maxwell's correction to ampere
    E['z'][ 1:-1, 1:-1] += (                                        \
        +(B['y'][ 1:-1, 1:-1] - B['y'][ 1:-1,  :-2])*dus[1]         \
        -(B['x'][ 1:-1, 1:-1] - B['x'][  :-2, 1:-1])*dus[0]);
    E['y'][  :  , 1:-1] += -((B['z'][  :  , 1:-1]-B['z'][  :  ,  :-2])*dus[1])
    E['x'][ 1:-1,  :  ] +=  ((B['z'][ 1:-1,  :  ]-B['z'][  :-2,  :  ])*dus[0])

    #E boundary condition
    current_edge = sincol2d(t,width,sindat);
    print(f"...max edge -> {current_edge.max()}");
    E['y'][ :, 0] = current_edge;
    E['y'][ :,-1] = 0;
    
    E['z'][ 0, :] = 0;
    E['z'][-1, :] = 0;


    #maxwell, faraday
    B['z'][  :-1,  :-1] -= (
        +(E['y'][  :-1, 1:  ] - E['y'][  :-1,  :-1])*dus[1]
        -(E['x'][ 1:  ,  :-1] - E['x'][  :-1,  :-1])*dus[0]);
    B['y'][  :  ,  :-1] -= ((E['z'][  :  , 1:  ]-E['z'][  :  ,  :-1])*dus[0]);
    B['x'][  :-1,  :  ] -= ((E['z'][ 1:  ,  :  ]-E['z'][  :-1,  :  ])*dus[1]);
    #implicitly left to zero

    if i in outputs:
        output(t);
        if doplot: plot();
pass
plt.show()
# print(f"outputting to {outname}");
# out = dict();
# fst = acc[0];
# for dim in 'xyz':
#    out[f'E{dim}'] = np.array([d[f'E{dim}'] for d in acc]);
#    out[f'B{dim}'] = np.array([d[f'B{dim}'] for d in acc]);
# out['t'] = np.array([d['t'] for d in acc]);
# out['ys'],out['xs']  =  posx;

# out['yhs'],out['xhs'] = posxh;
# out['lm']  = lm;
# np.savez(outname,**out);
