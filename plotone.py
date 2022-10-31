#!/usr/bin/env python
'''
usage:
  ./plotone.py [options] <input> <output>
  ./plotone.py [options] (--show|-S) <input> 

options:
  --help -h           Help.
  --verbose -v        Verbose.
  --field=Q -Q        Field to plot. [default: Ey]
  --lims=L -L L       Cbar limits. [default: (-1.0,1.0)]
  --linthresh=L       Linear threshold. [default: 1e-3]
  --slice=S -t S      Set time slice. [default: 0]
  --equal -E          Make plot equal.
  --show -S           Show, don't plot.
'''
from docopt import docopt;
import numpy as np;
opts = docopt(__doc__,help=True);
s = int(opts['--slice']);
with np.load(opts['<input>']) as f:
    d = f[opts['--field']][s,:,:]
    xs ,ys  = f['xs'], f['ys'];
    xhs,yhs = f['xhs'],f['yhs'];
dims = dict(
    Ex = (xhs, ys),
    Ey = (xs , yhs),
    Ez = (xs , ys),
    Bx = (xs , yhs),
    By = (xhs, ys),
    Bz = (xhs, yhs),)[opts['--field']];
import matplotlib.pyplot as plt;
from lspplot.pc import pc;
from pys import parse_ftuple;

lims = parse_ftuple(opts['--lims'],length=2);
linth= float(opts['--linthresh']);
pc(d.T,p=dims,cmap='RdYlBu',lims=lims,log=True,linthresh=linth);
plt.grid();
if opts['--equal']:
    plt.axis('equal');
plt.savefig(opts['<output>'],dpi=200);
