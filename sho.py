#!/usr/bin/env python

#just libraries
import numpy as np;
import matplotlib.pyplot as plt;

#physical quantities
T     = 1.0       #oscillation period
omega = 2*np.pi/T #angular frequency
dt    = 0.01      #timestep
n     = 100       #number of steps

x0    = 1.0     #initial position
vhalf = 0.0     #initial velocity, start at rest

# storage

xn    = np.zeros(n); #simulated positions
vn    = np.zeros(n); #simulated velocities
                     #remember, these are at half steps
x = x0;
v = vhalf;

omega_sq = omega**2  #not that it matters in 2022, but
                     #save the square to avoid
                     #recalculating it every time.

print(f"simulate simple harmonic oscillator for "
      f"{n} steps with timestep {dt:.2e}");

for i in range(n): #c-like loop
    xn[i] = x; 
    vn[i] = v; #minor nitpick, list with append is
               #"more pythonic" but slower, grumble grumble.    
    x += v*dt;
    v += -omega_sq*x*dt;

print("done. plotting...");
# here we generate our "analytic solution" which is cosine!
# surprise, surprise!
times = np.arange(n)*dt;
half_times = times + dt/2.0;
x_analytic = np.cos(times*omega)*x0 + np.sin(times*omega)*(
    (vhalf/omega+x0*np.sin(omega*dt/2))/np.cos(omega*dt/2));
v_analytic =-np.sin(half_times*omega)*omega*x0 + np.cos(omega*half_times)*(
    (vhalf + x0*omega*np.sin(omega*dt/2))/np.cos(omega*dt/2));

#plot them
ax1=plt.subplot(211);
plt.scatter(times, xn, marker='x',
            label='simulated positions')
plt.plot(times, x_analytic, c='orange',
         label='analytic position solution');
plt.legend();
plt.subplot(212,sharex=ax1);
plt.scatter(half_times, vn, marker='+',
            label='simulated velocities')
plt.plot(half_times, v_analytic, c='orange',
         label='analytic velocity solution');
plt.legend();
plt.xlabel('time');
plt.show();
