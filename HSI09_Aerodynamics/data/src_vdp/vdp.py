# coding:utf-8
#!/usr/bin/env python3

# van der Pol oscillation

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/15

import numpy as np
#import matplotlib.pyplot as plt

"""
強制振動
Forced vibration
"""

rad2deg = 180.0/np.pi
deg2rad = np.pi/180.0


# Functions
def f(mu, delta, a_dot, x, v, t):
#    inertia   =  6.07e-4
#    amplitude = -0.5
#    k=amplitude*np.sin(x)/inertia
    return mu*(1-(x/delta)**2)*v - a_dot*x


if __name__ == '__main__':

  #設定ファイルの読み込み

  # Damping coef
  mu_input    = 5.0
  # Angle switching devergene/attenuation
  delta_input = 9.0
  # Cm_st
  a_dot_input = 500.0

  # Start time, s
  t_min       = 0.0
  # End time, s
  t_max       = 10.0
  # Maximum iteration
  iter_max    = 100000

  # initial
  x0          = 1.e-4
  v0          = 1.e-2

# Set parameters
  mu    = mu_input
  delta = (delta_input)*deg2rad
  a_dot = a_dot_input
  del_t = (t_max-t_min)/float(iter_max)


# Start solving van der Pol equation by 4 stage Runge-Kutta method
  x, v = x0, v0

  tpoints = np.arange(t_min, t_max, del_t)
  xpoints = []
  vpoints = []
  mpoints = []

  for t in tpoints:
    xpoints.append(x)
    vpoints.append(v)

    k1v = f(mu, delta, a_dot, x, v, t)*del_t
    k1x = v*del_t

    k2v = f(mu, delta, a_dot, x+k1x/2, v+k1v/2, t+del_t/2 )*del_t
    k2x = (v+k1v/2)*del_t 

    k3v = f(mu, delta, a_dot, x+k2x/2, v+k2v/2, t+del_t/2 )*del_t
    k3x = (v+k2v/2)*del_t 

    k4v = f(mu, delta, a_dot, x+k3x, v+k3v, t+del_t )*del_t
    k4x = (v+k3v)*del_t 

    v += (k1v + 2*k2v + 2*k3v + k4v)/6
    x += (k1x + 2*k2x + 2*k3x + k4x)/6

    m = f(mu, delta, a_dot, x, v, t)
    mpoints.append(m)

  xpoints_deg = [n*rad2deg for n in xpoints]
  vpoints_deg = [n*rad2deg for n in vpoints]

  # Output: Tecplot format
  file_result = 'tecplot_vdp.dat'  
  header  = 'variables=Time[s],AoA[deg],AngularVelocity[deg/s],m'
  cp_data = np.c_[ tpoints,
                   xpoints_deg,
                   vpoints_deg,
                   mpoints
                  ]
  np.savetxt(file_result, cp_data, header=header, delimiter='\t', comments='' )
