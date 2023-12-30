# coding:utf-8
#!/usr/bin/env python3

# Pressure coefficient

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/04

import numpy as np

tiny_value = 1.e-50


def cp_isentropic(gamma, mach):
 
  cp = 2.0/(gamma*mach**2)*( (1.0+(gamma-1)/2.0*mach**2)**(gamma/(gamma-1.0)) - 1.0 )

  return cp


def cp_rayleigh(gamma, mach):
  
  var1 = ( (gamma+1.0)/2.0*mach**2 )**( gamma/(gamma-1.0) ) 
  var2 = ( (gamma+1.0)/(2.0*gamma*mach**2 - gamma+1.0 ) )**( 1.0/(gamma-1.0) )
  cp = 2.0/(gamma*mach**2)*( var1*var2 - 1.0 )

  return cp


def main():

  gamma = 1.66

  # isentropic
  mach_isent = np.linspace(0.0, 5, 100)
  for i in range(0,len(mach_isent)):
    if mach_isent[i] <= 0.0:
      mach_isent[i] = mach_isent[i] + tiny_value

  cp_isent = cp_isentropic(gamma, mach_isent)

 # Output
  filename_tmp = 'tecplot_cp_isent.dat'
  header = 'Variables = mach, cp \n zone t=Isentropic i= '+str(len(mach_isent))+' f=point'
  output_data = np.c_[ mach_isent, 
                       cp_isent
                      ]
  delimiter = '\t'
  comments = ''
  np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )


  # Shock wave
  mach_shock = np.linspace(0.5, 10, 200)
  cp_shock = cp_rayleigh(gamma, mach_shock)

 # Output
  filename_tmp = 'tecplot_cp_shock.dat'
  header = 'Variables = mach, cp \n zone t=Shockwave i= '+str(len(mach_shock))+' f=point'
  output_data = np.c_[ mach_shock, 
                       cp_shock
                      ]
  delimiter = '\t'
  comments = ''
  np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()