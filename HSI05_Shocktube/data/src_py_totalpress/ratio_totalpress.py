# coding:utf-8
#!/usr/bin/env python3

# Ratio of total pressure across shock wave

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/04

import numpy as np

tiny_value = 1.e-50


def cal_ratio_totalpres(gamma, mach):
 
  var1 = (  (gamma+1.0)/(2.0*gamma*mach**2 - (gamma-1.0)) )**( 1.0/(gamma-1.0) ) 
  var2 = ( ((gamma+1.0)*mach**2)/((gamma-1.0)*mach**2 + 2.0 ) )**( gamma/(gamma-1.0) )
  ratio_totalpres = var1*var2

  return ratio_totalpres


def main():

  gamma = 1.40

  # isentropic
  mach = np.linspace(1.0, 10, 200)
  for i in range(0,len(mach)):
    if mach[i] <= 0.0:
      mach[i] = mach[i] + tiny_value

  ratio_totalpres = cal_ratio_totalpres(gamma, mach)

 # Output
  filename_tmp = 'tecplot_totalpress.dat'
  header = 'Variables = mach, p2_by_p1 \n zone t=flow i= '+str(len(mach))+' f=point'
  output_data = np.c_[ mach, 
                       ratio_totalpres 
                      ]
  delimiter = '\t'
  comments = ''
  np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()