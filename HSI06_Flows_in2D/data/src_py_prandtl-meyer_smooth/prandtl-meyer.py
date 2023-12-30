# coding:utf-8
#!/usr/bin/env python3

# Prandtl-Meyer function

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/07

import numpy as np
  

def cal_pm_function(gamma, mach):
 
  gamma_tmp = (gamma+1.0)/(gamma-1.0)
  nu = np.sqrt(gamma_tmp) * np.arctan( np.sqrt(1.0/gamma_tmp*(mach**2-1.0) ) ) - np.arctan( np.sqrt(mach**2-1.0) )

  return nu


def main():

  gamma = 1.40
  mach  = np.linspace(1.0, 50.0, 500)
#  mach  = np.linspace(1.0, 10.0, 901)
  nu = cal_pm_function(gamma, mach)
  nu_degree = nu*180.0/np.pi

#  print('Mach number [-], Prandlt-Meyer function [degree], gamma=1.4')
#  for n in range(0,len(mach)):
#    print( '{:.3f}'.format(mach[n]), ',  ', '{:.6f}'.format(nu_degree[n]) )

  filename_tmp = 'tecplot_pm.dat'
  header = 'Variables = mach, nu  \n zone t=Prandtl-Meyer i= '+str(len(mach))+' f=point'
  output_data = np.c_[ mach, 
                       nu_degree
                      ]
  delimiter = '\t'
  comments = ''
  #np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()