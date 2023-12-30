# coding:utf-8
#!/usr/bin/env python3

# Ratio of total pressure across shock wave

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/04

import numpy as np

tiny_value = 1.e-50


def cal_entropy(gamma, ratio_press):
 
  gamma_ratio_tmp = (gamma+1.0)/(gamma-1.0)

  var1 = (ratio_press)**(1.0/(gamma-1.0))
  var2 = ((ratio_press+gamma_ratio_tmp)/(gamma_ratio_tmp*ratio_press+1.0))**(gamma/(gamma-1.0))
  entropy = np.log(var1*var2)

  return entropy


def main():

  gamma = 1.40
  ratio_press = np.linspace(0.20, 10.0, 50)
  entropy = cal_entropy(gamma, ratio_press)

 # Output
  filename_tmp = 'tecplot_entropy.dat'
  header = 'Variables = p2_by_p1, entropy_change \n zone t=flow i= '+str(len(ratio_press))+' f=point'
  output_data = np.c_[ ratio_press, 
                       entropy
                      ]
  delimiter = '\t'
  comments = ''
  np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()