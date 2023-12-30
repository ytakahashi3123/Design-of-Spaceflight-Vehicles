# coding:utf-8
#!/usr/bin/env python3

# Oblique shock wave: relations between shock angle and deflection angle

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/06

import numpy as np
import os as os
import shutil as shutil


def make_directory(dir_path):

  if not os.path.exists(dir_path):
    os.mkdir(dir_path)
  
  return

def check_file_exist(dir_path):

  if os.path.exists(dir_path):
    flag_file_exist = True
  else: 
    flag_file_exist = False

  return flag_file_exist
  

def cal_theta(gamma, mach, beta):
 
  tan_theta = 2.0/np.tan(beta)*(mach**2 * np.sin(beta)**2 - 1.0) / (mach**2 * (gamma+ np.cos(2.0*beta))+2)
  theta = np.arctan(tan_theta)

  return theta


def cal_relation_mach(gamma, mach, beta, theta, mach2):
 
  mach2_tmp = ( 1.0+((gamma-1.0)/2.0) * mach**2 * np.sin(beta)**2 ) / (gamma*mach**2 * np.sin(beta)**2 - ((gamma-1.0)/2.0) )
  mach2_square = mach2_tmp/(np.sin(theta-beta)**2)

  for i in range(0, len(beta)):
    if mach2_square[i] <= 0.0 :
      mach2[i] = 0.0
    else :
      mach2[i] = np.sqrt(mach2_square[i])

  return mach2


def main():

  gamma = 1.40
  mach1 = [1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0, 1000.0]

  beta_degree = np.linspace(1.0, 90.0, 100)
  beta = beta_degree*np.pi/180.0

  for n in range(0, len(mach1)):
    mach = mach1[n]

    theta = cal_theta(gamma, mach, beta)
    theta_deree = theta*180.0/np.pi

    array_tmp = len(beta_degree)
    mach2 = np.zeros(array_tmp).reshape(array_tmp)
    mach2 = cal_relation_mach(gamma, mach, beta, theta, mach2)

  # Output
    dirname = './mach'+str(mach)
    if not check_file_exist(dirname) :
      make_directory(dirname)

    casename='Oblique_mach'+str(mach)

    filename_tmp = dirname+'/'+'tecplot_deflection.dat'
    header = 'Variables = beta, theta, mach2  \n zone t='+casename+' i= '+str(len(beta_degree))+' f=point'
    output_data = np.c_[ beta_degree,
                        theta_deree,
                        mach2
                        ]
    delimiter = '\t'
    comments = ''
    np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()