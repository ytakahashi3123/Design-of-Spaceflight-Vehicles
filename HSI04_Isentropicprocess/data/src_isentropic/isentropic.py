# coding:utf-8
#!/usr/bin/env python3

# Nozzle flow relations

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/10

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
  

def cal_machnumber(gamma, area_ratio):
 
  # Calculate Mach number from specific heat ratio and area ration between throat and nozzle
  # accroding to nozzle flow relations
 
  # Temporal variables
  gamma_power  = (gamma+1.0)/(2.0*(gamma-1.0))
  gamma_power2 = (-gamma+3.0)/(2.0*(gamma-1.0))
  gamma_p = (gamma+1.0)/2.0
  gamma_m = (gamma-1.0)/2.0

  # Calculation parameters for Newton-Raphson method
  iter_max = 100
  omega    = 1.0
  eps = 1.e-8

  # Initial mach number
  mach     = 1.1
  
  for n in range(0,iter_max):
    var_tmp = ( 1.0+gamma_m*mach**2 )/gamma_p
    f  =   1.0/mach    * var_tmp**gamma_power - area_ratio
    df = - 1.0/mach**2 * var_tmp**gamma_power +  var_tmp**gamma_power2
    mach_prev    = mach
    mach_update  = mach - f/df
    mach         = mach_prev + omega* (mach_update-mach_prev)
    relres       = abs(mach - mach_prev)/mach
    print(n,mach,relres)
    if relres <= eps :
      break

  return mach


def isentropic_relation(gamma, mach):
  
  # To calculate pressure, density velocity, and temperature ratios between stagnation and certain point
  # from Mach number based on isentropic relations

  # T/T0, p/p0, rho/rho0

  var_tmp = (1.0+(gamma-1.0)/2.0*mach**2)
  tempeature_ratio = var_tmp**(-1.0)
  pressure_ratio   = var_tmp**(-gamma/(gamma-1.0))
  density_ratio    = var_tmp**(-1.0/(gamma-1.0))

  return tempeature_ratio, pressure_ratio, density_ratio 


def main():

  # Parameter setting
  gamma = 1.2
  mach_number = np.linspace(0.1, 20, 500)

  # --make directory
  dirname = 'case_gamma_'+str(gamma)
  if not check_file_exist(dirname) :
    make_directory(dirname)
  casename = 'gamma'+str(gamma)

  # --Get mach number and isentropic relations
  mach_series = []
  temp_ratio_series = []
  pres_ratio_series = []
  dens_ratio_series = []
  for n in range(0,len(mach_number)):
    mach_tmp = mach_number[n]
    mach_series.append(mach_tmp)

    tempeature_ratio, pressure_ratio, density_ratio = isentropic_relation(gamma, mach_tmp)
    temp_ratio_series.append(tempeature_ratio)
    pres_ratio_series.append(pressure_ratio)
    dens_ratio_series.append(density_ratio)

  mach_series_ndarray = np.array(mach_series)
  temp_ratio_ndarray  = np.array(temp_ratio_series)
  pres_ratio_ndarray  = np.array(pres_ratio_series)
  dens_ratio_ndarray  = np.array(dens_ratio_series)

  # --Output
  filename_tmp = dirname+'/'+'tecplot_nozzleflow.dat'
  header = 'Variables = mach, temperature_ratio(TbyT0),pressure_ratio(pbyp0),density_ratio(rhobyrho0) \n zone t='+str(casename)+' i= '+str(len(mach_series))+' f=point'
  output_data = np.c_[ mach_series_ndarray,
                       temp_ratio_ndarray,
                       pres_ratio_ndarray,
                       dens_ratio_ndarray
                      ]
  delimiter = '\t'
  comments = ''
  np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()