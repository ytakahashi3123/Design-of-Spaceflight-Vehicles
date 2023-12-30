# coding:utf-8
#!/usr/bin/env python3

# Normal shock wave in nozzle flow

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2022/07/14

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
  

def set_nozzle_geom(a, b, x):
  y = a*x + b
  return y


def set_nozzle_area_2d(a, b, x):
  # 2D nozzle
  area = 2.0*(a*x + b)
  return area


def cal_machnumber_arbitrary(gamma, area_ratio, mach_arb):
 
  # Calculate Mach number from specific heat ratio and area ration between two arbitrary sections
  # accroding to nozzle flow relations
 
  # Temporal variables
  gamma_power  = (gamma+1.0)/(2.0*(gamma-1.0))
  gamma_power2 = (-gamma+3.0)/(2.0*(gamma-1.0))
  gamma_p = (gamma+1.0)/2.0
  gamma_m = (gamma-1.0)/2.0

  # Calculation parameters for Newton-Raphson method
  iter_max = 100
  omega    = 1.0
  eps      = 1.e-8

  # Initial mach number
  if mach_arb >= 1.0:
    mach = 1.1
  else:
    mach = 0.1
  
  for n in range(0,iter_max):
    var_tmp = ( 1.0 + gamma_m*mach**2 )/( 1.0 + gamma_m*mach_arb**2 )
    f  =   mach_arb/mach    * var_tmp**gamma_power - area_ratio
    df = - mach_arb/mach**2 * var_tmp**gamma_power + (gamma_p*mach_arb)/(1.0 + gamma_m*mach_arb**2 ) * var_tmp**gamma_power2
    mach_prev    = mach
    mach_update  = mach - f/df
    mach         = mach_prev + omega* (mach_update-mach_prev)
    relres       = abs(mach - mach_prev)/mach
    #print(n,mach,relres)
    if relres <= eps :
      break

  return mach


def isentropic_relation_arbitrary(gamma, mach_arb, mach):
  
  # To calculate pressure, density velocity, and temperature ratios between two arbitrary sections
  # from Mach number based on isentropic relations

  # _1: unknow, _arb: known
  # T1/T_arb, p1/p_arb, rho1/rho_arb

  var_tmp = (1.0+(gamma-1.0)/2.0*mach_arb**2)/(1.0+(gamma-1.0)/2.0*mach**2)
  tempeature_ratio = var_tmp**(1.0)
  pressure_ratio   = var_tmp**(gamma/(gamma-1.0))
  density_ratio    = var_tmp**(1.0/(gamma-1.0))

  return tempeature_ratio, pressure_ratio, density_ratio 


def main():

  # Parameter setting
  gamma = 1.40
  pres_stag = 10.e3 # Stagnation pressure, Pa
  pres_back = 300.0 # Back pressure, Pa

  # make directory
  dirname = 'case_gamma_'+str(gamma)
  if not check_file_exist(dirname) :
    make_directory(dirname)
  casename = 'gamma'+str(gamma)

  # Assume that the nozzle is the conical nozzle (nozzle wall line follows: function y=ax+b)
  # --Parameters: nozzle geometry
  y_nozzleexit = 0.05     # m
  y_throat     = 0.0025   # m
  x_nozzleexit = 0.2      # m

  b_nozzle =  y_throat
  a_nozzle = (y_nozzleexit-y_throat)/x_nozzleexit

  # Set x coordinate
  x_nozzle = np.linspace(0.0, 0.2, 100)

  # Nozzle wall line (2D nozzle)
  y_nozzle = set_nozzle_geom(a_nozzle, b_nozzle, x_nozzle)

  # Nozzle area (2D nozzle)
  area_throat = 2.0*y_throat
  area_nozzle = set_nozzle_area_2d(a_nozzle, b_nozzle, x_nozzle)
  area_ratio  = area_nozzle/area_throat
  area_exit   = area_nozzle[-1]

  print('Area ratio',area_ratio)

  # Pressure at throat
  pres_throat = pres_stag/((gamma+1.0)/2.0)**(gamma/(gamma-1.0))

  # ノズル内部の垂直衝撃波発生位置をx_shockとする。x_shock=x_throatのとき、よどみ圧力/背圧の関係はどうなるか？
  # 流れはスロートでチョークしているものとする。
  x_shock = x_nozzle

  num_array = len(x_nozzle)
  mach_array_tmp = np.zeros(num_array*num_array).reshape(num_array,num_array)
  pres_array_tmp = np.zeros(num_array*num_array).reshape(num_array,num_array)

  for n in range(0,len(x_shock)):

    # Flow properties behind the shock wave
    area_ratio_tmp = area_nozzle[n]/area_throat
    mach_tmp = cal_machnumber_arbitrary(gamma, area_ratio_tmp, 1.0)
    #tempeature_ratio, pressure_ratio, density_ratio = isentropic_relation(gamma, mach_tmp)
    tempeature_ratio, pressure_ratio, density_ratio = isentropic_relation_arbitrary(gamma, 1.0, mach_tmp)
    pres1 = pressure_ratio*pres_throat

    pres2 = pres1* (1.0+2.0*gamma/(gamma+1.0)*(mach_tmp**2-1.0) )
    mach2 = np.sqrt( (1.0+(gamma-1.0)/2.0*mach_tmp**2)/(gamma*mach_tmp**2-(gamma-1.0)/2.0) )

    for m in range(0,len(x_nozzle)):
            
      # Isentropic relation in region betwee throat-->shock wave
      if x_nozzle[m] <= x_shock[n] :
        area_ratio_tmp = area_nozzle[m]/area_throat
        mach_tmp = cal_machnumber_arbitrary(gamma, area_ratio_tmp, 1.0)
        tempeature_ratio, pressure_ratio, density_ratio = isentropic_relation_arbitrary(gamma, 1.0, mach_tmp)
        pres_tmp = pressure_ratio*pres_throat
      elif x_nozzle[m] > x_shock[n] :
        area_ratio_tmp = area_nozzle[m]/area_nozzle[n]
        mach_tmp = cal_machnumber_arbitrary(gamma, area_ratio_tmp, mach2)
        tempeature_ratio, pressure_ratio, density_ratio = isentropic_relation_arbitrary(gamma, mach2, mach_tmp)
        pres_tmp = pressure_ratio*pres2

      mach_array_tmp[n,m] = mach_tmp
      pres_array_tmp[n,m] = pres_tmp

  # Output: Tecplot
  filename_tmp = dirname+'/'+'tecplot_nozzleflow.dat'
  num_data = len(x_nozzle)
  newline_code = '\n'
  blank_code=' '
  file = open(filename_tmp, "w")
  for n in range(0,len(x_shock)):
    casename = 'shock_x_'+str(x_shock[n])
    header   = 'Variables = x[m], mach, pressure[Pa], \n zone t='+str(casename)+' i= '+str(num_data )+' f=point'
    file.write( header + newline_code )
    txt_tmp = ''
    for m in range(0,num_data):
      txt_tmp = txt_tmp + str( x_nozzle[m] ) + blank_code + str( mach_array_tmp[n,m] ) + blank_code + str( pres_array_tmp[n,m] ) + newline_code
    file.write( txt_tmp + newline_code )
  file.close()
  #output_data = np.c_[ x_nozzle,
  #                     mach_array_tmp[80,:],
  #                     pres_array_tmp[80,:]
  #                    ]
  #delimiter = '\t'
  #comments = ''
  #np.savetxt(filename_tmp, output_data, header=header, delimiter=delimiter, comments=comments )

  return


if __name__ == '__main__':

  # main 
  main()

  exit()