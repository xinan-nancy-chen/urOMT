function rn = load_rhon(dir_str)
%this function loads rho_n called from spdm block of driver_spdm.m from path
%given by dir_str
rn = load(dir_str);
rn = rn.rho_n;
