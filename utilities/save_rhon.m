function save_rhon(dir_str,rn)
%this function saves rho_n called from spdm block of driver_spdm.m to path
%given by dir_str
  rho_n = rn;
save(dir_str,'rho_n');
end


