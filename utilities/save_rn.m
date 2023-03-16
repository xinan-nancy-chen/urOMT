function save_rn(dir_str,rn)
%this function saves r called from spdm block of driver_spdm.m to path
%given by dir_str
  r = rn;
save(dir_str,'r');
end


