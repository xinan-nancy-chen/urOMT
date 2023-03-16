function save_un(dir_str,un)
%this function saves u_n called from spdm block of driver_spdm.m to path
%given by dir_str
  u = un;
save(dir_str,'u');
end
