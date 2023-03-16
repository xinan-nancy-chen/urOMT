function un = load_un(dir_str)
%this function loads u_n called from spdm block of driver_spdm.m from path
%given by dir_str
un = load(dir_str);
un = un.u;
