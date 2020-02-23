foreach ii { 50.0 48.5 47.0 45.5 44.0 42.5 41.0 39.5 } { 
    mol new $ii/ionized.parm7
    mol addfile $ii/ionized.rst7
    mol addfile $ii/06_Prod.rst type netcdf
}
