rm -f code_ifem code_ifem_dbg code_ifem_asan res_form_fin_*.dat essai.dat exact2.dat interpol2.dat fort.* nodes10.dat xi.dat tetra.dat
gfortran -O2 -o code_ifem *.f && ./code_ifem
