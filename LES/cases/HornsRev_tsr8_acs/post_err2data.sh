cat err1 | grep 'C_Thrust' > log.cthrust
cat err1 | grep 'Uref=' > log.uref
python post_coeff_single.py 
