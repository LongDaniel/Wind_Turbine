echo "Extract C_Thrust from 'err1' to 'log.cthrust'"
cat err1 | grep 'C_Thrust' > log.cthrust
echo "Extract Uref from 'err1' to 'log.uref'"
cat err1 | grep 'Uref=' > log.uref
echo "Plot"
python post_coeff_single.py 
