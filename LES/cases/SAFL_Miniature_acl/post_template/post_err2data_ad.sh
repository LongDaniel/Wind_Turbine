echo "Extract C_Thrust from 'err1' to 'log.cthrust'"
cat err1 | grep 'ud_instant' > log.ud
echo "Plot"
python post_coeff_single_ad.py 
