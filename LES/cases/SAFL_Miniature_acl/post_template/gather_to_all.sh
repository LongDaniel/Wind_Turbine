fn="val_if1_nac_18"

#cp -t ../admr03_all -r POST_U_* export post_png2mp4.sh post_*.inp
for i in $(seq 10 1 33)
do
#  cp -r ../${fn}_$i/export/ ./
  
  echo ${fn}_$i

#  echo 1D1
#  cp -r ../${fn}_$i/POST_UI_1D1_0001 ./
#  cp -r ../${fn}_$i/POST_UI_1D1_0005 ./
#  cp -r ../${fn}_$i/POST_UI_1D1_0008 ./
#  
#  echo 2D2
#  cp -r ../${fn}_$i/POST_UI_2D2_0002 ./
#    
#  echo 2D3
#  cp -r ../${fn}_$i/POST_UI_2D3_0001 ./
#  cp -r ../${fn}_$i/POST_UI_2D3_0002 ./
#  cp -r ../${fn}_$i/POST_UI_2D3_0003 ./
#  cp -r ../${fn}_$i/POST_UI_2D3_0004 ./


  echo ${fn}_$i
  cd ../${fn}_$i
  cp -t ./ ../${fn}_all/post_coeff_single.py ../${fn}_all/post_mplconfig.py
  cp -t ./ ../${fn}_all/post_err2data.sh
  sh post_err2data.sh
  cp coeff.dat ../${fn}_all/coeff_$i.dat
  cp uref.dat ../${fn}_all/uref_$i.dat
  cd ../${fn}_all

#
#  cd ../${fn}_$i
#  cp Nacelle_*.dat ../${fn}_all/
#  cd ../${fn}_all

  
  mkdir -p shortcuts
  cd shortcuts  
  ln -s ../../${fn}_$i/DAT_* ./
  ln -s ../../${fn}_$i/line_*_nf.dat ./
  ln -s ../../${fn}_$i/restart_hos.* ./
  ln -s ../../${fn}_$i/restart_param_hos.* ./
  ln -s ../../${fn}_$i/surfnac_* ./
  ln -s ../../${fn}_$i/surface_* ./
 
  cd ../../${fn}_$i
  cp -t ../${fn}_all/shortcuts/ fort.11 fort.11.second LES.IN *.inp *.sh *.py
  cp -t ../${fn}_all/shortcuts/ acldata* acsdata* admrdata* CD* CL* control.dat FOIL* nacelle* Urefdata* 
  cp -t ../${fn}_all/shortcuts/ restart_uref.dat 
  cd ../${fn}_all
done
