fn="admr06"

#cp -t ../admr03_all -r POST_U_* export post_png2mp4.sh post_*.inp
for i in $(seq 5 1 6)
do
  #cp -r ../${fn}_$i/export/ ./
  
  #echo ${fn}_$i
  #echo 1D1_0001
  #cp -r  ../${fn}_$i/POST_U_1D1_0001 ./
  #echo 1D1_0005
  #cp -r ../${fn}_$i/POST_U_1D1_0005 ./
  #echo 1D1_0008
  #cp -r ../${fn}_$i/POST_U_1D1_0008 ./
  #
  #echo 2D1_0001
  #cp -r ../${fn}_$i/POST_U_2D1_0001 ./
  #echo 2D1_0005
  #cp -r ../${fn}_$i/POST_U_2D1_0005 ./
  #echo 2D1_0008
  #cp -r ../${fn}_$i/POST_U_2D1_0008 ./
  #
  #echo 2D2_0004
  #cp -r ../${fn}_$i/POST_U_2D2_0004 ./
  #
  #echo 2D3_0001
  #cp -r ../${fn}_$i/POST_U_2D3_0001 ./

  echo ${fn}_$i
  cd ../${fn}_$i
  cp -t ./ ../${fn}_all/post_coeff_single.py ../${fn}_all/post_mplconfig.py
  cp -t ./ ../${fn}_all/post_err2data.sh
  sh post_err2data.sh
  cp coeff.dat ../${fn}_all/coeff_$i.dat
  cp uref.dat ../${fn}_all/uref_$i.dat
  cd ../${fn}_all
done
