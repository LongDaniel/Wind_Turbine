module load tecplot
tec360 -b -p inflow05_02.mcr > log.tec 2>&1 &
tec360 -b -p inflow05_03.mcr > log.tec 2>&1 &

#for i in $(seq -f "%010g" 15000 100 20000)
#do
#  echo ''"${i}"''
#  sed 's/PLYUTS/'"${i}"'/g' q.mcr.template > q.mcr
#  tec360 -b -mesa q.mcr
#done
