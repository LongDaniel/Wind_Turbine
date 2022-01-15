#module load tecplot
mkdir -p export
tec360 -b -p hosles04_1.mcr > log.tec 2>&1 &

#for i in $(seq -f "%010g" 15000 100 20000)
#do
#  echo ''"${i}"''
#  sed 's/PLYUTS/'"${i}"'/g' q.mcr.template > q.mcr
#  tec360 -b -mesa q.mcr
#done
