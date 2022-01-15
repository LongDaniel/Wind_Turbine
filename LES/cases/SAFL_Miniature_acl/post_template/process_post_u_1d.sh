#!/bin/bash

for i in $(seq 1 1 15)
do
  #cp fort.12.template fort.12.first
  #sed 's/0\ #\ ISTART/1\ #\ ISTART/g' fort.12.template > fort.12.second
  #echo $i
  #echo "$i"
  echo "s/isec\ =\ ISEC/isec\ =\ $i/g"
  index=$(printf "%04d" $i)
  echo $index
  sed "s/isec\ =\ ISEC/isec\ =\ $i/g" post_u_1d_s1.py.template > POST_U_1D1_$index/post_u_1d_s1.py
  cd POST_U_1D1_$index
  python post_u_1d_s1.py > post_u_1d_s1.log 2>&1
  cd ..
done
