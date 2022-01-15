echo "Step 1: generate fort.11.template, control.dat, Turbine.inp (python)"
python setup_case.py

echo "Step 2: generate fort.11.first, fort.11.second (shell sed)"
cp fort.11.template fort.11.first
sed 's/0\ #\ ISTART/1\ #\ ISTART/g' fort.11.template > fort.11.second

echo "Step 3: please mannually modify core number in run_first.sh and run_second.sh"

echo "Step 4: please mannually run copy.sh"
