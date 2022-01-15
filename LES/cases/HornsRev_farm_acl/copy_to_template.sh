rm -r ../case_template
mkdir ../case_template
cp -t ../case_template/ acldata000 CD00 CL00 *.sh control.dat FOIL00 fort.* fst LES.IN turbine_dy_param.dat Turbine.inp Urefdata000
cp -r ModelTurbine ../case_template
cp -r python ../case_template
cp -t ../case_template/ *.lay *.mcr *.py nacelle0* *.stl *.inp
cp -t ../case_template/ restart* grid.h5 z_grid.inp

# actuator surface model
cp -t ../case_template/ acsdata000
cp -t ../case_template/ main_hos

# remove unnecessary files
rm ../case_template/restart*hos*0
