rm -r ../case_template
mkdir ../case_template
cp -t ../case_template/ acldata000 CD* CL* *.sh control.dat FOIL* fort.* fst LES.IN turbine_dy_param.dat *.inp Urefdata000 *.template
cp -r ModelTurbine ../case_template
cp -r python ../case_template
cp -t ../case_template/ *.lay *.mcr *.py readme

cp -t ../case_template/ restart* z_grid.inp grid.h5
cp -t ../case_template/ post*

# actuator surface model
cp -t ../case_template/ acsdata000 admrdata000

cp -t ../case_template/ nacelle0*
