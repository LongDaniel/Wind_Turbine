#SOURCEDIR=/home/shenl/plyu/code/note_LES
SOURCEDIR=/home/plyu/new_home_test/note_LES
#SOURCEDIR=/home/pin/LES_HOS_Turbine/new_source_codes/note_LES
#SOURCEDIR=/u/zengx372/Pin/ShenCode/note_LES
EXECDIR=$(pwd)
cd $SOURCEDIR/src
make -j4 > log.make 2>&1
tail log.make
cp ./fst $EXECDIR 
echo 'Executable file moved to pwd.'
cp ./log.make $EXECDIR/log.make.les
cd $EXECDIR


#SOURCEDIR=/home/shenl/plyu/code/note_HOS
#SOURCEDIR=/home/plyu/ShenCode/note_HOS
#SOURCEDIR=/home/pin/LES_HOS_Turbine/new_source_codes/note_HOS
#SOURCEDIR=/u/zengx372/Pin/ShenCode/note_HOS
#EXECDIR=$(pwd)
#cd $SOURCEDIR/bin
#make -j4 > log.make 2>&1
#tail log.make
#cp ./main_hos $EXECDIR 
#echo 'Executable file moved to pwd.'
#cp ./log.make $EXECDIR/log.make.hos
#cd $EXECDIR
