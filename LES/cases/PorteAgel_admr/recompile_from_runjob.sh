#SOURCEDIR=/home/shenl/plyu/code/note_LES
SOURCEDIR=/home/plyu/ShenCode/note_LES
EXECDIR=$(pwd)
cd $SOURCEDIR/bin
#make clean
make -j4 > log.make 2>&1
tail log.make
cp ./fst $EXECDIR 
echo 'Executable file moved to pwd.'
cp ./log.make $EXECDIR
cd $EXECDIR
