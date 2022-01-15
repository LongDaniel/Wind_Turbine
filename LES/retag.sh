rm tags
rm TAGS
#ctags --fortran-kinds=+i * 
cd 2decomp
rm tags
rm TAGS
#ctags --fortran-kinds=+i *.*90
cd ..
ctags -R --fortran-kinds=+i --exclude=backup
