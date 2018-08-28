
all: sce
	./sce

sce: $(wildcard *.f90)
	./quick

sd_fortran.tar.gz: $(wildcard *.f90) $(wildcard *.py) quick CmakeLists.txt README.md Makefile
	tar -cvzf $@ $^

clean:
	rm sce
	rm -rf bld/
	rm sd_fortran.tar.gz
	*.txt
