# Makefile

all: CG_programs yaml links

CG_programs: CG BiCGStab

CG:
	g++ CG.cpp -o CG.exe


BiCGStab:
	g++ BiCGStab.cpp -o BiCGStab.exe

links:
	g++ links.cpp -o links.exe

clean:
	rm *.exe
