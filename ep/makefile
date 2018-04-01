all:
	rootcint -f Dict.cpp -c ../qcQ.h ../LinkDef.h
	g++ -o ep eventplane.cpp Dict.cpp `root-config --cflags --glibs`

	rm Dict.*
