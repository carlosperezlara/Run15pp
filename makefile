all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	g++ -o Run_PiZero PiZero.cpp AT_PiZero.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
