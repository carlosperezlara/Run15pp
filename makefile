all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	g++ -o Run_PiZeroMass PiZeroMass.cpp AT_PiZeroMass.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
