all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	g++ -o Run_BBC_EPC BBC_EPC.cpp AT_BBC_EPC.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	#g++ -o Run_PIDFlow PIDFlow.cpp AT_PIDFlow.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
