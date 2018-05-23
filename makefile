all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	#g++ -o epcalib epcalib.cpp AT_EventPlaneCalibrator.cpp AT_ReadTree.cpp Analysis.cpp qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	g++ -o pid pid.cpp AT_PIDFlow.cpp AT_ReadTree.cpp Analysis.cpp qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
