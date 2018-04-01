all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	#g++ -o pizeropizero pizeropizero.cpp Dict.cpp `root-config --cflags --glibs`
	#g++ -o pizeromass pizeromass.cpp Dict.cpp `root-config --cflags --glibs`
	#g++ -o tracks tracks.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o ep eventplane.cpp Dict.cpp `root-config --cflags --glibs`

	rm Dict.*
