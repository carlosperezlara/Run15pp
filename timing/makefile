all:
	rootcint -f Dict.cpp -c LinkDef2.h
	g++ -o adc adc.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tdc2 tdc2.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tdc3 tdc3.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tof1 tof1.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tof2 tof2.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tof3 tof3.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tof4 tof4.cpp Dict.cpp `root-config --cflags --glibs`
	g++ -o tof5 tof5.cpp Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
