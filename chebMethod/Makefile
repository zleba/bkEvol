#ROOT libraries
rootCFLAGS := $(shell root-config --cflags)
rootLibs := $(shell root-config --libs)


OBJ=obj/cheb.o obj/integration.o obj/main.o obj/nodes.o obj/Spline.o obj/alphaSpline.o
cheb: ${OBJ}
	g++ -g -O2 $^ -o $@  $(rootLibs)

testmain: obj/testmain.o obj/tests.o obj/cheb.o obj/nodes.o  obj/Spline.o obj/alphaSpline.o
	g++ -g  -O2 $^ -o $@

obj/%.o: src/%.cpp
	g++ -g -O2 -Iinc  $(rootCFLAGS)   -c $^ -o $@
