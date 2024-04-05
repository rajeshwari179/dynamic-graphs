kcore: main.cpp
	g++ -std=c++0x main.cpp graph_methods.h core_maintenance.h -o kcores
clean :
	rm kcores