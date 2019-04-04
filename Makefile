CXX = g++
CXXFLAGS = -g -O3
II = -lnetcdf_c++
NC = ncview

simv2.x: simv2.cpp 
	echo $@
	$(CXX) $(CXXFLAGS) -o simv2.x simv2.cpp $(II) -I.

run:
	./simv2.x < param.txt
	$(NC) output.nc

clean:
	rm -f *.o *.x *.nc
