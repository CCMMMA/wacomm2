NETCDF_LIBS = -lnetcdf_c++4 -lnetcdf
WACOMM = main.cpp
CXX = g++
CXXFLAGS = -Wall -ggdb
CPPFLAGS = -I$(BOOST_ROOT)/include -I$(NETCDF)/include -I$(NETCDF_CXX4)/include
LDFLAGS = -L$(BOOST_ROOT)/lib -L$(NETCDF)/lib -L$(NETCDF_CXX4)/lib

wacomm.exe: $(WACOMM) jsoncpp.cpp source.cpp particle.cpp
	        $(CXX) $(CPPFLAGS) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ $(NETCDF_LIBS)
