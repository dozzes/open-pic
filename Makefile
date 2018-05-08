INC = -I/home/dosipyan/opic/3rdparty/boost_1_52_0 \
	-I/usr/include/lua5.1 \
	-I/home/dosipyan/opic/3rdparty/luabind-0.8.1

LIB = -L/home/dosipyan/opic/3rdparty/luabind-0.8.1/stage -L/usr/lib

CC = g++

CFLAGS = -std=c++11 -c -Wall $(INC) -ffast-math -O3 -fopenmp

LDFLAGS = $(LIB) -fopenmp -llua5.1 -lluabind -ldl

SOURCES = src/main.cpp \
	src/bind_to_lua.cpp \
	src/check_particle.cpp \
	src/gather_scatter.cpp \
	src/grid_filters.cpp \
	src/particle_groups.cpp \
	src/save_particles.cpp \
	src/use_lua.cpp \
	src/call_lua_function.cpp \
	src/config.cpp \
	src/grid.cpp \
	src/save_grid.cpp \
	src/io_utilities.cpp \
	src/marker_particles.cpp \
	src/particles.cpp \
	src/simulate.cpp

OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = opic

all: $(SOURCES) $(EXECUTABLE)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
