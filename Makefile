
OPIC_DIR=/home/dosipyan/build_place/opic/open-pic

#FILES = $(shell ls $(OPIC_DIR)/3rdparty/lua-5.1/etc)
#echo $(FILES)

INC = -I$(OPIC_DIR)/3rdparty/boost_1_52_0 \
	-I$(OPIC_DIR)/3rdparty/lua-5.1/src \
	-I$(OPIC_DIR)/3rdparty/lua-5.1/etc \
	-I$(OPIC_DIR)/3rdparty/luabind-0.8.1/lib/include

LIB = -L$(OPIC_DIR)/3rdparty/luabind-0.8.1/stage -L/usr/lib

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
