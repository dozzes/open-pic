CC=g++

CFLAGS=-std=c++11 -c -Wall -ffast-math -O3 -fopenmp

LDFLAGS=-fopenmp -llua -lluabind -ldl

SOURCES=main.cpp \
	bind_to_lua.cpp \
	check_particle.cpp \
	gather_scatter.cpp \
	grid_filters.cpp \
	particle_groups.cpp \
	save_particles.cpp \
	use_lua.cpp \
	call_lua_function.cpp \
	config.cpp \
	grid.cpp \
	io_utilities.cpp \
	marker_particles.cpp \
	particles.cpp \
	simulate.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=opic

all: $(SOURCES) $(EXECUTABLE)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
