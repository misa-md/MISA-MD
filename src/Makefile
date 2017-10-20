CC=		mpicxx
CCFLAGS=	-g -c

LINK=		mpicxx
LINKFLAGS=	-g -o

BUILD_DIR = build/
SRCS = main.cpp atom.cpp config.cpp createatom.cpp crystal_md.cpp domain.cpp domaindecomposition.cpp eam.cpp input.cpp integrator.cpp InterpolationObject.cpp latparticledata.cpp mpi_utils.cpp particledata.cpp  simulation.cpp
OBJS = $(patsubst %.cpp,$(BUILD_DIR)%.o,$(SRCS))

all:Crystal-MD

vpath %.cpp ./
vpath %.o $(BUILD_DIR)

Crystal-MD: $(OBJS)
	$(LINK) $(LINKFLAGS) $@ $^

$(BUILD_DIR)%.o:%.cpp
	$(CC) $(CCFLAGS) -o $@ $<

clean:
	rm $(BUILD_DIR)*.o

