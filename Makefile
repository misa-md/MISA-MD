CC=		mpicxx
CCFLAGS=	-g -c

LINK=		mpicxx
LINKFLAGS=	-g -o

BUILD_DIR = build/
SRCS = main.cpp simulation.cpp domain.cpp domaindecomposition.cpp atom.cpp input.cpp integrator.cpp eam.cpp InterpolationObject.cpp particledata.cpp latparticledata.cpp createatom.cpp
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

