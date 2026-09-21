CXX = g++
CXXFLAGS = -std=c++20 -Wall -O2
INCLUDES = -I/ucrt64/include
LIBS = -L/ucrt64/lib -lraylib -lopengl32 -lgdi32 -lwinmm

SRC := $(shell find . -name "*.cpp")
TARGET = AtomicSim
LOGTARGET := $(wildcard atoms-*.csv)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)
	rm -f $(LOGTARGET)
