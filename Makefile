CXX = g++
CXXFLAGS = -std=c++20 -Wall
INCLUDES = -I/ucrt64/include
LIBS = -L/ucrt64/lib -lraylib -lopengl32 -lgdi32 -lwinmm

SRC := $(shell find . -name "*.cpp")
TARGET = ParticleSim

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)
	python clear-logs.py
