TARGET = simulation
CC = g++
CFLAGS = -Wall -Wextra -std=c++11 -fopenmp
SRCS = main.cpp mnt.cpp fonctions.cpp
OBJS = $(SRCS:.cpp=.o)
LIBS = -fopenmp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

