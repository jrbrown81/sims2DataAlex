CC = g++
LD = g++
RM = rm -rf

CXXFLAGS = -O3 -funroll-all-loops -Wno-write-strings
ROOTINCLUDEFLAGS = $(shell root-config --cflags) 

LDFLAGS = -lMinuit -lFoam
ROOTLIBFLAGS = $(shell root-config --glibs)

SOURCES = $(wildcard *.cpp)

EXEC1 = sims2Data_v5

.PHONY: all clean $(EXEC1) 

all: $(EXEC1) 
	@echo done.

$(EXEC1): sims2Data_v5.o
	@echo linking $@...
	@$(CC) -o $@ $^ $(LDFLAGS) $(ROOTLIBFLAGS)


sims2Data_v5.o: *.cpp
	@echo compiling $(notdir $<)...
	@$(CC) -o $@ -c $^ $(CXXFLAGS) $(ROOTINCLUDEFLAGS)



clean:
	@echo cleaning...
	@$(RM) *.o $(EXEC1) 
