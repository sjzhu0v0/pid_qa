DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow/include
FLAGS_INCLUDE=-I$(DIR_BASE)/include -I$(DIR_BASE)/macro
FLAGS_ROOT=$(shell root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit

all: \
	macro/SkimTreeReading.exe

macro/SkimTreeReading.exe: macro/SkimTreeReading.C
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)