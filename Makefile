DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow/
FLAGS_INCLUDE=-I$(DIR_BASE)/include -I$(DIR_BASE)/macro
FLAGS_ROOT=$(shell root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit -Wno-narrowing -O2

all: \
	macro/SkimTreeReading.exe \
	draw/ComparisonSeparationPower.exe

macro/SkimTreeReading.exe: macro/SkimTreeReading.C
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)

draw/ComparisonSeparationPower.exe: draw/ComparisonSeparationPower.C
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)