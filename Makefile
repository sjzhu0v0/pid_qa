DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow/
DIR_BASE2= /home/szhu/work/alice/analysis/PairFlow
FLAGS_INCLUDE=-I$(DIR_BASE)/include -I$(DIR_BASE2)/include -Wno-narrowing -O2
FLAGS_ROOT=$(shell root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit

all: \
	macro/SkimTreeReading.exe \
	draw/ComparisonSeparationPower.exe \
	draw/RunDependentProfile.exe \
	process/SeparationPower.exe

macro/SkimTreeReading.exe: macro/SkimTreeReading.C
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)

draw/ComparisonSeparationPower.exe: draw/ComparisonSeparationPower.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)

draw/RunDependentProfile.exe: draw/RunDependentProfile.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT)

process/SeparationPower.exe: process/SeparationPower.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)