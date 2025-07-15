#pragma once
#include "DIVMLSKernel.h"
#include "DIVMLS.h"

class DIVMLSKernel {
public:
	DIVMLSKernel();
	~DIVMLSKernel();
public:
	void simulation(DIVMLS divmls);
};