#pragma once
#include<cuda_runtime.h>
#include<helper_cuda.h>
#include<memory>

#include"TArray.h"

#define RES_SMOKE	128

struct NODE {
	int		_neighbor[9];
	int		_constrain[4];
	double	_div[9];
	double	_P_Vec[(2 + 1) * 4];				//[(dim + divCoeffCount) * constrainCount]
	double	_B_Vec[(2 + 1) * 4][2 * 6];			//[(dim + divCoeffCount) * constrainCount][dim * coeffCount]
	double	_Bt_Vec[2 * 6][(2 + 1) * 4];					//[dim * coeffCount][(dim + divCoeffCount) * constrainCount]
	double	_W_Vec[(2 + 1) * 4][(2 + 1) * 4];	//[(dim + divCoeffCount) * constrainCount][(dim + divCoeffCount) * constrainCount]
	double	_W_Mat[2 * 6][(2 + 1) * 4];					//[dim * coeffCount][(dim + divCoeffCount) * constrainCount]
	double	_B_Mat[2][2 * 6];						//[dim][dim * coeffCount]
};

struct SVD {
	double _A[12 * 12];
	double _b[12];
	double _V[12 * 12];
	double _S[12];
	double _x[12];
};

class DIVMLS {
public:
	DIVMLS();
	~DIVMLS();

public:
	int _sizeNoEdge = RES_SMOKE * RES_SMOKE;
	int _size = (RES_SMOKE + 2) * (RES_SMOKE + 2);
	int _dim = 2;
	int _coeffCount = 6;
	int _divCoeffCount = 1;
	int _constrainCount = 4;
	double _dt;
	double3 _scale;

public:
	//Host Data
	double3*				_velocityH;
	int3*					_indexH;
	//Device Data
	NODE*					_nodes;
	SVD*					_svd;
	double3*				_velocityD;
	int3*					_indexD;
	double3*				_posD;
	int3*					_btIndexD;
	double3*				_btPosD;
	double3*				_resultD;

public:
	void initMLS();
	void freeMLS();
	void zeroMLS();
	void copyHtoD(double* ux, double* uy);
	void copyDtoH(double* dx, double* dy);
};