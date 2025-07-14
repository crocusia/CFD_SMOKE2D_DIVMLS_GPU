#include "DIVMLS.h"
DIVMLS::DIVMLS() {}
DIVMLS::~DIVMLS() {}

void DIVMLS::initMLS() {
	_scale = make_double3(1.0 / double(RES_SMOKE + 2), 1.0 / double(RES_SMOKE + 2), 0.0);
	_velocityH = new double3[_size];
	_indexH = new int3[_size];

	checkCudaErrors(cudaMalloc((void**)&_velocityD, _size * sizeof(double3)));
	checkCudaErrors(cudaMalloc((void**)&_indexD, _size * sizeof(int3)));
	checkCudaErrors(cudaMalloc((void**)&_posD, _size * sizeof(double3)));
	checkCudaErrors(cudaMalloc((void**)&_btIndexD, _size * sizeof(int3)));
	checkCudaErrors(cudaMalloc((void**)&_btPosD, _size * sizeof(double3)));
	checkCudaErrors(cudaMalloc((void**)&_resultD, _size * sizeof(double3)));
	checkCudaErrors(cudaMalloc((void**)&_nodes, _size * sizeof(NODE)));
	checkCudaErrors(cudaMalloc((void**)&_svd, _size * sizeof(SVD)));
}

void DIVMLS::zeroMLS() {
	checkCudaErrors(cudaMemset(_indexD, 0, _size * sizeof(int3)));
	checkCudaErrors(cudaMemset(_posD, 0, _size * sizeof(double3)));
	checkCudaErrors(cudaMemset(_btIndexD, 0, _size * sizeof(int3)));
	checkCudaErrors(cudaMemset(_btPosD, 0, _size * sizeof(double3)));
	checkCudaErrors(cudaMemset(_resultD, 0, _size * sizeof(double3)));
	checkCudaErrors(cudaMemset(_nodes, 0, _size * sizeof(NODE)));
	checkCudaErrors(cudaMemset(_svd, 0, _size * sizeof(SVD)));
}

void DIVMLS::freeMLS() {
	delete[]	_velocityH;
	delete[]	_indexH;
	checkCudaErrors(cudaFree(_velocityD));
	checkCudaErrors(cudaFree(_indexD));
	checkCudaErrors(cudaFree(_posD));
	checkCudaErrors(cudaFree(_btIndexD));
	checkCudaErrors(cudaFree(_btPosD));
	checkCudaErrors(cudaFree(_resultD));
	checkCudaErrors(cudaFree(_nodes));
}

void DIVMLS::copyHtoD(double* ux, double* uy) {
	//dataSet
	for (int i = 0; i < RES_SMOKE + 2; i++) {
		for (int j = 0; j < RES_SMOKE + 2; j++) {
			int idx = i + (RES_SMOKE + 2) * j;
			double3 velocity = make_double3(ux[idx], uy[idx], 0.0);
			int3 pos = make_int3(i, j, 0);
			_velocityH[idx] = velocity;
			_indexH[idx] = pos;
		}
	}
	//copyData
	checkCudaErrors(cudaMemcpy(_velocityD, _velocityH, _size * sizeof(double3), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(_indexD, _indexH, _size * sizeof(int3), cudaMemcpyHostToDevice));
}

void DIVMLS::copyDtoH(double* dx, double* dy) {
	checkCudaErrors(cudaMemcpy(_velocityH, _resultD, _size * sizeof(double3), cudaMemcpyDeviceToHost));
	for (int i = 0; i < _size; i++) {
		dx[i] = _velocityH[i].x;
		dy[i] = _velocityH[i].y;
		//printf("%f, %f\n", dx[i], dy[i]);
	}
}
