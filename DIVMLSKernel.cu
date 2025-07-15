#include "DIVMLSKernel.h"

DIVMLSKernel::DIVMLSKernel() {}
DIVMLSKernel::~DIVMLSKernel() {}

int calcMin(int n, int blockSize) {
	if (n > blockSize) { return blockSize; }
	return n;
}

int iDivUp(int a, int b)
{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

void computeGridSize(int n, int blockSize, int& numBlocks, int& numThreads)
{
	numThreads = calcMin(blockSize, n);
	numBlocks = iDivUp(n, numThreads);
}

__global__ void backTracingD(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		double h0 = divmls._dt * RES_SMOKE;

		//calculate Original Pos
		double px = divmls._indexD[tid].x * divmls._scale.x;
		double py = divmls._indexD[tid].y * divmls._scale.y;
		divmls._posD[tid] = make_double3(px, py, 0.0);

		//backTrace Index
		double x = divmls._indexD[tid].x - h0 * divmls._velocityD[tid].x;
		double y = divmls._indexD[tid].y - h0 * divmls._velocityD[tid].y;
		if (x < 0.5)				x = 0.5;
		if (x > RES_SMOKE + 0.5)	x = RES_SMOKE + 0.5;
		if (y < 0.5)				y = 0.5;
		if (y > RES_SMOKE + 0.5)	y = RES_SMOKE + 0.5;

		int idx = int(x);
		int idy = int(y);
		divmls._btIndexD[tid] = make_int3(idx, idy, 0);

		//backTrace Pos
		x = x * divmls._scale.x;
		y = y * divmls._scale.y;
		divmls._btPosD[tid] = make_double3(x, y, 0.0);

		//set NeighborNode Index
		int n = 0;
		int id = 0;
		int index = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				id = i + 3 * j;
				index = (idx + i - 1) + (RES_SMOKE + 2) * (idy + j - 1);
				divmls._nodes[tid]._neighbor[id] = index;
			}
		}
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				id = i + 2 * j;
				index = (idx + i) + (RES_SMOKE + 2) * (idy + j);
				divmls._nodes[tid]._constrain[id] = index;
			}
		}
	}
}

__global__ void computeDivD(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < divmls._size) {
		int idx = divmls._btIndexD[tid].x;
		int idy = divmls._btIndexD[tid].y;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				int id = i + 3 * j;
				double v0x = 0.0;
				double v1x = 0.0;
				double v2y = 0.0;
				double v3y = 0.0;
				int index0 = divmls._nodes[tid]._neighbor[(i - 1) + 3 * j];	//좌
				int index1 = divmls._nodes[tid]._neighbor[(i + 1) + 3 * j];	//우
				int index2 = divmls._nodes[tid]._neighbor[i + 3 * (j + 1)];	//상
				int index3 = divmls._nodes[tid]._neighbor[i + 3 * (j - 1)];	//하
				if ((i == 0 && j == 0) || (i == 2 && j == 0) || (i == 0 && j == 2) || (i == 2 && j == 2)) {
					if (i == 0 && j == 0) {
						v0x = 0.0;
						v1x = divmls._velocityD[index1].x;
						v2y = divmls._velocityD[index2].y;
						v3y = 0.0;
					}
					if (i == 2 && j == 0) {
						v0x = divmls._velocityD[index0].x;
						v1x = 0.0;
						v2y = divmls._velocityD[index2].y;
						v3y = 0.0;
					}
					if (i == 0 && j == 2) {
						v0x = 0.0;
						v1x = divmls._velocityD[index1].x;
						v2y = 0.0;
						v3y = divmls._velocityD[index3].y;
					}
					if (i == 2 && j == 2) {
						v0x = divmls._velocityD[index0].x;
						v1x = 0.0;
						v2y = 0.0;
						v3y = divmls._velocityD[index3].y;
					}
				}
				else if (i == 0) {
					v0x = 0.0;
					v1x = divmls._velocityD[index1].x;
					v2y = divmls._velocityD[index2].y;
					v3y = divmls._velocityD[index3].y;
				}
				else if (j == 0) {
					v0x = divmls._velocityD[index0].x;
					v1x = divmls._velocityD[index1].x;
					v2y = divmls._velocityD[index2].y;
					v3y = 0.0;
				}
				else if (i == 2) {
					v0x = divmls._velocityD[index0].x;
					v1x = 0.0;
					v2y = divmls._velocityD[index2].y;
					v3y = divmls._velocityD[index3].y;
				}
				else if (j == 2) {
					v0x = divmls._velocityD[index0].x;
					v1x = divmls._velocityD[index1].x;
					v2y = 0.0;
					v3y = divmls._velocityD[index3].y;
				}
				else {
					v0x = divmls._velocityD[index0].x;
					v1x = divmls._velocityD[index1].x;
					v2y = divmls._velocityD[index2].y;
					v3y = divmls._velocityD[index3].y;
				}
				divmls._nodes[tid]._div[id] = -0.5 * (v1x - v0x + v3y - v2y) / RES_SMOKE;
			}
		}
		divmls._resultD[tid].x = divmls._nodes[tid]._div[4];
	}
}

__global__ void precomputeBVec(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		int _dim = divmls._dim;
		int _divCoeffCount = divmls._divCoeffCount;
		int _constrainCount = divmls._constrainCount;
		int _coeffCount = divmls._coeffCount;
		int index = 0;

		for (int i = 0; i < (_dim + _divCoeffCount) * _constrainCount; i++) {
			int id = i / (_dim + _divCoeffCount);
			int cid = divmls._nodes[tid]._constrain[id];
			
			double cx = divmls._posD[cid].x; 
			double cy = divmls._posD[cid].y;
			double Basis[6] = { 1.0, cx, cy, cx * cx, cy * cy, cx * cy };
			double divBasis[12] = { 0.0,1.0,0.0,2.0 * cx,0.0,cy,	// derivative x
				0.0,0.0,1.0,0.0,2.0 * cy,cx };	// derivative y

			for (int j = 0; j < _dim * _coeffCount; j++) {
				index = i % (_dim + _divCoeffCount);
				if (index == 0) { //x에 대해서
					if (j < _coeffCount) {
						divmls._nodes[tid]._B_Vec[i][j] = Basis[j];
					}
					else {
						divmls._nodes[tid]._B_Vec[i][j] = 0.0;
					}
				}
				else if (index == 1) { //y에 대해서
					if ((j >= _coeffCount) && (j < (2 * _coeffCount))) {
						divmls._nodes[tid]._B_Vec[i][j] = Basis[j - _coeffCount];
					}
					else {
						divmls._nodes[tid]._B_Vec[i][j] = 0.0;
					}
				}
				else { //div에 대해서
				 // divergence constraine
					divmls._nodes[tid]._B_Vec[i][j] = divBasis[j];
				}
			}

			for (int i = 0; i < (_dim + _divCoeffCount) * _constrainCount; i++) {
				for (int j = 0; j < _dim * _coeffCount; j++) {
					divmls._nodes[tid]._Bt_Vec[j][i] = divmls._nodes[tid]._B_Vec[i][j];
				}
			}
		}
	}
}

__global__ void precomputePVec(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		int _dim = divmls._dim;
		int _divCoeffCount = divmls._divCoeffCount;
		int _constrainCount = divmls._constrainCount;
		int _coeffCount = divmls._coeffCount;
		
		int index = 0;
		int id = 0;
		int cid = 0;

		for (int i = 0; i < (_dim + _divCoeffCount) * _constrainCount; i++) {
			index = i % (_dim + _divCoeffCount);
			id = i / (_dim + _divCoeffCount);
			cid = divmls._nodes[tid]._constrain[id];
			if (index == 0) {
				divmls._nodes[tid]._P_Vec[i] = divmls._velocityD[cid].x;
			}
			else if (index == 1) {
				divmls._nodes[tid]._P_Vec[i] = divmls._velocityD[cid].y;
			}
			else if (index == 2) {
				double div = 0.0;
				if (id == 0) {
					div = divmls._nodes[tid]._div[4];
				}
				if (id == 1) {
					div = divmls._nodes[tid]._div[5];
				}
				if (id == 2) {
					div = divmls._nodes[tid]._div[7];
				}
				if (id == 3) {
					div = divmls._nodes[tid]._div[8];
				}
				divmls._nodes[tid]._P_Vec[i] = div; 
			}
		}
	}
}

__global__ void computeWeightingMatD(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		int _dim = divmls._dim;
		int _divCoeffCount = divmls._divCoeffCount;
		int _constrainCount = divmls._constrainCount;
		int _coeffCount = divmls._coeffCount;

		double r;
		double h;
		h = 1.0 / (RES_SMOKE + 2);
		double3 btPos = divmls._btPosD[tid];
		
		for (int i = 0; i < (_dim + _divCoeffCount) * _constrainCount; i++) {
			int cid = divmls._nodes[tid]._constrain[i / (_dim + _divCoeffCount)];
			double3 cPos = divmls._posD[tid];
			double3 diff = make_double3(btPos.x - cPos.x, btPos.y-cPos.y, btPos.z-cPos.z); 
			r = fabs(sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z));

			for (int j = 0; j < (_dim + _divCoeffCount) * _constrainCount; j++) {
				if (i == j) {
					divmls._nodes[tid]._W_Vec[i][j] = 1.0 / ((r * r * r * r) + 0.00000001);
				}
				else {
					divmls._nodes[tid]._W_Vec[i][j] = 0.0;
				}
			}
		}
	}
}

__global__ void computePolyInterpD(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		int _dim = divmls._dim;
		int _divCoeffCount = divmls._divCoeffCount;
		int _constrainCount = divmls._constrainCount;
		int _coeffCount = divmls._coeffCount;
		int index = 0;

		for (int i = 0; i < _dim; i++) {
			double cx = divmls._btPosD[tid].x;
			double cy = divmls._btPosD[tid].y;
			double basis[6] = { 1.0, cx, cy, cx * cx, cy * cy, cx * cy };
			for (int j = 0; j < _dim * _coeffCount; j++) {
				index = i % _dim;
				if (index == 0) {
					if (j < _coeffCount) {
						divmls._nodes[tid]._B_Mat[i][j] = basis[j];
					}
					else {
						divmls._nodes[tid]._B_Mat[i][j] = 0.0;
					}
				}
				else if (index == 1) {
					if ((j >= _coeffCount) && (j < (2 * _coeffCount))) {
						divmls._nodes[tid]._B_Mat[i][j] = basis[j - _coeffCount];
					}
					else {
						divmls._nodes[tid]._B_Mat[i][j] = 0.0;
					}
				}
			}
		}
	}
}

__device__ double sign(double a, double b) {
	return (b >= 0.0 ? fabs(a) : -fabs(a));
}

__device__ double dmax(double a, double b)
{
	double dmaxarg1(a), dmaxarg2(b);
	return (dmaxarg1 > dmaxarg2 ? dmaxarg1 : dmaxarg2);
}

__device__ int imin(int a, int b)
{
	int iminarg1(a), iminarg2(b);
	return (iminarg1 < iminarg2 ? iminarg1 : iminarg2);
}

__device__ double dsqr(double a)
{
	double dsqrarg(a);
	return (dsqrarg == 0.0 ? 0.0 : dsqrarg * dsqrarg);
}

__device__ double dpythag(double a, double b)
{
	double absa(fabs(a)), absb(fabs(b));
	if (absa > absb) return absa * sqrt(1.0 + dsqr(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + dsqr(absa / absb)));
}

__device__ bool svd_decomp(double _A[], double _S[], double _V[]) {
	//특이값 분해 - KimJonghyun
}

__device__ void svd_backsub(double _A[], double _S[], double _V[], double _b[], double _x[]) {
	int     m(12), n(12);
	int     jj, j, i;
	double  s;
	double  tmp[12];

	// Calculate U^T * B
	for (j = 0; j < n; j++)
	{
		s = 0.0;
		if (_S[j])  // non-zero result only if S_j is nonzero
		{
			for (i = 0; i < m; i++) s += _A[i * 12 + j] * _b[i];
			s /= _S[j];
		}
		tmp[j] = s;
	}

	// Matrix multiply by V to get answer
	for (j = 0; j < n; j++)
	{
		s = 0.0;
		for (jj = 0; jj < n; jj++) s += _V[j*12+jj] * tmp[jj];
		_x[j] = s;
	}
}

__global__ void solveD(DIVMLS divmls) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		int _dim = divmls._dim;
		int _divCoeffCount = divmls._divCoeffCount;
		int _constrainCount = divmls._constrainCount;
		int _coeffCount = divmls._coeffCount;
		double value = 0.0;
		for (int i = 0; i < _dim * _coeffCount; i++) {
			for (int j = 0; j < (_dim + _divCoeffCount) * _constrainCount; j++) {
				divmls._nodes[tid]._W_Mat[i][j] = 0.0;
				for (int k = 0; k < (_dim + _divCoeffCount) * _constrainCount; k++) {
					divmls._nodes[tid]._W_Mat[i][j] += divmls._nodes[tid]._Bt_Vec[i][k] * divmls._nodes[tid]._W_Vec[k][j];
				}
			}
		}

		for (int i = 0; i < _dim * _coeffCount; i++) {
			for (int j = 0; j < _dim * _coeffCount; j++) {
				value = 0.0;
				for (int k = 0; k < (_dim + _divCoeffCount) * _constrainCount; k++) {
					value += divmls._nodes[tid]._W_Mat[i][k] * divmls._nodes[tid]._B_Vec[k][j];
				}
				divmls._svd[tid]._A[i*12+j] = value;
			}
		}

		for (int i = 0; i < _dim * _coeffCount; i++) {
			value = 0.0;
			for (int j = 0; j < (_dim + _divCoeffCount) * _constrainCount; j++) {
				value += divmls._nodes[tid]._W_Mat[i][j] * divmls._nodes[tid]._P_Vec[j];
			}
			divmls._svd[tid]._b[i] = value;
		}

		bool success = svd_decomp(divmls._svd[tid]._A, divmls._svd[tid]._S, divmls._svd[tid]._V);

		svd_backsub(divmls._svd[tid]._A, divmls._svd[tid]._S, divmls._svd[tid]._V, divmls._svd[tid]._b, divmls._svd[tid]._x);
		double v[2];
		for (int i = 0; i < _dim; i++) {
			v[i] = 0.0;
			for (int j = 0; j < _dim * _coeffCount; j++) {
				v[i] += divmls._nodes[tid]._B_Mat[i][j] * divmls._svd[tid]._x[j];
			}
		}
		//divmls._resultD[tid] = make_double3(divmls._svd[tid]._x[0], divmls._svd[tid]._x[1], 0.0);
		divmls._resultD[tid] = make_double3(v[0], v[1], 0.0);
	}

}

__device__ double AngleBetweenVectors(double3 v1, double3 v2) {
	double dot, cross;
	double3 tmp;
	dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	tmp = make_double3(v1.y*v2.z-v1.z*v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
	cross = sqrt(tmp.x*tmp.x+tmp.y*tmp.y+tmp.z*tmp.z);
	return abs(atan2(cross, dot) * (180/3.14));
}

__global__ void computeVelocity(DIVMLS divmls)
{
	//원본 데이터와의 각도 계산
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < divmls._size) {
		double3 original = divmls._velocityD[tid];
		double3 result = divmls._resultD[tid];
		double angle = AngleBetweenVectors(original, result);
		if (angle > 120) {
			divmls._resultD[tid] = make_double3(0.0, 0.0, 0.0);
		}
	}
}

void DIVMLSKernel::simulation(DIVMLS divmls) {
	int numThreads, numBlocks;
	computeGridSize(divmls._size, 256, numBlocks, numThreads);
	printf("numT : %d, numB : %d\n", numThreads, numBlocks);

	//1.역추적
	backTracingD <<<numBlocks, numThreads >>> (divmls);
	//2.발산계산
	computeDivD << <numBlocks, numThreads >> > (divmls);
	//3.MLS 보간을 위한 데이터 전처리
	precomputeBVec << <numBlocks, numThreads >>> (divmls);
	precomputePVec << <numBlocks, numThreads >>> (divmls);
	//4.거리 기반 가중치
	computeWeightingMatD << <numBlocks, numThreads >> > (divmls);
	//5.MLS 보간을 위한 행렬 구성
	computePolyInterpD << <numBlocks, numThreads >> > (divmls);
	
	//6.SVD 특이값 분해로 MLS 기반 보간된 속도 벡터 계산 (by.KimJongHyun)
	solveD << <numBlocks, numThreads >> > (divmls);
	
	//7.보간된 속도 유효성 검사
	computeVelocity << <numBlocks, numThreads >> > (divmls);
}