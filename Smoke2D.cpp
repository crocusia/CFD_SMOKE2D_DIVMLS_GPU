#include "Smoke2D.h"

void Smoke2D::Simulation(void)
{
	double t = 0.0;
	dt = 0.0;
	int iter = 0;
	while (t < dt2) {
		dt = CFL();
		if (t + dt > dt2) {
			dt = dt2 - t;
		}
		//밀도 및 속도 설정
		Sourcing();
		//속도 업데이트
		UpdateVelocity();
		//밀도 업데이트
		UpdateDensity();
		
		iter++;
		t += dt;
	}
	m_Frame++;
}

//가우스-자이델 기반 확산 연산
void Smoke2D::Solver(double *x, double *x0, int alpha, int beta, int type, int times)
{
	for (int k = 0; k < times; k++) {
		for (int i = 1; i <= RES_SMOKE; i++) {
			for (int j = 1; j <= RES_SMOKE; j++) {
				x[IX2(i, j)] = (x0[IX2(i, j)] + alpha * (x[IX2(i - 1, j)] + x[IX2(i + 1, j)] + x[IX2(i, j - 1)] + x[IX2(i, j + 1)])) / beta;
			}
		}
		SetBoundary(type, x);
	}
}

//발산 제거를 위한 압력 연산
void Smoke2D::Project(double *ux, double *uy, double *p, double *div)
{
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			div[IX2(i, j)] = -0.5 * (ux[IX2(i + 1, j)] - ux[IX2(i - 1, j)] + uy[IX2(i, j + 1)] - uy[IX2(i, j - 1)]) / RES_SMOKE;
			p[IX2(i, j)] = 0;
		}
	}
	SetBoundary(0, div);
	SetBoundary(0, p);

	Solver(p, div, 1.0, 4.0, 0, 20);

	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			ux[IX2(i, j)] -= 0.5 * RES_SMOKE * (p[IX2(i + 1, j)] - p[IX2(i - 1, j)]);
			uy[IX2(i, j)] -= 0.5 * RES_SMOKE * (p[IX2(i, j + 1)] - p[IX2(i, j - 1)]);
		}
	}
	SetBoundary(1, ux);
	SetBoundary(2, uy);
}

//이류
void Smoke2D::Advect(double *d, double *d0, double *ux, double *uy, int bType)
{
	int i0, j0, i1, j1;
	double x, y, s0, t0, s1, t1, h0;

	h0 = dt * RES_SMOKE;
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			x = i - h0 * ux[IX2(i, j)];
			y = j - h0 * uy[IX2(i, j)]; 

			if (x < 0.5)				x = 0.5;
			if (x > RES_SMOKE + 0.5)	x = RES_SMOKE + 0.5;
			i0 = (int)x; i1 = i0 + 1; 
			if (y < 0.5)				y = 0.5;
			if (y > RES_SMOKE + 0.5)	y = RES_SMOKE + 0.5;
			j0 = (int)y; j1 = j0 + 1;

			s1 = x - i0; s0 = 1 - s1; 
			t1 = y - j0; t0 = 1 - t1;
		
			d[IX2(i, j)] = s0 * (t0 * d0[IX2(i0, j0)] + t1 * d0[IX2(i0, j1)]) + s1 * (t0 * d0[IX2(i1, j0)] + t1 * d0[IX2(i1, j1)]);
		}
	}
	SetBoundary(bType, d);
}

//GPU 기반 DCMLS 속도장 연산
void Smoke2D::advectMLS(double* dx, double* dy, double* ux, double* uy) {
	_divMLS.zeroMLS();
	_divMLS.copyHtoD(ux, uy);
	_kernel.simulation(_divMLS);
	_divMLS.copyDtoH(dx, dy);
}

void Smoke2D::UpdateVelocity(void)
{
	double halfrdx = 0.5 * RES_SMOKE;
	int size = (RES_SMOKE + 2) * (RES_SMOKE + 2);
	double a = dt * VISC * RES_SMOKE* RES_SMOKE;
	
	//x축 속도 확산
	SWAP(u0, u);
	Solver(u, u0, a, 1.0 + 4.0 * a, 1, 20);
	//y축 속도 확산
	SWAP(v0, v);
	Solver(v, v0, a, 1.0 + 4.0 * a, 2, 20);
	
	//확산된 속도장의 발산 제거
	Project(u, v, u0, v0);
	
	SWAP(u0, u);
	SWAP(v0, v);

	// 1-1.이중선형보간법 기반 속도장 이류
	Advect(u, u0, u0, v0, 1);
	Advect(v, v0, u0, v0, 2);

	// 1-2.DCMLS보간법 기반 속도장 이류 (by. GPU)
	advectMLS(mu, mv, u0, v0);

	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			int index = i + (RES_SMOKE + 2) * j;

			Vec3<double> uv(u[index], v[index], 0.0);
			auto uv_mag = uv.Length();

			Vec3<double> mls_uv(mu[index], mv[index], 0.0);
			mls_uv = mls_uv * dt * d[index] * ALPHA * 10.0;

			//DCMLS 기반 와류 가중치 계산
			Vec3<double> uv_unit = uv; uv_unit.Normalize();
			Vec3<double> mls_unit = mls_uv; mls_unit.Normalize();
			//3.두 속도 벡터의 각도 차(유사도) 연산
			mls_vort[index] = uv_unit.Dot(mls_unit) * 0.5 + 0.5;

			//2-1.DCMLS 속도를 외력으로 적용
			uv += mls_uv;
			//2-2.정규화로 크기는 유지하고 방향만 변환
			uv.Normalize();
			uv *= uv_mag;

			//속도 업데이트
			u[index] = uv.x();
			v[index] = uv.y();
		}
	}

	//와류 크기 연산
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			int index = i + (RES_SMOKE + 2) * j;
			vort[IX2(i, j)] = halfrdx * ((v[IX2(i + 1, j)] - v[IX2(i - 1, j)]) - (u[IX2(i, j + 1)] - u[IX2(i, j - 1)]));
			vort_mag[IX2(i, j)] = fabs(vort[IX2(i, j)]);
		}
	}
	SetBoundary(0, vort);
	SetBoundary(0, vort_mag);

	//와류 기울기 (변화율) 연산 및 정규화 by.jhkim
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			vort_u[IX2(i, j)] = halfrdx * (vort_mag[IX2(i + 1, j)] - vort_mag[IX2(i - 1, j)]);
			vort_v[IX2(i, j)] = halfrdx * (vort_mag[IX2(i, j + 1)] - vort_mag[IX2(i, j - 1)]);
			double len = sqrt(vort_u[IX2(i, j)] * vort_u[IX2(i, j)] + vort_v[IX2(i, j)] * vort_v[IX2(i, j)]);
			if (len < VORTICITY_EPS) {
				vort_u[IX2(i, j)] = 0.0;
				vort_v[IX2(i, j)] = 0.0;
			}
			else {
				vort_u[IX2(i, j)] /= len;
				vort_v[IX2(i, j)] /= len;
			}
		}
	}
	SetBoundary(0, vort_u);
	SetBoundary(0, vort_v);
	
	//4.mls_vort(Stable Fluid 속도 벡터와 DCMLS 속도 벡터의 유사도) 기반 와류 강도 조절
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			int index = i + (RES_SMOKE + 2) * j;
			auto w = mls_vort[index] * 0.05;
			u[IX2(i, j)] += dt * (VORTICITY + w) * (vort_v[IX2(i, j)] * vort[IX2(i, j)]);
			v[IX2(i, j)] += dt * (VORTICITY + w) * (-vort_u[IX2(i, j)] * vort[IX2(i, j)]);
		}
	}

	SetBoundary(1, u);
	SetBoundary(2, v);

	for (int i = 0; i < size; i++) {
		mv[i] = mu[i] = 0.0;
	}

	Project(u, v, u0, v0);
}

//업데이트 된 속도장에 따라 밀도 이류
void Smoke2D::UpdateDensity(void)
{
	if (DIFF > VORTICITY_EPS) {
		double a = dt * DIFF * RES_SMOKE * RES_SMOKE;
		SWAP(d0, d);
		Solver(d, d0, a, 1.0 + 4.0 * a, 1, 20);
	}
	SWAP(d0, d);
	Advect(d, d0, u, v, 0);
}
