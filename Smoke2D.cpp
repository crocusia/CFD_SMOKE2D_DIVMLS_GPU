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
		//�е� �� �ӵ� ����
		Sourcing();
		//�ӵ� ������Ʈ
		UpdateVelocity();
		//�е� ������Ʈ
		UpdateDensity();
		
		iter++;
		t += dt;
	}
	m_Frame++;
}

//���콺-���̵� ��� Ȯ�� ����
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

//�߻� ���Ÿ� ���� �з� ����
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

//�̷�
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

//GPU ��� DCMLS �ӵ��� ����
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
	
	//x�� �ӵ� Ȯ��
	SWAP(u0, u);
	Solver(u, u0, a, 1.0 + 4.0 * a, 1, 20);
	//y�� �ӵ� Ȯ��
	SWAP(v0, v);
	Solver(v, v0, a, 1.0 + 4.0 * a, 2, 20);
	
	//Ȯ��� �ӵ����� �߻� ����
	Project(u, v, u0, v0);
	
	SWAP(u0, u);
	SWAP(v0, v);

	// 1-1.���߼��������� ��� �ӵ��� �̷�
	Advect(u, u0, u0, v0, 1);
	Advect(v, v0, u0, v0, 2);

	// 1-2.DCMLS������ ��� �ӵ��� �̷� (by. GPU)
	advectMLS(mu, mv, u0, v0);

	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			int index = i + (RES_SMOKE + 2) * j;

			Vec3<double> uv(u[index], v[index], 0.0);
			auto uv_mag = uv.Length();

			Vec3<double> mls_uv(mu[index], mv[index], 0.0);
			mls_uv = mls_uv * dt * d[index] * ALPHA * 10.0;

			//DCMLS ��� �ͷ� ����ġ ���
			Vec3<double> uv_unit = uv; uv_unit.Normalize();
			Vec3<double> mls_unit = mls_uv; mls_unit.Normalize();
			//3.�� �ӵ� ������ ���� ��(���絵) ����
			mls_vort[index] = uv_unit.Dot(mls_unit) * 0.5 + 0.5;

			//2-1.DCMLS �ӵ��� �ܷ����� ����
			uv += mls_uv;
			//2-2.����ȭ�� ũ��� �����ϰ� ���⸸ ��ȯ
			uv.Normalize();
			uv *= uv_mag;

			//�ӵ� ������Ʈ
			u[index] = uv.x();
			v[index] = uv.y();
		}
	}

	//�ͷ� ũ�� ����
	for (int i = 1; i <= RES_SMOKE; i++) {
		for (int j = 1; j <= RES_SMOKE; j++) {
			int index = i + (RES_SMOKE + 2) * j;
			vort[IX2(i, j)] = halfrdx * ((v[IX2(i + 1, j)] - v[IX2(i - 1, j)]) - (u[IX2(i, j + 1)] - u[IX2(i, j - 1)]));
			vort_mag[IX2(i, j)] = fabs(vort[IX2(i, j)]);
		}
	}
	SetBoundary(0, vort);
	SetBoundary(0, vort_mag);

	//�ͷ� ���� (��ȭ��) ���� �� ����ȭ by.jhkim
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
	
	//4.mls_vort(Stable Fluid �ӵ� ���Ϳ� DCMLS �ӵ� ������ ���絵) ��� �ͷ� ���� ����
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

//������Ʈ �� �ӵ��忡 ���� �е� �̷�
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
