#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define DX 0.01
#define DT 0.005
#define DU 0.0002
#define DV 0.0001
#define STEPMAX 10000
#define INTV 100

double u[N+2] = {0.0};
double v[N+2] = {0.0};

void writeVtk(int count) {
	FILE *fp;
	char fname[256];
	double uv[N+2][256] = {0};
	int i, j;

	for (i=0; i<N+2; i++) {
		j = (int)(u[i]/0.005)+10;
		uv[i][j] = 1.0;
		j = (int)(v[i]/0.005)+10;
		uv[i][j] = -1.0;
	}
	uv[0][255] = 1.0;
	uv[N+1][255] = -1.0;

	sprintf(fname, "data_%08d.vtk", count);
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 4.1 \n");
	fprintf(fp, "REACTION_DIFFUSION_1D_2C\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS \n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", N+2, 256);
	fprintf(fp, "ORIGIN 0 0 0\n");
	fprintf(fp, "SPACING %f %f 0 \n", DX, DX/5);
	fprintf(fp, "POINT_DATA %d\n", (N+2)*256);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(j=0; j<256; j++) {
		for(i=0; i<N+2; i++) {
			fprintf(fp, "%f\n", uv[i][j]);
		}
	}

	fclose(fp);
}

double func_u(double x, double y) {
	double a = 0.1;
	double b = 1.0;
	double c = 0.01;
	return 1.0/a*( x*(1.0-x) - b*y*(x-c)/(x+c) );
}

double func_v(double x, double y) {
	return x-y;
}

void init() {
	int i;
	for (i = 0; i < N+2; i++) {
		if (i == N/2) {
			u[i] = 0.1;  // 初期値 u を中心に設定
			v[i] = 0.2;  // 初期値 v を中心に設定
		} else {
			u[i] = 0.0;
			v[i] = 0.0;
		}
	}
	return;
}

void pbc() {
	u[0] = u[N];
	u[N+1] = u[1];
	v[0] = v[N];
	v[N+1] = v[1];
	return;
}

int main() {
	int i;
	int step;
	double u0[N+2];
	double v0[N+2];

	// Initialize the system
	init();

	// Output the initial system
	step = 0;
	printf("%d\n", step);
	writeVtk(step);

	// Start main loop
	for (step=1; step<=STEPMAX; step++) {
		// Save the current u[:] and v[:] values as u0[:] and v0[:], respectively
		for (i=0; i<N+2; i++) {
			u0[i] = u[i];
			v0[i] = v[i];
		}

		// Update u[:] and v[:]
		for (i=1; i<N+1; i++) {
			u[i] = u0[i] + DU/(DX*DX)*(u0[i+1] + u0[i-1] - 2.0*u0[i])*DT + func_u(u0[i], v0[i])*DT;
			v[i] = v0[i] + DV/(DX*DX)*(v0[i+1] + v0[i-1] - 2.0*v0[i])*DT + func_v(u0[i], v0[i])*DT;
		}

		// Periodic boundary conditions
		pbc();

		// Output current state
		printf("%d\n", step);
		if (step % INTV == 0) writeVtk(step);
	}

	return 0;
}
