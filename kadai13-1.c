#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define DX 0.1
#define DT 0.002
#define D 1.0
#define STEPMAX 10000
#define INTV 100

double u[N+2] = {0.0};

void writeVtk(int count) {
	FILE *fp;
	char fname[256];
	double uu[N+2][150] = {0};
	int i, j;

	for (i=0; i<N+2; i++) {
		j = (int)(u[i]/0.01)+10;
		uu[i][j] = 1.0;
	}

	sprintf(fname, "data_%08d.vtk", count);
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 4.1 \n");
	fprintf(fp, "REACTION_DIFFUSION_1D\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS \n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", N+2, 150);
	fprintf(fp, "ORIGIN 0 0 0\n");
	fprintf(fp, "SPACING %f %f 0 \n", DX, DX);
	fprintf(fp, "POINT_DATA %d\n", (N+2)*150);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(j=0; j<150; j++) {
		for(i=0; i<N+2; i++) {
			fprintf(fp, "%f\n", uu[i][j]);
		}
	}

	fclose(fp);
}

double func_u(double x) {
	return x*(1.0-x);
}

void init() {
	int i;
	for (i = 0; i < N + 2; i++) {
		if (i == N / 2) {
			u[i] = 0.1; // 中心のみ u = 0.1
		} else {
			u[i] = 0.0; // 他は u = 0.0
		}
	}
	return;
}

void pbc() {
	u[0] = u[N];       // 左端を右端に合わせる
	u[N+1] = u[1];     // 右端を左端に合わせる
	return;
}

int main() {
	int i;
	int step;
	double u0[N+2];

	// Initialize the system
	init();

	// Output the inital system
	step = 0;
	printf("%d\n", step);
	writeVtk(step);

	// Start main loop
	for (step=1; step<=STEPMAX; step++) {
	
		// Save current u as u0
		for (i=0; i<N+2; i++) {
			u0[i] = u[i];
		}

		// Update u
		for (i=1; i<N+1; i++) {
			u[i] = u0[i] + D/(DX*DX)*(u0[i+1]+u0[i-1]-2.0*u0[i])*DT + func_u(u0[i])*DT;
		}

		// Periodic boundary conditions
		pbc();

		// Output
		printf("%d\n", step);
		if ( step%INTV == 0 ) writeVtk(step);
	}

	return 0;
}
