#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 256
#define DX 0.01
#define DT 0.005
#define DU 0.0005
#define DV 0.000005
#define STEPMAX 10000
#define INTV 100

double u[N+2][N+2] = {0.0};
double v[N+2][N+2] = {0.0};

// Random number generator [0,1)
double urand(void) {
	return rand()/(RAND_MAX+1.0);
}

void writeVtk(int count) {
	FILE *fp;
	char fname[256];
	int i, j;

	sprintf(fname, "data_%08d.vtk", count);
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 4.1 \n");
	fprintf(fp, "REACTION_DIFFUSION_2D\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS \n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", N+2, N+2);
	fprintf(fp, "ORIGIN 0 0 0\n");
	fprintf(fp, "SPACING %f %f 0 \n", DX, DX);
	fprintf(fp, "POINT_DATA %d\n", (N+2)*(N+2));
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(j=0; j<N+2; j++) {
		for(i=0; i<N+2; i++) {
			fprintf(fp, "%f\n", u[i][j]);
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
	int i, j;
	double u_initial = 0.13651;

	for (i = 0; i < N+2; i++) {
		for (j = 0; j < N+2; j++) {
			// ‰Šú’l‚ÉƒmƒCƒY‚ð‰Á‚¦‚é
			u[i][j] = u_initial + (2.0 * urand() - 1.0) * 0.01;  // }1% ƒmƒCƒY
			v[i][j] = u_initial + (2.0 * urand() - 1.0) * 0.01;  // }1% ƒmƒCƒY
		}
	}
	return;
}

void pbc() {
	int i, j;

	// ¶‰E‚Ì‹«ŠE
	for (j = 0; j < N+2; j++) {
		u[0][j] = u[N][j];
		u[N+1][j] = u[1][j];
		v[0][j] = v[N][j];
		v[N+1][j] = v[1][j];
	}

	// ã‰º‚Ì‹«ŠE
	for (i = 0; i < N+2; i++) {
		u[i][0] = u[i][N];
		u[i][N+1] = u[i][1];
		v[i][0] = v[i][N];
		v[i][N+1] = v[i][1];
	}

	return;
}

int main() {
	int i, j, step;
	double uu[N+2][N+2];
	double vv[N+2][N+2];

	// Initialize the system
	init();

	// Apply periodic boundary conditions
	pbc();

	// Output the initial system
	step = 0;
	printf("%d\n", step);
	writeVtk(step);

	// Start main loop
	for (step=1; step<=STEPMAX; step++) {
		// Save the current u[:][:] and v[:][:] as uu[:][:] and vv[:][:]
		for (i=0; i<N+2; i++) {
			for (j=0; j<N+2; j++) {
				uu[i][j] = u[i][j];
				vv[i][j] = v[i][j];
			}
		}

		// Update u and v
		for (i=1; i<N+1; i++) {
			for (j=1; j<N+1; j++) {
				double du_diff = DU/(DX*DX) * (uu[i+1][j] + uu[i-1][j] + uu[i][j+1] + uu[i][j-1] - 4.0*uu[i][j]);
				double dv_diff = DV/(DX*DX) * (vv[i+1][j] + vv[i-1][j] + vv[i][j+1] + vv[i][j-1] - 4.0*vv[i][j]);
				u[i][j] = uu[i][j] + du_diff*DT + func_u(uu[i][j], vv[i][j])*DT;
				v[i][j] = vv[i][j] + dv_diff*DT + func_v(uu[i][j], vv[i][j])*DT;
			}
		}

		// Periodic boundary conditions
		pbc();

		// Output current state
		printf("%d\n", step);
		if (step % INTV == 0) writeVtk(step);
	}

	return 0;
}
