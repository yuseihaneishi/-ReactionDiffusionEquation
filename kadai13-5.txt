#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 256
#define DX 0.01
#define DT 1.0
#define DU 2.0e-5
#define DV 1.0e-5
#define STEPMAX 10000
#define INTV 100
#define ALPHA 0.022
#define BETA 0.051

double u[N+2][N+2] = {0.0};
double v[N+2][N+2] = {0.0};

void writeVtk(int count) {
	FILE *fp;
	char fname[256];
	int i, j;

	sprintf(fname, "data_%08d.vtk", count);
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 4.1 \n");
	fprintf(fp, "GRAY_SCOTT\n");
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

// Differential equation part for u in Eq. (9)
double func_u(double u_val, double v_val) {
	return -u_val * v_val * v_val + ALPHA * (1.0 - u_val);
}

// Differential equation part for v in Eq. (10)
double func_v(double u_val, double v_val) {
	return u_val * v_val * v_val - (ALPHA + BETA) * v_val;
}

// Random number [0,1) generator
double urand(void) {
	return rand()/(RAND_MAX+1.0);
}

// Initialize the system
void init() {
	int i, j;

	// Set the entire system to initial values
	for (i = 0; i < N+2; i++) {
		for (j = 0; j < N+2; j++) {
			u[i][j] = 1.0;
			v[i][j] = 0.0;
		}
	}

	// Set the center 20x20 grid to a different state
	for (i = N/2 - 10; i < N/2 + 10; i++) {
		for (j = N/2 - 10; j < N/2 + 10; j++) {
			u[i][j] = 0.5 + (urand() - 0.5) * 0.1;  // Add noise }0.05
			v[i][j] = 0.25 + (urand() - 0.5) * 0.1; // Add noise }0.05
		}
	}

	return;
}

// Periodic boundary conditions
void pbc() {
	int i, j;

	// Left-right boundaries
	for (j = 0; j < N+2; j++) {
		u[0][j] = u[N][j];
		u[N+1][j] = u[1][j];
		v[0][j] = v[N][j];
		v[N+1][j] = v[1][j];
	}

	// Top-bottom boundaries
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

	// Output the initial state
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
				u[i][j] = uu[i][j] + du_diff * DT + func_u(uu[i][j], vv[i][j]) * DT;
				v[i][j] = vv[i][j] + dv_diff * DT + func_v(uu[i][j], vv[i][j]) * DT;
			}
		}

		// Apply periodic boundary conditions
		pbc();

		// Output
		printf("%d\n", step);
		if (step % INTV == 0) writeVtk(step);
	}

	return 0;
}
