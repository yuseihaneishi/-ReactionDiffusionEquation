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
#define ALPHA 0.035
#define BETA 0.065

double u[N+2][N+2] = {0.0};
double v[N+2][N+2] = {0.0};

void writeVtk(int count) {
    FILE *fp;
    char fname[256];
    int i, j;

    sprintf(fname, "data_%08d.vtk", count);
    fp = fopen(fname, "w");
    fprintf(fp, "# vtk DataFile Version 4.1\n");
    fprintf(fp, "GRAY_SCOTT\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", N+2, N+2);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %f %f 0\n", DX, DX);
    fprintf(fp, "POINT_DATA %d\n", (N+2)*(N+2));
    fprintf(fp, "SCALARS U double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(j = 0; j < N + 2; j++) {
        for(i = 0; i < N + 2; i++) {
            fprintf(fp, "%f\n", u[i][j]);
        }
    }

    fclose(fp);
}

double urand(void) {
    return rand() / (RAND_MAX + 1.0);
}

double func_u(double u, double v) {
    return -u * v * v + ALPHA * (1.0 - u);
}

double func_v(double u, double v) {
    return u * v * v - (ALPHA + BETA) * v;
}

void init() {
    int i, j;
    for (i = 0; i < N + 2; i++) {
        for (j = 0; j < N + 2; j++) {
            u[i][j] = 1.0 + 0.1 * urand();
            v[i][j] = 0.0 + 0.1 * urand();
        }
    }
    for (i = N / 2 - 10; i <= N / 2 + 10; i++) {
        for (j = N/2 - 10; j <= N/2 + 10; j++) {
            u[i][j] = 0.5;
            v[i][j] = 0.25;
        }
    }
}

void pbc() {
    int i;
    for (i = 1; i < N + 1; i++) {
        u[i][0] = u[i][N];
        u[i][N+1] = u[i][1];
        u[0][i] = u[N][i];
        u[N+1][i] = u[1][i];

        v[i][0] = v[i][N];
        v[i][N+1] = v[i][1];
        v[0][i] = v[N][i];
        v[N+1][i] = v[1][i];
    }
}

int main() {
    int i, j, step;
    double uu[N+2][N+2], vv[N+2][N+2];

    init();
    pbc();

    step = 0;
    printf("%d\n", step);
    writeVtk(step);

    for (step = 1; step <= STEPMAX; step++) {
        for (i = 0; i < N+2; i++) {
            for (j = 0; j < N+2; j++) {
                uu[i][j] = u[i][j];
                vv[i][j] = v[i][j];
            }
        }

        for (i = 1; i < N+1; i++) {
            for (j = 1; j < N+1; j++) {
                u[i][j] = uu[i][j] + DU * (uu[i+1][j] + uu[i-1][j] + uu[i][j+1] + uu[i][j-1] - 4.0 * uu[i][j]) * DT / (DX * DX)
                          + func_u(uu[i][j], vv[i][j]) * DT;
                v[i][j] = vv[i][j] + DV * (vv[i+1][j] + vv[i-1][j] + vv[i][j+1] + vv[i][j-1] - 4.0 * vv[i][j]) * DT / (DX * DX)
                          + func_v(uu[i][j], vv[i][j]) * DT;
            }
        }

        pbc();

        if (step % INTV == 0) {
            printf("%d\n", step);
            writeVtk(step);
        }
    }
    return 0;
}