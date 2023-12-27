#include "pso.h"
#include <math.h>

Matrix createMatrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double*)malloc(cols * sizeof(double));
    }
    return mat;
}

void freeMatrix(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        free(mat.data[i]);
    }
    free(mat.data);
}

void printMatrix(Matrix mat){
    int rows = mat.rows;
    int cols = mat.cols;
    for (int i=0; i<rows; i++)
    {
        printf("[");
        for (int j=0; j<cols; j++)
        {
            printf("  %.6f", mat.data[i][j]);
        }
        printf("]\n");
    }
    return;
}

void printMatrix(const char* str, Matrix mat)
{
    printf("%s:\n",str);
    printMatrix(mat);
    return;
}

Matrix get_row(Matrix mat, int row_idx)
{
    Matrix res = createMatrix(1, mat.cols);
    for (int i=0; i<mat.cols; i++)
    {
        res.data[0][i] = mat.data[row_idx][i]; 
    }
    return res;
}

Matrix matrixAdd(Matrix mat1, Matrix mat2) {
    Matrix result = createMatrix(mat1.rows, mat1.cols);
    for (int i = 0; i < mat1.rows; i++) {
        for (int j = 0; j < mat1.cols; j++) {
            result.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
        }
    }
    return result;
}

Matrix matrixSubtract(Matrix mat1, Matrix mat2) {
    Matrix result = createMatrix(mat1.rows, mat1.cols);
    for (int i = 0; i < mat1.rows; i++) {
        for (int j = 0; j < mat1.cols; j++) {
            result.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
        }
    }
    return result;
}

Matrix matrixScalarMultiply(Matrix mat, double scalar) {
    Matrix result = createMatrix(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            result.data[i][j] = mat.data[i][j] * scalar;
        }
    }
    return result;
}

Matrix matrixElementWiseMultiply(Matrix mat1, Matrix mat2) {
    Matrix result = createMatrix(mat1.rows, mat1.cols);
    for (int i = 0; i < mat1.rows; i++) {
        for (int j = 0; j < mat1.cols; j++) {
            result.data[i][j] = mat1.data[i][j] * mat2.data[i][j];
        }
    }
    return result;
}

Matrix matrixCopy(Matrix mat) {
    Matrix result = createMatrix(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            result.data[i][j] = mat.data[i][j];
        }
    }
    return result;
}

Matrix matrixTranspose(Matrix mat) {
    Matrix result = createMatrix(mat.cols, mat.rows);
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            result.data[j][i] = mat.data[i][j];
        }
    }
    return result;
}

Matrix matrixNullaryExpr(int rows, int cols, double (*func)()) {
    Matrix result = createMatrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.data[i][j] = func();
        }
    }
    return result;
}

double randomDouble() {
    return (double)rand() / RAND_MAX;
}

Matrix vectorLinSpaced(int size, double start, double end) {
    Matrix result = createMatrix(1, size);
    double step = (end - start) / (size - 1);
    double value = start;
    for (int i = 0; i < size; i++) {
        result.data[0][i] = value;
        value += step;
    }
    return result;
}

double fitnessFunction(Matrix position) {
    // TODO: Implement your fitness function here
    return 0.0;
}


void CLPSO(int Dimension, int Max_Gen, int Particle_Number, Matrix Rmin, Matrix Rmax, Matrix& gbest, double& gbestval) {
    int ps = Particle_Number;
    int D = Dimension;
    int me = Max_Gen;

    if (Rmin.cols == 1) {
        Matrix temp1 = matrixCopy(Rmin);
        Matrix temp2 = matrixCopy(Rmax);
        Rmin = matrixTranspose(matrixNullaryExpr(D, D, randomDouble));
        Rmin = matrixNullaryExpr(ps, D, randomDouble);
        Rmax = matrixNullaryExpr(ps, D, randomDouble);
        freeMatrix(temp1);
        freeMatrix(temp2);
    }

    Matrix mv = matrixScalarMultiply(matrixSubtract(Rmax, Rmin), 0.2);
    printMatrix("mv = ", mv);
    // 

    // Rmin_M, Rmax_M, Vmin_M, Vmax_M
    /*
    Matrix Rmin_M = matrixCopy(Rmin);
    Matrix Rmax_M = matrixCopy(Rmax);
    Matrix Vmin_M = matrixScalarMultiply(mv, -1.0);
    Matrix Vmax_M = matrixCopy(Vmin_M);
    */
    Matrix Rmin_M = createMatrix(ps, Rmin.cols);
    Matrix Rmax_M = createMatrix(ps, Rmin.cols);
    Matrix Vmin_M = createMatrix(ps, mv.cols);
    Matrix Vmax_M = createMatrix(ps, mv.cols);
    for (int i=0; i<ps; i++)
    {
        for(int j=0; j<Rmin.cols; j++)
        {
            Rmin_M.data[i][j] = Rmin.data[0][j];
            Rmax_M.data[i][j] = Rmax.data[0][j];
        }
        for(int k; k<mv.cols; k++)
        {
            Vmin_M.data[i][k] = -mv.data[0][k];
            Vmax_M.data[i][k] = mv.data[0][k];
        }
    }

    Matrix w = matrixSubtract(vectorLinSpaced(me, 0.9, 0.9), matrixScalarMultiply(vectorLinSpaced(me, 1.0, me), (0.7 / me)));
    double c1 = 0.8;
    double c2 = 1.49;

    Matrix pos = matrixAdd(Rmin_M, matrixElementWiseMultiply(matrixSubtract(Rmax_M, Rmin_M), matrixNullaryExpr(ps, D, randomDouble)));
    Matrix vel = matrixAdd(Vmin_M, matrixElementWiseMultiply(matrixSubtract(Vmax_M, Vmin_M), matrixNullaryExpr(ps, D, randomDouble)));

    Particle* particles = (Particle*)malloc(ps * sizeof(Particle));
    
    // 指针赋值操作, 不是值复制, 是地址复制, position操作数组改变时,pos.data同时改变, 应该改成值复制
    /*
    for (int i = 0; i < ps; i++) {
        particles[i].position = pos.data[i];
        particles[i].velocity = vel.data[i];
        particles[i].fitness = fitnessFunction(matrixCopy(pos));
        particles[i].pBest = particles[i].position;
        particles[i].pBestFitness = particles[i].fitness;
    }
    */
    //-----------------------------------------------------------//
    for (int i=0; i<ps; i++){
        particles[i].position = get_row(pos, i);
        particles[i].velocity = get_row(vel, i);
        particles[i].fitness = fitnessFunction(matrixCopy(pos));
        particles[i].pBest    = particles[i].position;
        particles[i].pBestFitness = particles[i].fitness;
    } 
    //-----------------------------------------------------------//

    gbestval = particles[0].pBestFitness;
    int minIndex = 0;
    for (int i = 1; i < ps; i++) {
        if (particles[i].pBestFitness < gbestval) {
            gbestval = particles[i].pBestFitness;
            minIndex = i;
        }
    }
    /* 这里好像不用转置, 转置后变成 4*1 的向量了，后面第252行的减法会出现问题
    gbest = matrixCopy(matrixTranspose(matrixNullaryExpr(1, D, randomDouble)));
    */
    gbest = matrixNullaryExpr(1, D, randomDouble);

    /* PBest 改成了Matrix, 代码同时更新
    for (int i = 0; i < D; i++) {
        gbest.data[0][i] = particles[minIndex].pBest[i];
    }
    */
    for (int i = 0; i < D; i++)
    {
        gbest.data[0][i] = particles[minIndex].pBest.data[0][i];
    }
    // 到这还没有问题，问题出在最后一个for 循环中
 
    for (int i = 1; i < me; i++) 
    {
        for (int k = 0; k < ps; k++) {
            Matrix r1 = matrixNullaryExpr(1, D, randomDouble);
            Matrix r2 = matrixNullaryExpr(1, D, randomDouble);

            Matrix temp1 = matrixSubtract(matrixCopy(particles[k].pBest), matrixCopy(particles[k].position));
            /* gbest 改成了引用类型， 不在使用*gbest索引, 因为都是Matrix 类型， 感觉可以不用copy函数
            Matrix temp2 = matrixSubtract(matrixCopy(*gbest), matrixCopy(particles[k].position));
            */
            Matrix temp2 = matrixSubtract(gbest, particles[k].position);

            particles[k].velocity = matrixAdd(matrixAdd(matrixScalarMultiply(particles[k].velocity, w.data[0][i]),
                matrixScalarMultiply(temp1, c1)), matrixScalarMultiply(temp2, c2));

            for (int j = 0; j < D; j++) {
                /* velocity 现在是Matrix，
                particles[k].velocity[j] = fmax(fmin(particles[k].velocity[j], mv.data[0][j]), -mv.data[0][j]);
                */
                particles[k].velocity.data[0][j] = fmax(fmin(particles[k].velocity.data[0][j], mv.data[0][j]), -mv.data[0][j]);
            }

            particles[k].position = matrixAdd(particles[k].position, particles[k].velocity);

            for (int j = 0; j < D; j++) {
                /* position 现在是Matrix， 
                particles[k].position[j] = fmax(fmin(particles[k].position[j], Rmax_M.data[k][j]), Rmin_M.data[k][j]);
                */
                particles[k].position.data[0][j] = fmax(fmin(particles[k].position.data[0][j], Rmax_M.data[k][j]), Rmin_M.data[k][j]);
            }
            particles[k].fitness = fitnessFunction(matrixCopy(particles[k].position));

            if (particles[k].fitness <= particles[k].pBestFitness) {
                for (int j = 0; j < D; j++) {
                    /* pBest 现在是Matrix,
                    particles[k].pBest[j] = particles[k].position[j];
                    */
                    particles[k].pBest.data[0][j] = particles[k].position.data[0][j];
                }
                particles[k].pBestFitness = particles[k].fitness;
            }
            /* 由于gbestval类型改为 引用
            if (particles[k].pBestFitness < *gbestval) {
                for (int j = 0; j < D; j++) {
                    (*gbest).data[0][j] = particles[k].pBest[j];
                }
                *gbestval = particles[k].pBestFitness;
            }
            */
            if (particles[k].pBestFitness < gbestval){
                for (int j=0; j<D; j++){
                    gbest.data[0][j] = particles[k].pBest.data[0][j];
                }
                gbestval = particles[k].pBestFitness;
            }
            freeMatrix(temp1);
            freeMatrix(temp2);
            freeMatrix(r1);
            freeMatrix(r2);
        }
    }
    
    free(particles);
}

