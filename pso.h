/*** 
 * @Author: zicen xiong
 * @Date: 2023-12-27 
 */
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double** data;
    int rows;
    int cols;
} Matrix;

typedef struct {
    Matrix position;
    Matrix velocity;
    Matrix pBest;
    double fitness;
    double pBestFitness;
} Particle;

/*** 
 * @description: 创建Matrix
 * @param {int} rows
 * @param {int} cols
 * @return {*}
 */
Matrix createMatrix(int rows, int cols);

/*** 
 * @description: 释放内存
 * @param {Matrix} mat
 * @return {*}
 */
void freeMatrix(Matrix mat);

/*** 
 * @description: 打印mat数据信息
 * @param {Matrix} mat
 * @return {*}
 */
void printMatrix(Matrix mat);
void printMatrix(const char* str, Matrix mat);

/*** 
 * @description: 获取矩阵mat的第row_idx行
 * @param {Matrix} mat
 * @param {int} row_idx
 * @return {Matrix} 
 */
Matrix get_row(Matrix mat, int row_idx);


/*** 
 * @description: 矩阵加法
 * @param {Matrix} mat1
 * @param {Matrix} mat2
 * @return {*}
 */
Matrix matrixAdd(Matrix mat1, Matrix mat2);

/*** 
 * @description: 矩阵连接
 * @param {Matrix} mat1
 * @param {Matrix} mat2
 * @return {*}
 */
Matrix matrixSubtract(Matrix mat1, Matrix mat2);

/*** 
 * @description: 矩阵数乘
 * @param {Matrix} mat
 * @param {double} scalar
 * @return {*}
 */
Matrix matrixScalarMultiply(Matrix mat, double scalar);

/*** 
 * @description: 矩阵乘法(元素相乘)
 * @param {Matrix} mat1
 * @param {Matrix} mat2
 * @return {*}
 */
Matrix matrixElementWiseMultiply(Matrix mat1, Matrix mat2);

/*** 
 * @description: 矩阵copy
 * @param {Matrix} mat
 * @return {*}
 */
Matrix matrixCopy(Matrix mat);

/*** 
 * @description: 矩阵转置
 * @param {Matrix} mat
 * @return {*}
 */
Matrix matrixTranspose(Matrix mat);

/*** 
 * @description: 生成随机数组 func=randondouble
 * @param {int} rows
 * @param {int} cols
 * @return {*}
 */
Matrix matrixNullaryExpr(int rows, int cols, double (*func)());

/*** 
 * @description: 生成随机的double 
 * @return {*}
 */
double randomDouble();

/*** 
 * @description: linspace函数
 * @param {int} size
 * @param {double} start
 * @param {double} end
 * @return {*}
 */
Matrix vectorLinSpaced(int size, double start, double end);

/*** 
 * @description: PSO算法的目标函数
 * @param {Matrix} position
 * @return {*}
 */
double fitnessFunction(Matrix position);

/*** 
 * @description: PSO 算法
 * @param {int} Dimension  变量数量
 * @param {int} Max_Gen    迭代最大次数
 * @param {int} Particle_Number 种群个体数量？
 * @param {Matrix} Rmin    下限
 * @param {Matrix} Rmax    上限
 * @param {Matrix&} gbest    返回值, 一维数组(1,n)
 * @param {double&} gbestval 返回值, 最优值, 单double 变量, 建议使用引用
 * @return {*}
 */
void CLPSO(int Dimension, int Max_Gen, int Particle_Number, Matrix Rmin, Matrix Rmax, Matrix& gbest, double& gbestval);