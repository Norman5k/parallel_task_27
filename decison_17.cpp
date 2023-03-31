#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "iostream"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
struct timeval tv1,tv2,dtv;
struct timezone tz;
void time_start() { gettimeofday(&tv1, &tz); }
long time_stop()
{ gettimeofday(&tv2, &tz);
  dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
  dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
  if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
  return dtv.tv_sec*1000+dtv.tv_usec/1000;
}
#ifndef __CUDACC__  
#define __CUDACC__
#endif

__global__ void kernel(double* dev_A,double* dev_B, double* dev_X, double* dev_P, int* dev_cnt, double eps, int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i >= n) 
		return;
	dev_P[i] = dev_X[i];
	double var = 0;
	for (int j = 0; j < n; j++)
		if (j != i) var += (dev_A[i*n+j] * dev_X[j]);
	dev_X[i] = (dev_B[i] - var) / dev_A[i*n + i];
	if(fabs(dev_X[i] - dev_P[i]) <= eps)
		dev_cnt[i] = 1;
}

double okr(double x, double eps) // Функция округления числа до заданной точности
{
	int i = 0;
	double neweps = eps;
	while (neweps < 1)
	{
		i++;
		neweps *= 10;
	}
	int okr = pow((double) 10, i);
	x = (int) (x * okr + 0.5) / (double) okr;

	return x;
}

bool diagonal(double* a, int n) // Проверка на диагональное преобладание
{
	int i, j;
	double sum;
	for (i = 0; i < n; i++)
	{
		sum = 0;
		for (j = 0; j < n; j++) sum += abs(a[i * n + j]);
		sum -= abs(a[i * n + i]);
		if (sum > a[i * n + i])
			return false;
	}
	return true;
}

void GenerateRand(double* a, double* b, int n) // Генерация матриц A и B
{
	double p = 100.0;
	for (int i = 0; i < n * n; i++)
	{
		if (i < n)
		{
			b[i] = (p * p * rand() / (RAND_MAX + 1.0)) + p * n * n;	
		}
		a[i] = (i % (n + 1) == 0) ? (p * rand() / (RAND_MAX + 1.0)) + (p - 1) * n : (p * rand() / (RAND_MAX + 1.0));
	}
}

int main(int argc, char* argv[])
{
	FILE* file_result;
	FILE* file_start;
	char filename_start[] = "start_matrix.txt";
	char filename_result[] = "result.txt";
	while (true)
	{
		srand(time(NULL));
		int n = 0;
		double eps, * b, * x, * p, * a, t, t2;
		printf("Введите ширину квадратной матрицы: ");
		scanf("%d", &n);
		printf("Введите точность вычислений: ");
		scanf("%lf", &eps);
		int size = 0;
		if(n <= 100)
			size = n;
		else
			size = 100;
		b = (double*) malloc(sizeof(double) * n);
		x = (double*) malloc(sizeof(double) * n);
		p = (double*) malloc(sizeof(double) * n);
		a = (double*) malloc(sizeof(double) * n * n);
		int i, j, iter = 0;
		do{
			GenerateRand(a,b,n); // Генерируем значения матриц A и B, которые дадут однозначное решение
		}while(!diagonal(a, n));
		if ((file_start = fopen(filename_start, "w")) == NULL) // Записываем начальные значения матриц A,B в файл
		{
			printf("Could not open file");
			return 1;
		}
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{ 
				fprintf(file_start,"A[%d][%d]=%lf, ",i,j, a[i * n + j]);
			}
			fprintf(file_start,"|B[%d]=%lf\n",i,b[i]);
		}
		fclose(file_start);
		/* Начало последовательной части*/
		for (int i = 0; i < n; i++) // Начальное приближение
			x[i] = 1;
		int count = 0;
		time_start();
			do
			{
				count = 0;
				for (int i = 0; i < n; i++)
					p[i] = x[i];
				for (int i = 0; i < n; i++)
				{
					double var = 0;
					for (int j = 0; j < n; j++)
						if (j != i) var += (a[i*n + j] * x[j]);

					x[i] = (b[i] - var) / a[i*n+i];
					if (fabs(x[i] - p[i]) <= eps)
						count++;
				}
				iter++;
			} while (count != n);
		printf("Time for line program: %ld ms\n", time_stop());
		printf("Iterations for line program: %d\n", iter);
		/* Конец последовательной части*/	
		/* Начало параллельной части*/
		double* x_new = (double*)malloc(n*sizeof(double));
		int* cnt = (int*)malloc(n*sizeof(int));
		cudaError_t cudaStatus;
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		}
		// Выделение памяти на устройстве под глобальные переменные.
		double* dev_A = 0;
		cudaStatus = cudaMalloc((void**) &dev_A, n * n * sizeof(double));
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc dev_A failed!");
		}
		double* dev_B = 0;
		cudaStatus = cudaMalloc((void**) &dev_B, n * sizeof(double));
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc dev_B failed!");
		}
		double* dev_X = 0;
		cudaStatus = cudaMalloc((void**) &dev_X, n * sizeof(double));
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc dev_X failed!");
		}
		double* dev_P = 0;
		cudaStatus = cudaMalloc((void**) &dev_P, n * sizeof(double));
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc dev_P failed!");
		}
		int * dev_cnt = 0;
		cudaStatus = cudaMalloc((void**) &dev_cnt,n * sizeof(int));
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc dev_cnt failed!");
		}
		// Копирование данных с хоста на устройство.
		cudaStatus = cudaMemcpy(dev_A, a, n * n * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy HTD dev_A failed!");
		}
		cudaStatus = cudaMemcpy(dev_B, b, n * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy HTD dev_B failed!");
		}
		cudaMemset(dev_X, 1, n * sizeof(double));
		cudaMemset(dev_P, 0, n * sizeof(double));
		cudaEvent_t Start, Stop;
		float time;
		cudaEventCreate(&Start);
		cudaEventCreate(&Stop);
		cudaEventRecord(Start, 0);
		dim3 blockDim = dim3(512), gridDim = dim3((n - 1) / 512 + 1);
		int sum = 0;
		iter = 0;
		do
		{
		sum = 0;
		iter++;
		cudaMemset(dev_cnt, 0, sizeof(int));
		kernel <<< gridDim, blockDim >>> (dev_A, dev_B, dev_X,dev_P,dev_cnt, eps, n);
		// Копирование данных с устройства на хост для проверки окончания итерационного метода
		cudaStatus = cudaMemcpy(cnt, dev_cnt, n * sizeof(int), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy DTH dev_cnt failed!");
			break;
		}
		for(int i = 0; i < n; i++)
			sum+=cnt[i];
		}while(sum != n);
		cudaEventRecord(Stop, 0);
		cudaEventSynchronize(Stop);
		cudaEventElapsedTime(&time, Start, Stop);
		printf("Time for parallel program: %f ms\n", time);
		printf("Iterations for parallel program: %d\n", iter);
		// Копирование окончательного решения с устройства на хост.
		cudaStatus = cudaMemcpy(x_new, dev_X, n * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy DTH dev_X failed!");
		}
		/*Конец параллельной части*/
		/*Запись результатов в файл*/
		if ((file_result = fopen(filename_result, "w")) == NULL)
			{
				printf("Could not open file");
			}
			fprintf(file_result, "Line_Result:\t\tParallel_Result:\t\tCompare:\n");
			for (i = 0; i < size; i++) 
			{ 
				t = okr(x[i], eps);
				t2 = okr(x_new[i], eps);
				fprintf(file_result, "x[%d] = %lf\t\t", i + 1, t);
				fprintf(file_result, "x[%d] = %lf\t\t", i + 1, t2);
				fprintf(file_result, "x[%d] = %lf\n", i + 1, fabs(t-t2));
			}
			fclose(file_result);
		/*Очистка памяти устройства и хоста*/	
		cudaEventDestroy(Start);
		cudaEventDestroy(Stop);
		cudaFree(dev_A);
		cudaFree(dev_B);
		cudaFree(dev_X);
		cudaFree(dev_P);
		cudaFree(dev_cnt);
		free(p);
		free(x);
		free(b);
		free(a);
		free(x_new);
		free(cnt);
		printf("Результаты работы сохранены в файл result.txt\n");
		printf("Введите \"0\" для выхода, либо другое число для перезапуска программы: ");
		int flag = -1;
		scanf("%d", &flag);
		if (flag == 0)
			break;
	}
}
