#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <CL/cl_platform.h>
#include <CL/cl.h>
#include <sys/time.h>

struct timeval tv1,tv2,dtv;
struct timezone tz;
void time_start() { gettimeofday(&tv1, &tz); }
double time_stop()
{ gettimeofday(&tv2, &tz);
  dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
  dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
  if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
  return dtv.tv_sec*1000+dtv.tv_usec/1000;
}
void StartMatrix(double*** P, int k, int n, int m, double T1, double T2) // Начальные приближения
{
	int i,j,l;
	for (i = 0; i < k; i++) // Заполнение нижней грани нулями
	{
		for (j = 0; j < n; j++)
		{
			P[0][i][j] = 0;
		}
	}

	for (i = 0; i < k; i++) // Заполнение двух сторон верхней грани
	{
		P[m - 1][i][0] = T1;
		P[m - 1][i][n - 1] = T2;
	}

	double kf = (T2 - T1) / (n - 1);
	for (i = 0; i < k; i++) // Заполнение верхней грани
	{
		for (j = 0; j < n; j++)
			P[m - 1][i][j] = T1 + j * kf;
	}
	
	
		for (j = 0; j < k; j++) // Заполнение оставшейся части параллелепипеда
			for (l = 0; l < n; l++)
			{
				for (i = 0; i < m; i++)
				{ 
					kf = (P[m - 1][j][l]) / (m - 1);
					P[i][j][l] = kf * i;
				}
			}
	
	
}

double Update(double*** P,int k,int n, int m, int i, int j, int l) // Функция для обновления значения в последовательной программе
{
	double tmp = P[i-1][j][l];
	int cnt = 1;
	if (i != m - 1)
	{ 
		tmp += P[i+1][j][l];
		cnt++;
	}
	if (l != n-1)
	{ 
		tmp += P[i][j][l + 1];
		cnt++;
	}
	if (l != 0)
	{ 
		tmp += P[i][j][l - 1];
		cnt++;
	}
	if (j != k-1)
	{ 
		tmp += P[i][j + 1][l];
		cnt++;
	}
	if (j != 0)
	{ 
		tmp += P[i][j - 1][l];
		cnt++;
	}
	return tmp/cnt;
}

void Func(double* P, double*** Pold,int k, int n, int m) // переход от 3-мерного к 1-мерному массиву
{
	int c = 0;
	int i,j,l;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			for (l = 0; l < n; l++)
				P[c++] = Pold[i][j][l];
		}
	}
}

int SumCnt(int* A, int n) // Подсчёт количества подходящих элементов и зануление массива
{
	int sum = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		sum += A[i];
		A[i] = 0;
	}
	return sum;
}

const char* source =
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"__kernel void UpdateP(__global double *P, __global int *Counters, int k, int n, int m, double Eps) \n" // Функция-ядро для обновления значения в параллельной программе
"{ \n"
"unsigned int tid = get_global_id(0);\n"
"unsigned int i = tid / (k*n); \n" // высота m
"unsigned int j = (tid - i * k * n) / n; \n" // длина k
"unsigned int l = tid - i * k * n - j * n; \n" // ширина n
"if ((i == m - 1 && (l == 0 || l == n - 1)) || i == 0) \n"
"return; \n"
"double tmp = P[tid - k * n];\n"
"int cnt = 1;\n"
"if (i != m - 1)\n"
"{\n"
"tmp += P[tid + k*n];\n"
"cnt++;\n"
"}\n"
"if (l != n-1)\n"
"{\n"
"tmp += P[tid + 1];\n"
"cnt++;\n"
"}\n"
"if (l != 0)\n"
"{\n"
"tmp += P[tid - 1];\n"
"cnt++;\n"
"}\n"
"if (j != k-1)\n"
"{\n"
"tmp += P[tid + n];\n"
"cnt++;\n"
"}\n"
"if (j != 0)\n"
"{\n"
"tmp += P[tid - n];\n"
"cnt++;\n"
"}\n"
"double old = P[tid];\n"
"P[tid] = tmp / cnt;\n"
"if (fabs(old - P[tid]) < Eps)\n"
"Counters[tid] = 1;\n"
"else\n"
"Counters[tid] = 0;}\n";
int main(int argc, char** argv)
{
	int flag = 1;
	do
	{ 
	FILE* file;
	FILE* file2;
	FILE* file3;
	char name[] = "Line_result.txt";
	char name2[] = "Parall_result.txt";
	char name3[] = "Compare_result.txt";
	double Eps = 0.0; // точность подсчёта
	int k = 0, n = 0, m = 0; // k - длина, n - ширина, m - высота
	double T1 = 0.0, T2 = 0.0;
	printf("Enter the dimensions of the box:\n");
	printf("k(length): ");
	scanf("%d", &k);
	printf("n(width): ");
	scanf("%d", &n);
	printf("m(height): ");
	scanf("%d", &m);
	printf("Enter starting temperatures:\n");
	printf("T1: ");
	scanf("%lf", &T1);
	printf("T2: ");
	scanf("%lf", &T2);
	printf("Enter precision:\n");
	printf("Eps: ");
	scanf("%lf", &Eps);
	const int SIZE = k * n * m;
	int cnt = 0;
	int iter = 0;
	double tmp = 0;
	int i,j,l;
	double*** P = (double***)malloc(m * sizeof(double**));
	for (i = 0; i < m; i++)
	{
		P[i] = (double**)malloc(k * sizeof(double*));
		for (j = 0; j < k; j++)
		{
			P[i][j] = (double*)malloc(n * sizeof(double));
		}
	}
	double* pp = (double*)malloc(k * n * m * sizeof(double));
	StartMatrix(P, k, n, m, T1, T2);
	Func(pp, P, k, n, m);
	int* AllCounters = (int*)malloc(sizeof(int) * k * n * m);
	for (i = 0; i < k * n * m; i++)
	{
		AllCounters[i] = 0;
	}
	time_start();
	do
	{
		cnt = 0;
		iter++;
		for (i = 1; i < m; i++)
		{
			for (j = 0; j < k; j++)
			{
				for (l = 0; l < n; l++)
				{
					if (i == m - 1 && (l == 0 || l == n - 1))
						continue;
					tmp = P[i][j][l];
					P[i][j][l] = Update(P, k, n, m, i, j, l);
					if (fabs(P[i][j][l] - tmp) < Eps)
						cnt++;
				}
			}
		}
	} while (cnt != (n * m * k - k * n - 2 * k));
	printf("Time for line program: %lf s\n", time_stop()/1000);
	printf("Iterations for line program: %d\n", iter);
	
	cl_int clerr;
	cl_uint qty_platforms = 0;
	cl_platform_id* platforms;
	cl_uint ui;
	cl_uint* qty_devices;
	cl_device_id** devices;
	cl_context context1;
	cl_command_queue queue1;
	cl_program program1;
	cl_kernel kernel1;
	size_t sizeWork[1] = { SIZE };
	cl_mem arg1;
	cl_mem arg2;

	
	clerr = clGetPlatformIDs(0, NULL, &qty_platforms);
	if (qty_platforms == 0)
	{
		fprintf(stderr, "OpenCL not found!\n");
		return 1;
	}
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n",clerr);
		return 1;
	}
	platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * qty_platforms);
	devices = (cl_device_id**)malloc(sizeof(cl_device_id*) * qty_platforms);
	qty_devices = (cl_uint*)malloc(sizeof(cl_uint) * qty_platforms);
	clerr = clGetPlatformIDs(qty_platforms, platforms, NULL);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n", clerr);
		return 1;
	}
	for (ui = 0; ui < qty_platforms; ui++){
		clerr = clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL, 0, NULL, &qty_devices[ui]);
		if (qty_devices[ui])
		{
			devices[ui] = (cl_device_id*)malloc(qty_devices[ui] * sizeof(cl_device_id));
			clerr = clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL, qty_devices[ui], devices[ui], NULL);
			if (clerr != CL_SUCCESS)
			{
				fprintf(stderr, "OpenCL error: %d\n", clerr);
				return 1;
			}
		}
	}
	clerr = CL_SUCCESS;
	context1 = clCreateContext(0, qty_devices[0], devices[0], NULL, NULL, &clerr);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n", clerr);
		return 1;
	}
#pragma warning(suppress : 4996) 
	queue1 = clCreateCommandQueue(context1, devices[0][0], 0, &clerr);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n", clerr);
		return 1;
	}
	program1 = clCreateProgramWithSource(context1, 1, &source, NULL, &clerr);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n", clerr);
		return 1;
	}
	clerr = clBuildProgram(program1, 0, NULL, NULL, NULL, NULL);
	if (clerr == CL_BUILD_PROGRAM_FAILURE)
	{
		size_t log_size;
		clGetProgramBuildInfo(program1, devices[0][0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		char* log = (char*)malloc(log_size);
		clGetProgramBuildInfo(program1, devices[0][0], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
		fprintf(stderr, "%s\n", log);
		return 1;
	}
	kernel1 = clCreateKernel(program1, "UpdateP", &clerr);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "OpenCL error: %d\n", clerr);
		return 1;
	}
	int PSize = SIZE * sizeof(cl_double);
	int CSize = SIZE * sizeof(cl_int);
	arg1 = clCreateBuffer(context1, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, PSize, pp, NULL);
	clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void*)&arg1);
	clSetKernelArg(kernel1, 2, sizeof(int), &k);
	clSetKernelArg(kernel1, 3, sizeof(int), &n);
	clSetKernelArg(kernel1, 4, sizeof(int), &m);
	clSetKernelArg(kernel1, 5, sizeof(double), &Eps);
	sizeWork[0] = n * m * k;
	time_start();
	iter = 0;
	do
	{
		arg2 = clCreateBuffer(context1, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, CSize, AllCounters, NULL);
		iter++;
		clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void*)&arg2);
		clerr = clEnqueueNDRangeKernel(queue1,
			kernel1,
			1,
			NULL,
			sizeWork,
			NULL, 0, NULL, NULL);
		if (clerr != CL_SUCCESS)
		{
			fprintf(stderr, "\nOpenCL error: %d", clerr);
			return 1;
		}
		clerr = clEnqueueReadBuffer(queue1,
			arg2,
			CL_TRUE,
			0,
			CSize,
			AllCounters,
			0, NULL, NULL);
		if (clerr != CL_SUCCESS)
		{
			fprintf(stderr, "\nOpenCL error: %d", clerr);
			return 1;
		}
		clReleaseMemObject(arg2);
		cnt = SumCnt(AllCounters, k * n * m);
	} while (cnt != SIZE - k*n - 2*k);
	
	printf("Time for parallel program: %lf s\n", time_stop()/1000);
	printf("Iter for parallel program: %d\n", iter);
	clerr = clEnqueueReadBuffer(queue1,
		arg1,
		CL_TRUE,
		0,
		PSize,
		pp,
		0, NULL, NULL);
	if (clerr != CL_SUCCESS)
	{
		fprintf(stderr, "\nOpenCL error: %d", clerr);
		return 1;
	}
	if(k*m*n <= 5000)
	{ 
	if ((file = fopen(name, "w")) == NULL)
	{
		printf("Could not open file");
	}
	else
	{
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < k; j++)
			{
				for (l = 0; l < n; l++)
					fprintf(file, "P[%d][%d][%d] = %lf ", i, j, l, P[i][j][l]);
				fprintf(file, "\n");
			}
			fprintf(file, "\n\n");
		}
	}
	fclose(file);
	
	if ((file2 = fopen(name2, "w")) == NULL)
	{
		printf("Could not open file");
	}
	else
	{
		int c = 0;
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < k; j++)
			{
				for (l = 0; l < n; l++)
					fprintf(file2, "P[%d][%d][%d] = %lf ", i, j, l, pp[c++]);
				fprintf(file2, "\n");
			}
			fprintf(file2, "\n\n");
		}
	}
	fclose(file2);

	if ((file3 = fopen(name3, "w")) == NULL)
	{
		printf("Could not open file");
	}
	else
	{
		int c = 0;
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < k; j++)
			{
				for (l = 0; l < n; l++)
					fprintf(file3, "P[%d][%d][%d] = %lf ", i, j, l, fabs(P[i][j][l] - pp[c++]));
				fprintf(file3, "\n");
			}
			fprintf(file3, "\n\n");
		}
	}
	fclose(file3);
	}
	else
	{
		if ((file = fopen(name, "w")) == NULL)
		{
			printf("Could not open file");
		}
		else
		{
			cnt = 0;
			for (i = 1; i < m && cnt < 5000; i++)
			{
				for (j = 0; j < k && cnt < 5000; j++)
				{
					for (l = 0; l < n && cnt < 5000; l++)
					{ 
						fprintf(file, "P[%d][%d][%d] = %lf ", i, j, l, P[i][j][l]);
						cnt++;
					}
					fprintf(file, "\n");
				}
				fprintf(file, "\n\n");
			}
		}
		fclose(file);

		if ((file2 = fopen(name2, "w")) == NULL)
		{
			printf("Could not open file");
		}
		else
		{
			int c = k*n;
			for (i = 1; i < m && c < k*n + 5000; i++)
			{
				for (j = 0; j < k && c < k * n + 5000; j++)
				{
					for (l = 0; l < n && c < k * n + 5000; l++)
						fprintf(file2, "P[%d][%d][%d] = %lf ", i, j, l, pp[c++]);
					fprintf(file2, "\n");
				}
				fprintf(file2, "\n\n");
			}
		}
		fclose(file2);

		if ((file3 = fopen(name3, "w")) == NULL)
		{
			printf("Could not open file");
		}
		else
		{
			int c = k * n;
			for (i = 1; i < m && c < k * n + 5000; i++)
			{
				for (j = 0; j < k && c < k * n + 5000; j++)
				{
					for (l = 0; l < n && c < k * n + 5000; l++)
						fprintf(file3, "P[%d][%d][%d] = %lf ", i, j, l, fabs(P[i][j][l] - pp[c++]));
					fprintf(file3, "\n");
				}
				fprintf(file3, "\n\n");
			}
		}
		fclose(file3);
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			free(P[i][j]);
		}
		free(P[i]);
	}
	free(P);
	clReleaseMemObject(arg1);
	clReleaseKernel(kernel1);
	clReleaseProgram(program1);
	clReleaseCommandQueue(queue1);
	clReleaseContext(context1);
	free(qty_devices);
	free(devices);
	free(platforms);
	free(pp);
	free(AllCounters);
	printf("Enter any number to continue ('0' - to complete the program)\n");
	scanf("%d", &flag);
	}while (flag);
	return 0;
	}
