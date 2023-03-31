#include <stdio.h>
#include <mpi.h>
long p(long n,long k)
{
	if (k > n) return p(n, n);
	else
		if (k > 0) return p(n, k - 1) + p(n - k, k);
		else
			if (n == 0) return 1;
			else
				return 0;
}
int main(int argc, char ** argv)
{
	double start_time;
	long long N;
	int a = 1;
	MPI_Init(&argc, &argv);
	int num_procs, proc_id;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	int* answers = (int*) malloc(sizeof(int) * num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	int flag = 1;
	while (flag)
	{
		if (proc_id == 0)
		{
			printf("Введите N: ");
			fflush(stdout);
			scanf("%llu", &N);
			start_time = MPI_Wtime();
		}
		MPI_Bcast(&N, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		int i = proc_id + 1;
		while (p(i, i - 1) < N)
		{
			i += num_procs;
		}
		MPI_Gather(&i, 1, MPI_INT, answers, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (proc_id == 0)
		{
			a = answers[0];
			int i;
			for (i = 0; i < num_procs; i++)
			{
				if (answers[i] < a) 
					a = answers[i];
			}
			printf("Parallel time: %lf s\n", MPI_Wtime() - start_time);
			printf("Parallel answer: %d\n", a);
			printf("p(%d) = %ld\np(%d) = %ld\n", a, p(a, a - 1),a - 1, p(a - 1, a - 2));
			a = 1;
			start_time = MPI_Wtime();
			while (p(a, a - 1) < N)
			{
				a++;
			}
			printf("Line time: %lf s\n", (MPI_Wtime() - start_time));
			printf("Line answer: %d\n", a);
			printf("p(%d) = %ld\np(%d) = %ld\n", a, p(a, a - 1),a-1,p(a-1,a-2));
			printf("Введите \"0\" для выхода, либо другое число для перезапуска программы: ");
			fflush(stdout);
			scanf("%d", &flag);
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
