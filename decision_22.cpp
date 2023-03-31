#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
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

struct Point
{
    double m;
    double r;
    double x;
    double y;
    double z;
    double ux;
    double uy;
    double uz;
};

void GenerateRand(struct Point* p, int N, FILE* fp)
{
    int a;
    int i;
    for (i = 0; i < N; i++)
    {
        p[i].m = (double)(rand() * rand() % 100000000) / 10 + 10000000.0; // массы назначаются в диапазоне от 10 000 000 до 19 999 999.9
        p[i].r = (double)(rand() % 100) / 10 + 1.0; // радиусы назначаются в диапазоне от 1 до 9.9
        p[i].x = (double)(rand() % N * 10) / 10 - N / 2; // координаты от -N/2 до N/2 
        p[i].y = (double)(rand() % N * 10) / 10 - N / 2;
        p[i].z = (double)(rand() % N * 10) / 10 - N / 2;
        p[i].ux = (double)(rand() % 1000) / 10 - 50.0; // скорости от -50 до 49.9 
        p[i].uy = (double)(rand() % 1000) / 10 - 50.0;
        p[i].uz = (double)(rand() % 1000) / 10 - 50.0;
    }
    int size = N;
    if (N > 1000)
        size = 1000;
    for (i = 0; i < size; i++)
    {
        fprintf(fp, "Object #%d:\n", i + 1);
        fprintf(fp, "m = %lf,r = %lf,x = %lf,y = %lf,z = %lf,ux = %lf,uy = %lf,uz = %lf\n", p[i].m, p[i].r, p[i].x, p[i].y, p[i].z, p[i].ux, p[i].uy, p[i].uz);
    }
    fclose(fp);
}

struct Point* NewPoint(struct Point* p, int N) // Увеличение размера массива точек
{
    struct Point* NewPoints = (struct Point*)malloc(sizeof(struct Point) * (N + 1));
    int i;
    for (i = 0; i < N; i++)
    {
        NewPoints[i] = p[i];
    }
    free(p);
    return NewPoints;
}

int* NewModify(int* arr, int size) // Увеличение размера массива модификаций
{
    int* NewModify = (int*)malloc(sizeof(int) * 2 * size);
    int i;
    for (i = 0; i < size; i++)
    {
        NewModify[i] = arr[i];
    }
    free(arr);
    return NewModify;
}

void SortModify(int* arr, int* arr_sec, int first, int last) // быстрая сортировка первого массива индексов модификаций
{
    int mid, tmp;
    int f = first, l = last;
    mid = arr[(f + l) / 2];
    do
    {
        while (arr[f] < mid) f++;
        while (arr[l] > mid) l--;
        if (f <= l)
        {
            tmp = arr[f];
            arr[f] = arr[l];
            arr[l] = tmp;
            tmp = arr_sec[f];
            arr_sec[f] = arr_sec[l];
            arr_sec[l] = tmp;
            f++;
            l--;
        }
    } while (f < l);
    if (first < l) SortModify(arr, arr_sec, first, l);
    if (f < last) SortModify(arr, arr_sec, f, last);
}

void SortModifySec(int* arr, int* arr_sec, int size) // сортировка второго массива модификаций
{
    int flag = 0;
    do
    {
        flag = 0;
        int i;
        for (i = 1; i < size; i++)
        {
            if (arr[i - 1] == arr[i] && arr_sec[i - 1] > arr_sec[i])
            {
                int tmp = arr[i - 1];
                arr[i - 1] = arr[i];
                arr[i] = tmp;
                tmp = arr_sec[i - 1];
                arr_sec[i - 1] = arr_sec[i];
                arr_sec[i] = tmp;
                flag++;
            }
        }
    } while (flag);

}

int main()
{
    const double G = 6.67 * pow(10, -11);
    while (1)
    {
        unsigned int N; // Количество тел
        unsigned int T; // Время симуляции
        unsigned int P; // Вероятность разделения тел
        int count = 0;
        int i, j, k;
        double temp = 0.0;
        double Fx = 0.0, Fy = 0.0, Fz = 0.0;
        srand(time(NULL));
        FILE* file_points;
        FILE* file_result;
        FILE* file_parall_result;
        FILE* file_compare;
        char filename_parall_result[] = "result_parall_points.txt";
        char filename_result[] = "result_line_points.txt";
        char filename_compare[] = "result_compare.txt";
        char filename_points[] = "start_points.txt";
        printf("Введите количество тел: ");
        scanf("%u", &N);
        unsigned int* N_p = (unsigned int*)malloc(sizeof(unsigned int));
        *N_p = N;
        struct Point* points = (struct Point*)malloc(sizeof(struct Point) * N);
        struct Point* parall_points = (struct Point*)malloc(sizeof(struct Point) * N);
        printf("Введите время симуляции: ");
        scanf("%u", &T);
        printf("Введите вероятность разделения тел: ");
        scanf("%u", &P);
        if ((file_points = fopen(filename_points, "w")) == NULL)
        {
            printf("Could not open file");
            return 1;
        }
        GenerateRand(points, N, file_points); // Генерируем начальное состояние тел
        for (i = 0; i < N; i++)
        {
            parall_points[i] = points[i];
        }
        int* modify = (int*)malloc(N * sizeof(int));
        int* modify_sec = (int*)malloc(N * sizeof(int));
        int size_mod = N;
        /* Начало последовательной части*/
        time_start();
        for (i = 0; i < T; i++)
        {
            for (j = 0; j < N; j++) // перемещение в момент времени T
            {
                points[j].x += points[j].ux;
                points[j].y += points[j].uy;
                points[j].z += points[j].uz;
            }
            for (j = 0; j < N; j++) // проверка на столкновение 
            {
                for (k = j + 1; k < N; k++)
                {
                    if ((points[j].x - points[j].r > points[k].x - points[k].r && points[j].x - points[j].r > points[k].x + points[k].r) ||
                        (points[j].x + points[j].r < points[k].x - points[k].r && points[j].x + points[j].r < points[k].x + points[k].r))
                        continue;
                    if ((points[j].y - points[j].r > points[k].y - points[k].r && points[j].y - points[j].r > points[k].y + points[k].r) ||
                        (points[j].y + points[j].r < points[k].y - points[k].r && points[j].y + points[j].r < points[k].y + points[k].r))
                        continue;
                    if ((points[j].z - points[j].r > points[k].z - points[k].r && points[j].z - points[j].r > points[k].z + points[k].r) ||
                        (points[j].z + points[j].r < points[k].z - points[k].r && points[j].z + points[j].r < points[k].z + points[k].r))
                        continue;
                    if (count == size_mod)
                    {
                        modify = NewModify(modify, size_mod);
                        modify_sec = NewModify(modify_sec, size_mod);
                        size_mod *= 2;
                    }
                    modify[count] = j; // Добавляем в массивы индексов столкнувшихся тел
                    modify_sec[count] = k;
                    count++;
                }
            }
            for (j = 0; j < count; j++) // Расчёт скоростей столкнувшихся тел
            {
                temp = points[modify[j]].ux;
                points[modify[j]].ux = ((points[modify[j]].m - points[modify_sec[j]].m) * temp + 2 * points[modify_sec[j]].m * points[modify_sec[j]].ux) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости первого тела по x
                points[modify_sec[j]].ux = ((points[modify_sec[j]].m - points[modify[j]].m) * points[modify_sec[j]].ux + 2 * points[modify[j]].m * temp) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости второго тела по x
                temp = points[modify[j]].uy;
                points[modify[j]].uy = ((points[modify[j]].m - points[modify_sec[j]].m) * temp + 2 * points[modify_sec[j]].m * points[modify_sec[j]].uy) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости первого тела по y
                points[modify_sec[j]].uy = ((points[modify_sec[j]].m - points[modify[j]].m) * points[modify_sec[j]].uy + 2 * points[modify[j]].m * temp) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости второго тела по y
                temp = points[modify[j]].uz;
                points[modify[j]].uz = ((points[modify[j]].m - points[modify_sec[j]].m) * temp + 2 * points[modify_sec[j]].m * points[modify_sec[j]].uz) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости первого тела по x
                points[modify_sec[j]].uz = ((points[modify_sec[j]].m - points[modify[j]].m) * points[modify_sec[j]].uz + 2 * points[modify[j]].m * temp) / (points[modify[j]].m + points[modify_sec[j]].m); // расчёт скорости второго тела по x
                if (rand() % 100 + 1 <= P) // С заданной вероятностью первое тело разбивается на два
                {
                    points = NewPoint(points, N); // Перезапись массива с увеличением его размера
                    points[modify[j]].m = points[modify[j]].m / 2;
                    points[modify[j]].r = points[modify[j]].r / 2;
                    points[N].m = points[modify[j]].m;
                    points[N].r = points[modify[j]].r;
                    points[N].x = (points[modify[j]].x + points[modify[j]].r > points[modify_sec[j]].x + points[modify_sec[j]].r ? points[modify[j]].x + points[modify[j]].r + points[N].r + 0.5 : points[modify_sec[j]].x + points[modify_sec[j]].r + points[N].r + 0.5); // Координаты устанавливаем правее, глубже и выше от исходных точек
                    points[N].y = (points[modify[j]].y + points[modify[j]].r > points[modify_sec[j]].y + points[modify_sec[j]].r ? points[modify[j]].y + points[modify[j]].r + points[N].r + 0.5 : points[modify_sec[j]].y + points[modify_sec[j]].r + points[N].r + 0.5); // Чтобы избежать повторного столкновения в этот же момент времени
                    points[N].z = (points[modify[j]].z + points[modify[j]].r > points[modify_sec[j]].z + points[modify_sec[j]].r ? points[modify[j]].z + points[modify[j]].r + points[N].r + 0.5 : points[modify_sec[j]].z + points[modify_sec[j]].r + points[N].r + 0.5);
                    points[N].ux = points[modify[j]].ux; // Направление и модуль векторов по осям x и y задаём такие же, как у первой половинки  
                    points[N].uy = points[modify[j]].uy;
                    points[N].uz = -points[modify[j]].uz; // Модуль по оси z тот же как у первой половинки, направление - противоположное
                    N++;
                }
            }
            if (count > 0) // Обнуляем массивы модификаций
            {
                size_mod = N;
                free(modify);
                free(modify_sec);
                modify = (int*)malloc(sizeof(int) * size_mod);
                modify_sec = (int*)malloc(sizeof(int) * size_mod);
            }
            count = 0;
            for (j = 0; j < N; j++) // расчёт влияния сил тяготения
            {
                for (k = 0; k < N; k++)
                {
                    if (k == j) continue;
                    if (round(points[j].x) == round(points[k].x) || (int)points[j].x == (int)points[k].x)
                        Fx += 0;
                    else if (points[k].x > points[j].x)
                        Fx += (G * ((points[k].m * points[j].m) / pow(points[k].x - points[j].x, 2)));
                    else
                        Fx -= (G * ((points[k].m * points[j].m) / pow(points[k].x - points[j].x, 2)));
                    if (round(points[j].y) == round(points[k].y) || (int)points[j].y == (int)points[k].y)
                        Fy += 0;
                    else if (points[k].y > points[j].y)
                        Fy += (G * ((points[k].m * points[j].m) / pow(points[k].y - points[j].y, 2)));
                    else
                        Fy -= (G * ((points[k].m * points[j].m) / pow(points[k].y - points[j].y, 2)));
                    if (round(points[j].z) == round(points[k].z) || (int)points[j].z == (int)points[k].z)
                        Fz += 0;
                    else if (points[k].z > points[j].z)
                        Fz += (G * ((points[k].m * points[j].m) / pow(points[k].z - points[j].z, 2)));
                    else
                        Fz -= (G * ((points[k].m * points[j].m) / pow(points[k].z - points[j].z, 2)));
                }
                points[j].ux += (Fx / points[j].m);
                points[j].uy += (Fy / points[j].m);
                points[j].uz += (Fz / points[j].m);
                Fx = 0.0;
                Fy = 0.0;
                Fz = 0.0;
            }
        }
        printf("Time Line: %ld ms\n", time_stop());
        if ((file_result = fopen(filename_result, "w")) == NULL)
        {
            printf("Could not open file");
            return 1;
        }
        if (N > 500) // Если много тел, то пишем последнюю тысячу, так как в конце массива находятся новые сгенерированные тела
        {
            int cnt = 0;
            for (i = N - 1; cnt != 500; i--, cnt++) // Начинаем запись с конца
            {
                fprintf(file_result, "Object #%d:\n", i + 1);
                fprintf(file_result, "m = %lf,r = %lf,x = %lf,y = %lf,z = %lf,ux = %lf,uy = %lf,uz = %lf\n", points[i].m, points[i].r, points[i].x, points[i].y, points[i].z, points[i].ux, points[i].uy, points[i].uz);
            }
        }
        else
        {
            for (i = 0; i < N; i++) // Начинаем запись с начала
            {
                fprintf(file_result, "Object #%d:\n", i + 1);
                fprintf(file_result, "m = %lf,r = %lf,x = %lf,y = %lf,z = %lf,ux = %lf,uy = %lf,uz = %lf\n", points[i].m, points[i].r, points[i].x, points[i].y, points[i].z, points[i].ux, points[i].uy, points[i].uz);
            }
        }
        fclose(file_result);
        /* Конец последовательной части*/
        /* Начало параллельной части*/
        int threadId, numthreads;
        count = 0;
        time_start();
        for (i = 0; i < T; i++)
        {
#pragma omp parallel default(none) private(threadId,Fx,Fy,Fz,temp,i,j,k) shared(numthreads,parall_points, N_p, count, modify, modify_sec,P,size_mod) 
            {
                threadId = omp_get_thread_num();
                numthreads = omp_get_num_threads();
                Fx = 0.0;
                Fy = 0.0;
                Fz = 0.0;
                for (j = threadId; j < *N_p; j += numthreads) // перемещение в момент T
                {
                    parall_points[j].x += parall_points[j].ux;
                    parall_points[j].y += parall_points[j].uy;
                    parall_points[j].z += parall_points[j].uz;
                }
#pragma omp barrier // Для правильной проверки на столкновение необходимо рассчитать все x,y,z в данный момент T
                for (j = threadId; j < *N_p; j += numthreads) // проверка на столкновение 
                {
                    for (k = j + 1; k < *N_p; k++)
                    {
                        if ((parall_points[j].x - parall_points[j].r > parall_points[k].x - parall_points[k].r && parall_points[j].x - parall_points[j].r > parall_points[k].x + parall_points[k].r) ||
                            (parall_points[j].x + parall_points[j].r < parall_points[k].x - parall_points[k].r && parall_points[j].x + parall_points[j].r < parall_points[k].x + parall_points[k].r))
                            continue;
                        if ((parall_points[j].y - parall_points[j].r > parall_points[k].y - parall_points[k].r && parall_points[j].y - parall_points[j].r > parall_points[k].y + parall_points[k].r) ||
                            (parall_points[j].y + parall_points[j].r < parall_points[k].y - parall_points[k].r && parall_points[j].y + parall_points[j].r < parall_points[k].y + parall_points[k].r))
                            continue;
                        if ((parall_points[j].z - parall_points[j].r > parall_points[k].z - parall_points[k].r && parall_points[j].z - parall_points[j].r > parall_points[k].z + parall_points[k].r) ||
                            (parall_points[j].z + parall_points[j].r < parall_points[k].z - parall_points[k].r && parall_points[j].z + parall_points[j].r < parall_points[k].z + parall_points[k].r))
                            continue;
#pragma omp critical // Директива необходима для того, чтобы два потока одновременно не записали информацию в один участок памяти
                        {
                            if (count == size_mod)
                            {
                                modify = NewModify(modify, size_mod);
                                modify_sec = NewModify(modify_sec, size_mod);
                                size_mod *= 2;
                            }
                            modify[count] = j;
                            modify_sec[count] = k;
                            count++;
                        }
                    }
                }
#pragma omp barrier // Чтобы начать расчёт скоростей ударяющихся тел и, с некоторой вероятностью P, генерацию новых тел, необходимо дождаться проверки каждого тела
                if (threadId == 0)
                {
		    if(count >= 2)
		    {
                    	SortModify(modify, modify_sec, 0, count - 1); // Для того чтобы расчёт производился в том же порядке, что и в последовательной программе
                    	SortModifySec(modify, modify_sec, count); // необходима сортировка массивов с индексами столкнувшихся тел
		    }
                    for (j = 0; j < count; j++)
                    {
                        temp = parall_points[modify[j]].ux;
                        parall_points[modify[j]].ux = ((parall_points[modify[j]].m - parall_points[modify_sec[j]].m) * temp + 2 * parall_points[modify_sec[j]].m * parall_points[modify_sec[j]].ux) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости первого тела по x
                        parall_points[modify_sec[j]].ux = ((parall_points[modify_sec[j]].m - parall_points[modify[j]].m) * parall_points[modify_sec[j]].ux + 2 * parall_points[modify[j]].m * temp) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости второго тела по x
                        temp = parall_points[modify[j]].uy;
                        parall_points[modify[j]].uy = ((parall_points[modify[j]].m - parall_points[modify_sec[j]].m) * temp + 2 * parall_points[modify_sec[j]].m * parall_points[modify_sec[j]].uy) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости первого тела по y
                        parall_points[modify_sec[j]].uy = ((parall_points[modify_sec[j]].m - parall_points[modify[j]].m) * parall_points[modify_sec[j]].uy + 2 * parall_points[modify[j]].m * temp) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости второго тела по y
                        temp = parall_points[modify[j]].uz;
                        parall_points[modify[j]].uz = ((parall_points[modify[j]].m - parall_points[modify_sec[j]].m) * temp + 2 * parall_points[modify_sec[j]].m * parall_points[modify_sec[j]].uz) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости первого тела по x
                        parall_points[modify_sec[j]].uz = ((parall_points[modify_sec[j]].m - parall_points[modify[j]].m) * parall_points[modify_sec[j]].uz + 2 * parall_points[modify[j]].m * temp) / (parall_points[modify[j]].m + parall_points[modify_sec[j]].m); // расчёт скорости второго тела по x
                        if (rand() % 100 + 1 <= P)
                        {
                            parall_points = NewPoint(parall_points, *N_p); // Перезапись массива с увеличением его размера
                            parall_points[modify[j]].m = parall_points[modify[j]].m / 2;
                            parall_points[modify[j]].r = parall_points[modify[j]].r / 2;
                            parall_points[*N_p].m = parall_points[modify[j]].m;
                            parall_points[*N_p].r = parall_points[modify[j]].r;
                            parall_points[*N_p].x = (parall_points[modify[j]].x + parall_points[modify[j]].r > parall_points[modify_sec[j]].x + parall_points[modify_sec[j]].r ? parall_points[modify[j]].x + parall_points[modify[j]].r + parall_points[*N_p].r + 0.5 : parall_points[modify_sec[j]].x + parall_points[modify_sec[j]].r + parall_points[*N_p].r + 0.5); // Координаты устанавливаем правее, глубже и выше от исходных точек
                            parall_points[*N_p].y = (parall_points[modify[j]].y + parall_points[modify[j]].r > parall_points[modify_sec[j]].y + parall_points[modify_sec[j]].r ? parall_points[modify[j]].y + parall_points[modify[j]].r + parall_points[*N_p].r + 0.5 : parall_points[modify_sec[j]].y + parall_points[modify_sec[j]].r + parall_points[*N_p].r + 0.5); // Чтобы избежать повторного столкновения в этот же момент времени
                            parall_points[*N_p].z = (parall_points[modify[j]].z + parall_points[modify[j]].r > parall_points[modify_sec[j]].z + parall_points[modify_sec[j]].r ? parall_points[modify[j]].z + parall_points[modify[j]].r + parall_points[*N_p].r + 0.5 : parall_points[modify_sec[j]].z + parall_points[modify_sec[j]].r + parall_points[*N_p].r + 0.5);
                            parall_points[*N_p].ux = parall_points[modify[j]].ux; // Направление и модуль векторов по осям x и y задаём такие же, как у первой половинки 
                            parall_points[*N_p].uy = parall_points[modify[j]].uy;
                            parall_points[*N_p].uz = -parall_points[modify[j]].uz; // Модуль по оси z тот же как у первой половинки, направление - противоположное
                            (*N_p)++;
                        }
                    }
                    if (count > 0)
                    {
                        size_mod = (*N_p);
                        free(modify);
                        free(modify_sec);
                        modify = (int*)malloc(sizeof(int) * size_mod);
                        modify_sec = (int*)malloc(sizeof(int) * size_mod);
                    }
                    count = 0;
                }
#pragma omp barrier // Все потоки ждут мастер-поток, для того чтобы продолжить работу с актуальными скоростями тел и количеством тел N_p
                for (j = threadId; j < *N_p; j += numthreads) // расчёт влияния сил тяготения
                {
                    for (k = 0; k < *N_p; k++)
                    {
                        if (k == j) continue;
                        if (round(parall_points[j].x) == round(parall_points[k].x) || (int)parall_points[j].x == (int)parall_points[k].x)
                            Fx += 0;
                        else if (parall_points[k].x > parall_points[j].x)
                            Fx += (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].x - parall_points[j].x, 2)));
                        else
                            Fx -= (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].x - parall_points[j].x, 2)));
                        if (round(parall_points[j].y) == round(parall_points[k].y) || (int)parall_points[j].y == (int)parall_points[k].y)
                            Fy += 0;
                        else if (parall_points[k].y > parall_points[j].y)
                            Fy += (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].y - parall_points[j].y, 2)));
                        else
                            Fy -= (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].y - parall_points[j].y, 2)));
                        if (round(parall_points[j].z) == round(parall_points[k].z) || (int)parall_points[j].z == (int)parall_points[k].z)
                            Fz += 0;
                        else if (parall_points[k].z > parall_points[j].z)
                            Fz += (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].z - parall_points[j].z, 2)));
                        else
                            Fz -= (G * ((parall_points[k].m * parall_points[j].m) / pow(parall_points[k].z - parall_points[j].z, 2)));
                    }
                    parall_points[j].ux += (Fx / parall_points[j].m);
                    parall_points[j].uy += (Fy / parall_points[j].m);
                    parall_points[j].uz += (Fz / parall_points[j].m);
                    Fx = 0.0;
                    Fy = 0.0;
                    Fz = 0.0;
                }
#pragma omp barrier
            }
        }
        printf("Time Parallel: %ld ms\n", time_stop());
        if ((file_parall_result = fopen(filename_parall_result, "w")) == NULL)
        {
            printf("Could not open file");
            return 1;
        }
        if ((file_compare = fopen(filename_compare, "w")) == NULL)
        {
            printf("Could not open file");
            return 1;
        }
        if (N > 500) // При большом количестве тел записываем последние 1000 тел, начиная с конца
        {
            int cnt = 0;
            for (i = (*N_p) - 1; cnt != 500; i--, cnt++)
            {
                fprintf(file_parall_result, "Object #%d:\n", i + 1);
                fprintf(file_parall_result, "m = %lf,r = %lf,x = %lf,y = %lf,z = %lf,ux = %lf,uy = %lf,uz = %lf\n", parall_points[i].m, parall_points[i].r, parall_points[i].x,
                    parall_points[i].y, parall_points[i].z, parall_points[i].ux, parall_points[i].uy, parall_points[i].uz);
                if(P == 100 || P == 0) // Не имеет смысла сравнивать результаты, при случайной генерации новых тел
		{
			fprintf(file_compare, "Object #%d:\n", i + 1);
                	fprintf(file_compare, "delta m = %lf,delta r = %lf,delta x = %lf,delta y = %lf,delta z = %lf,delta ux = %lf,delta uy = %lf,delta uz = %lf\n",
                    fabs(parall_points[i].m - points[i].m), fabs(parall_points[i].r - points[i].r), fabs(parall_points[i].x - points[i].x), fabs(parall_points[i].y - points[i].y),
                    fabs(parall_points[i].z - points[i].z), fabs(parall_points[i].ux - points[i].ux), fabs(parall_points[i].uy - points[i].uy), fabs(parall_points[i].uz - points[i].uz));
		}
            }
        }
        else // Записываем все тела, стартуя с начала
        {
            for (i = 0; i < *N_p; i++)
            {
                fprintf(file_parall_result, "Object #%d:\n", i + 1);
                fprintf(file_parall_result, "m = %lf,r = %lf,x = %lf,y = %lf,z = %lf,ux = %lf,uy = %lf,uz = %lf\n", parall_points[i].m, parall_points[i].r, parall_points[i].x,
                    parall_points[i].y, parall_points[i].z, parall_points[i].ux, parall_points[i].uy, parall_points[i].uz);
                 if(P == 100 || P == 0) // Не имеет смысла сравнивать результаты, при случайной генерации новых тел
		 {
			fprintf(file_compare, "Object #%d:\n", i + 1);
                	fprintf(file_compare, "delta m = %lf,delta r = %lf,delta x = %lf,delta y = %lf,delta z = %lf,delta ux = %lf,delta uy = %lf,delta uz = %lf\n",
                    fabs(parall_points[i].m - points[i].m), fabs(parall_points[i].r - points[i].r), fabs(parall_points[i].x - points[i].x), fabs(parall_points[i].y - points[i].y),
                    fabs(parall_points[i].z - points[i].z), fabs(parall_points[i].ux - points[i].ux), fabs(parall_points[i].uy - points[i].uy), fabs(parall_points[i].uz - points[i].uz));
            	 }
	     }
        }
        fclose(file_compare);
        fclose(file_parall_result);
        /* Конец параллельной части*/
        free(points);
        free(parall_points);
        free(N_p);
        free(modify);
        free(modify_sec);
        printf("Введите \"0\" для выхода, либо другое число для перезапуска программы: ");
        int flag = -1;
        scanf("%d", &flag);
        if (flag == 0)
            break;

    }
    return 0;
}
