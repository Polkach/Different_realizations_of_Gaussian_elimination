using namespace std;
#include "kolbasa.h"

//Считывает из файла размер матрицы
int visa(const char * a)
{
int n=2;
ifstream ff(a);
ff >> n;
ff.close();
return n;
}

//Заносит матрицу в массив 'а'
void Input1(int skolko, double **a, int otkuda, const char * name)
{
	int i,j,n;
	srand(time(NULL));

//Генерирует случайную матрицу
	if (otkuda == 0)
	{
	    n = skolko;
        ofstream g("gen.txt");//еще нам пригодится
		for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i][j] = rand()%1000/100.0;
                g << a[i][j] << ' ';
            }
            g << endl;
        }
	}

//Генерирует матрицу по формуле
	if (otkuda == 1)
	{
	    n = skolko;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				a[i][j] = 1.0/(1.0+i+j);
	}

//Считывает матрицу из файла
	if(otkuda==2)
    {
        ifstream f(name);
        f >> n;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				f >> a[i][j];
        f.close();
    }
}

//Выводит пересечение первых 5 строк и столбцов матрицы
void Output1(int n, double **a)
{
    int i,j;
    for(i=0;i<min(n,5);i++)
    {
        for(j=0;j<min(n,5);j++)
        cout << setw(15) << a[i][j];
        cout << endl;
    }
}

//Повторно считываем матрицу
int Input2(int skolko, double **a, int otkuda, const char * name)
{
	int i,j,n;

	if (otkuda == 0)
	{
	    n = skolko;
        ifstream g("gen.txt");
		for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                g >> a[i][j];
            }
        }
        remove("gen.txt");
	}
	if (otkuda == 1)
	{
	    n = skolko;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				a[i][j] = 1.0/(1.0+i+j);
	}
	if(otkuda==2)
    {
        ifstream f(name);
        f >> n;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				f >> a[i][j];
        f.close();
    }
}

//Синхронизация потоков
void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

//Построение обратной матрицы
int InvMatrix(int n, double **a, double **x, int my_rank, int total_threads)
{
	int i, j, k;
	int first_row;
	int last_row;
	double tmp;

	for (i = 0; i < n; i++)
	{
//Если мы в нулевом потоке, то находим главный элемент по строке и нормируем по нему строку
		if (my_rank == 0)
		{
			k = i;
			for (j = i + 1; j < n; j++)
				if (fabs(a[k][i]) < fabs(a[j][i]))
					k = j;

			for (j = 0; j < n; j++)
			{
				tmp = a[i][j];
				a[i][j] = a[k][j];
				a[k][j] = tmp;
			}

			for (j = 0; j < n; j++)
			{
				tmp = x[i][j];
				x[i][j] = x[k][j];
				x[k][j] = tmp;
			}

			tmp = a[i][i];

			tmp = 1.0/tmp;
			for (j = i; j < n; j++)
				a[i][j] *= tmp;
			for (j = 0; j < n; j++)
				x[i][j] *= tmp;
		}
		synchronize(total_threads);

//Распределяем равномерно строки по потокам. Каждый поток вычисляет, используя свой уникальный
//номер, какие строки ему обрабатывать
		first_row = (n - i - 1) * my_rank;
		first_row = first_row/total_threads + i + 1;
		last_row = (n - i - 1) * (my_rank + 1);
		last_row = last_row/total_threads + i + 1;

//Каждый поток отнимает выбранную строку от своих
		for (j = first_row; j < last_row; j++)
		{
			tmp = a[j][i];
			for (k = i; k < n; k++)
				a[j][k] -= tmp * a[i][k];
			for (k = 0; k < n; k++)
				x[j][k] -= tmp * x[i][k];
		}
		synchronize(total_threads);
	}

//Снова рапсределяем строки по потокам. Каждый поток вычисляет, используя свой уникальный
//номер, какие строки ему обрабатывать
	first_row = n * my_rank;
	first_row = first_row/total_threads;
	last_row = n * (my_rank + 1);
	last_row = last_row/total_threads;

//Каждый поток зануляет столбец у своих строк
	for (k = first_row; k < last_row; k++)
		for (i = n - 1; i >= 0; i--)
		{
			tmp = x[i][k];
			for (j = i + 1; j < n; j++)
				tmp -= a[i][j] * x[j][k];
			x[i][k] = tmp;
		}

	return 0;
}

//Параллельное вычисление нормы
double Norma(int n, double **a, double **x, int my_rank, int total_threads) //считаем норму матриц А*Х-Е
{
	int i, j, k;
	int first_row;
	int last_row;
	double tmp;
	double rezult;

    synchronize(total_threads);

//Равномерно распределяем строки между потоками. Каждый поток вычисляет, используя свой уникальный
//номер, какие строки ему обрабатывать
    first_row = n * my_rank;
	first_row = first_row/total_threads;
	last_row = n * (my_rank + 1);
	last_row = last_row/total_threads;

//Каждый поток прибавляет к норме квадраты элементов только своего куска
	rezult = 0.0;
	for (i = first_row; i < last_row; i++)
		for (j = 0; j < n; j++)
		{
			tmp = 0.0;
			for (k = 0; k < n; k++)
				tmp += a[i][k] * x[k][j];
			if (i == j)
				tmp -= 1.0;
			rezult += tmp * tmp;
		}

	return rezult;
}

struct timeval tv1,tv2,dtv;

struct timezone tz;

//Засекаем начало
void time_start() { gettimeofday(&tv1, &tz); }

//Засекаем конец
long time_stop()
{ gettimeofday(&tv2, &tz);
dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
return dtv.tv_sec*1000+dtv.tv_usec/1000;
}
