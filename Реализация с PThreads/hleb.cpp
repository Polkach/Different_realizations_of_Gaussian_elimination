using namespace std;

#include "kolbasa.h"

typedef struct
{
	int n;
	double **a;
	double **x;
	int my_rank;
	int total_threads;
} ARGS;

long int thread_time = 0;
double nm = 0.0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

//Вызов построения обратной матрицы
void *Inversion(void *p_arg)
{
	ARGS *arg = (ARGS*)p_arg;
	long int t1;

	time_start();
	InvMatrix(arg->n, arg->a, arg->x, arg->my_rank, arg->total_threads);
	t1 = time_stop();

	pthread_mutex_lock(&mutex);
	thread_time += t1;
	pthread_mutex_unlock(&mutex);

	return NULL;
}

//Вызов нахождения нормы
void *Multiplication(void *p_arg)
{
	ARGS *arg = (ARGS*)p_arg;
	long int t1;

	time_start();
	nm += Norma(arg->n, arg->a, arg->x, arg->my_rank, arg->total_threads);
	t1 = time_stop();

	pthread_mutex_lock(&mutex);
	thread_time += t1;
	pthread_mutex_unlock(&mutex);

	return NULL;
}

int main(int argc, char * argv[])
{
	int i;
	int n;
	double **a;
	double **x;
	//ifstream input;
	int total_threads=1;
	pthread_t *threads;
	ARGS *args;

	pthread_t *threads2;

    int t=5,k,j;
    double p;
    unsigned int skolko=2,otkuda=4;

//Такое же считывание, как и в однопоточном методе. Стандартной функцией atoi все так же пользоваться нельзя.
//Только в этот раз дополнительно указываем количество потоков
    for(i=1;i<argc;i++)
    {
        if(argv[i][0]=='-' && argv[i][1]=='g' && argv[i][2]=='\0') otkuda = 0;
        if(argv[i][0]=='-' && argv[i][1]=='m' && argv[i][2]=='\0') otkuda = 1;
        if(argv[i][0]=='-' && argv[i][1]=='f' && argv[i][2]=='\0') otkuda = 2;
    }
    if(otkuda==4) {cout << "-g - generirovat matricu sluchaino (nujen razmer, po umolchaniu - 2)" <<endl
        <<"-m - generirovat matricu po formule (nujen razmer, po umolchaniu - 2)" << endl
        <<"-f - otktir iz faila (nujno nazvanie faila)" <<endl
        <<"p=(chislo) - kolichestvo potokov (po umolchaniu - 1)"<< endl; return -1;}
	//cout << "Otkuda " << otkuda<<endl;

    for(i=1;i<argc;i++)
    {
        if(argv[i][0]=='p' && argv[i][1]=='=')
        {
            t = 0;
            total_threads = 0;
            for(j=2;j<(int)strlen(argv[i]);j++)
            {
                total_threads = total_threads*10 + int(argv[i][j]) - 48;
                if(int(argv[i][j])<48 || int(argv[i][j])>57) t = 1;
            }
            if(t==0) break;
        }
        if(t==0) break;
    }
    if(t==1 || total_threads<=0) total_threads = 1;
	//cout << "threads " << total_threads<<endl;

//В этот раз мы не запоминаем матрицу на этом этапе, а запоминаем способ генерации и количество строк
//(либо же файл, откуда считывать). Эту информацию мы подадим на вход "открывающей" функции
    if(otkuda==0)
    {
        for(i=1;i<argc;i++)
        {
            t = 0;
            skolko = 0;
            for(j=0;j<(int)strlen(argv[i]);j++)
            {
                skolko = skolko*10 + int(argv[i][j]) - 48;
                if(int(argv[i][j])<48 || int(argv[i][j])>57) t = 1;
            }
            if(t==0) break;
        }
        if(t==1) skolko = 2;
	n = skolko;
    }

    if(otkuda==1)
    {
        for(i=1;i<argc;i++)
        {
            t = 0;
            skolko = 0;
            for(j=0;j<(int)strlen(argv[i]);j++)
            {
                skolko = skolko*10 + int(argv[i][j]) - 48;
                if(int(argv[i][j])<48 || int(argv[i][j])>57) t = 1;
            }
            if(t==0) break;
        }
        if(t==1) skolko = 2;
	n = skolko;
    }

    if(otkuda==2)
    {
        t = 0;
        for(i=1;i<argc;i++)
        {
            if(ifstream(argv[i]))
            {
                skolko = i;
                n = visa(argv[i]);
                t = 1;
            }
        }
        if(t==0){cout << "Ne ukazan fail" << endl; return -2;}
    }

	a = new double*[n];
    for(i=0;i<n;i++) a[i] = new double[n];

//Эта функция занесет в массив 'а' матрицу указанным способом.
    Input1(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);

    x = new double*[n];
    for(i=0;i<n;i++) x[i] = new double[n];
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            x[i][j] = (double)(i == j);

//Создаем массив потоков
    threads = new pthread_t[total_threads];
    args = new ARGS[total_threads];

//Если что-то не создалось, то предупреждаем об этом и завершаем работу
    if (!(a && x && threads && args))
	{
		cout << "Not enough memory!" << endl;

		if (a)
            {for(i=0;i<n;i++) delete[]a[i];
            delete[]a;}
        if (x)
            {for(i=0;i<n;i++) delete[]x[i];
            delete[]x;}
		if (threads)
            delete[]threads;
		if (args)
            delete[]args;

		return -3;
	}

//Промежуточный вывод матрицы
	cout <<endl <<"Matrica:" << endl<<endl;
	Output1(n, a);

//Инициируем состояния
	for (i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].a = a;
		args[i].x = x;
		args[i].my_rank = i;
		args[i].total_threads = total_threads;
	}

//Засекаем время и запускаем нахождение обратной матрицы
	time_start();
//Если какой-то поток не запустится, то предупреждаем об этом и завершаем программу
	for (i = 0; i < total_threads; i++)
		if (pthread_create(threads + i, 0, Inversion, args + i))
		{
			cout << "Can not create thread " << i  << '!' << endl;

			if (a)
                {for(i=0;i<n;i++) delete[]a[i];
                delete[]a;}
            if (x)
                {for(i=0;i<n;i++) delete[]x[i];
                delete[]x;}
            if (threads)
                delete[]threads;
            if (args)
                delete[]args;

			return -4;
		}

//Ждем окончание работы потоков
	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			cout << "Can not wait thread " << i  << '!' << endl;

            if (a)
                {for(i=0;i<n;i++) delete[]a[i];
                delete[]a;}
            if (x)
                {for(i=0;i<n;i++) delete[]x[i];
                delete[]x;}
            if (threads)
                delete[]threads;
            if (args)
                delete[]args;

			return -5;
		}

//Останавливаем таймер и выводим полученные результаты
	double Time = time_stop();
	cout <<endl<< "Obratnaja matrica:" << endl<<endl;
	Output1(n, x);

    cout << endl << "Vrem`a raboti: " << (double)Time/1000.0<<" s" <<endl
    << "Vrem`a na potokah: " << (double)thread_time/1000.0 <<" s" <<endl
    << "Vrem`a na odnom potoke: " << (double)thread_time/total_threads/1000.0 <<" s"<< endl;

//Снова заносим оригинальную матрицу в массив 'а'
	Input2(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);

//Теперь используем параллельные вычисления для более быстрого подсчета квадратичной нормы
//разности единичной матрицы и произведения оригинальной матрицы на полученную обратную.
//Порядок действий тот же
    thread_time = 0;
	for (i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].a = a;
		args[i].x = x;
		args[i].my_rank = i;
		args[i].total_threads = total_threads;
	}

    time_start();

    for (i = 0; i < total_threads; i++)
		if (pthread_create(threads + i, 0, Multiplication, args + i))
		{
			cout << "Can not create thread " << i  << '!' << endl;

			if (a)
                {for(i=0;i<n;i++) delete[]a[i];
                delete[]a;}
            if (x)
                {for(i=0;i<n;i++) delete[]x[i];
                delete[]x;}
            if (threads2)
                delete[]threads;
            if (args)
                delete[]args;

			return -6;
		}

	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			cout << "Can not wait thread " << i  << '!' << endl;

            if (a)
                {for(i=0;i<n;i++) delete[]a[i];
                delete[]a;}
            if (x)
                {for(i=0;i<n;i++) delete[]x[i];
                delete[]x;}
            if (threads2)
                delete[]threads;
            if (args)
                delete[]args;

			return -7;
		}

	Time = time_stop();
    cout << endl << "Norma: " << sqrt(nm) << endl;

    cout << endl << "Vrem`a podsheta normi: " << (double)Time/1000.0<<" s" <<endl
    << "Vrem`a na potokah: " << (double)thread_time/1000.0 <<" s" <<endl
    << "Vrem`a na odnom potoke: " << (double)thread_time/total_threads/1000.0 <<" s"<< endl;
    cout<<endl;

//Освобождаем за собой память
	for(i=0;i<n;i++) delete[]a[i];
	delete[]a;
	for(i=0;i<n;i++) delete[]x[i];
	delete[]x;
	delete[]threads;
	delete[]args;

	return 0;
}
