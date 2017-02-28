#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <sys/time.h>
using namespace std;
int me;//Номер потока
int num;//Всего потоков

//Если х очень близок к нулю, то заменяем на ноль. Спасает от ошибок
double ravno(double x)
{
    if(fabs(x)<1e-7){
        return 0.0;}
    else {return x;}
}

//Считывает из файла размер матрицы
int visa(const char * a)
{
int n=2;
ifstream ff(a);
ff >> n;
ff.close();
return n;
}

//Генерация матрицы по формуле, либо считывание матрицы из файла
void Input(int skolko, double **a, int otkuda, const char * name)
{
    int i,j,n;
    double p;

//Строки будут равномерно распределены среди потоков. Причем в каждый поток попадают подряд
//идущие строки
    if (otkuda == 1)
    {
        n = skolko;
        for (i = 0; i < n/num; i++)
            for (j = 0; j < n; j++)
                //В формуле учитываем, что для потока строки начинаются с нулевой, хотя
                //в оригинальной матрице эта строка (i+me*n/num)-ая
                a[i][j] = min(n-j,n-i-me*n/num);//1.0/(1.0+i+me*n/num+j);
    }
    if(otkuda==2)
    {
        ifstream f(name);
        f >> n;
        //До своего блока поток просто игнорирует строки
        for(i=0;i<me*n/num;i++)
            for(j=0;j<n;j++)
            f>>p;
        //Свой блок поток сохраняет себе в память
        for (i = 0; i < n/num; i++)
            for (j = 0; j < n; j++)
                f >> a[i][j];
        f.close();
    }
}

//Многопоточный вывод матрицы на экран
void Output(int n, double **a)
{
    int i,j,k = 0;

//Выводим только начало матрицы
    for(i=0;i<min(n,6);i++)
    {
    //Синхронизируем все потоки перед выводом очередной строки
        MPI_Barrier(MPI_COMM_WORLD);
        //Так как как мы распределили строки равномерно между потоками, причем в каждом потоке подряд идущие
        //строки, то сначала выводит свои строки первый поток. Если у него меньше чем 6 строк, то выводит второй
        //поток. И так далее
        if(me==k)
        {
            for(j=0;j<min(n,6);j++)
                cout<<setw(12)<< a[i-me*n/num][j];
            cout << endl;
        }
        if(i-k*n/num==n/num-1) k++;
    }
}

//Сведение оригинальной матрицы к диагональному виду и построение обратной из присоединенной матрицы.
//Обмен информацией между потоками занимает O(n) памяти, где n - размер матрицы. функция зануляет
//столбец под номером k (кроме диагонального элемента) с выбором главного элемента по строке
int Inv(int n, double **a,double **x,int k)
{
//Массивы, с помощью который будет передваться информация между потоками
    double *ms;
    ms = new double[n];
    double *mf;
    mf = new double[n];
    int i;

//Каждый поток вставляет в массив ms значения своего столбца, который диагонализруется в данный момент. На
//остальный места каждый поток вставляет ноль
    for(i=0;i<n;i++) ms[i] = 0.0;
    for(i=me*n/num;i<(me+1)*n/num;i++) ms[i] = ravno(a[i-me*n/num][k]);

//Потоки суммируют полученные столбцы ms в столбец mf. Получится, что столбец mf совпадает с рассматриваемым
//столбцом
    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0;i<num;i++)
    MPI_Reduce(ms, mf, n, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   // if(me==0) for(i=0;i<n;i++)
    //    cout << mf[i] << ' ';
    //if(me==0) cout<<endl;

//Ищем главный элемент в столбце
    int M = k;
    for(i=k+1;i<n;i++)
        if(fabs(mf[i]) > fabs(mf[M]))
        M = i;

//Если главный элемент равен нулю, то у матрицы нет обратной
    if(fabs(mf[M])<1e-7)
    {
        delete[]ms;
        delete[]mf;
        return -1;
    }

//Освобождаем память
    delete[]ms;
    delete[]mf;

    int t1,t2;
    t1 = k/(n/num);//номер потока, где содержится строка с тем же номером, что столбец, который мы диагонализируем
    t2 = M/(n/num);//номер потока с главным элементом

//Новые массивы для обмена информацией между потоками
    double *s1;
    s1 = new double[4*n];
    double *s2;
    s2 = new double[4*n];

//Все потоки зануляют созданные массивы
    for(i=0;i<4*n;i++)
    {
        s1[i]=0.0;
    }

//Поток, содержащий элемент на диагонали, заносит в массивы соответсвующие строки
    if(me==t1)
        for(i=n;i<2*n;i++)
        {
            s1[i] = ravno(a[k-me*n/num][i-n]);
            s1[i+2*n] = ravno(x[k-me*n/num][i-n]);
        }

//Поток с главным элементом заносит в массивы соответсвующие строки
    if(me==t2)
        for(i=0;i<n;i++)
        {
            s1[i] = ravno(a[M-me*n/num][i]);
            s1[i+2*n] = ravno(x[M-me*n/num][i]);
        }

//Снова суммируем все массивы в массив s2
    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0;i<num;i++)
    MPI_Reduce(s1, s2, 4*n, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //if(me==0) for(i=0;i<4*n;i++)
    //    cout << s2[i] << ' ';
    //if(me==0) cout<<endl;

//Все потоки нормируют строку с главным элементом
    double c = s2[k];
    for(i=0;i<n;i++)
    {
        s2[i] = ravno(s2[i]/c);
        s2[i+2*n] = ravno(s2[i+2*n]/c);
    }
   // if(me==0) for(i=0;i<4*n;i++)
    //    cout << s2[i] << ' ';
    //if(me==0) cout<<endl;

//Потоки с главным элементом и элементом на диагонали перестанавливают строки
    if(me==t2)
        for(i=n;i<2*n;i++)
        {
            a[M-me*n/num][i-n] = s2[i];
            x[M-me*n/num][i-n] = s2[i+2*n];
        }
    if(me==t1)
        for(i=0;i<n;i++)
        {
            a[k-me*n/num][i] = s2[i];
            x[k-me*n/num][i] = s2[i+2*n];
        }

//Теперь все потоки отнимают строку на диагонали от своих строк, тем самым зануляя
//столбец, с которым мы работаем (кроме диагонального элемента)
    int j;
    for(i=0;i<n/num;i++)
    {
        if(i+me*n/num!=k)
        {
            c = a[i][k];
            for(j=0;j<n;j++)
            {
                a[i][j] = ravno(a[i][j]-s2[j]*c);
                x[i][j] = ravno(x[i][j]-s2[j+2*n]*c);
            }
        }
    }

    //Output(n,a);
    //MPI_Barrier(MPI_COMM_WORLD);
    //if(me==0) cout<<endl;
    //Output(n,x);
    //MPI_Barrier(MPI_COMM_WORLD);
    //if(me==0) cout<<endl;

//Освобождаем память
    delete[]s1;
    delete[]s2;

    return 1;
}

//Функция перемножает матрицы а и х, отнимает от произведения единичную матрицу
//и считает квадратичную норму полученной матрицы. Здесь вычисления производятся на одном потоке,
//так как по условию эту часть распараллеливать не нужно было
double Norma(double **a,double **x,int n)
{
    double * s;
    double * f;
//Кажый поток вносит в массив для обмена свои части матриц
    s = new double[2*n*n];
    f = new double[2*n*n];
    int i,j;
    for(i=0;i<2*n*n;i++) s[i] = 0;
    for(i=0;i<n;i++)
        for(j=0;j<n/num;j++) s[n*(me*(n/num)+j)+i] = a[j][i];
    for(i=0;i<n;i++)
        for(j=0;j<n/num;j++) s[n*(me*(n/num)+j)+i+n*n] = x[j][i];

//Все потоки получают полные матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0;i<num;i++)
    MPI_Reduce(s, f, 2*n*n, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    double ** m1;
    double ** m2;
    m1 = new double*[n];
    for(i=0;i<n;i++) m1[i] = new double[n];
    m2 = new double*[n];
    for(i=0;i<n;i++) m2[i] = new double[n];

    for(i=0;i<n*n;i++)
    {
        m1[i/n][i%n] = f[i];
        m2[i/n][i%n] = f[i+n*n];
    }

//Считаем норму
    double t,Norm = 0.0;
    int k;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            t=0.0;
            for(k=0;k<n;k++)
                t = t+m1[i][k]*m2[k][j];
            if(i==j)
                Norm = Norm+(t-1.0)*(t-1.0);
            else
                Norm = Norm+t*t;
        }

//Чистим память и выдаем значение нормы
    for(i=0;i<n;i++) delete[]m1[i];
    delete[]m1;
    for(i=0;i<n;i++) delete[]m2[i];
    delete[]m2;
    delete[]s;
    delete[]f;
    return sqrt(Norm);
}

int main(int argc, char * argv[])
{
//Каждый поток инициализирует свой номер и число потоков
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    int i;
    int n;
    double **a;
    double **x;
    ifstream input;
    double startwtime = 0.0;
    double endwtime;
    int t,j;
    unsigned int skolko=2,otkuda=4;
    bool need_norm=false;

//В этот раз обойдемся без случайной генерации матрицы
    for(i=1;i<argc;i++)
    {
        if(argv[i][0]=='-' && argv[i][1]=='m' && argv[i][2]=='\0') otkuda = 1;
        if(argv[i][0]=='-' && argv[i][1]=='f' && argv[i][2]=='\0') otkuda = 2;
        if(argv[i][0]=='-' && argv[i][1]=='n' && argv[i][2]=='\0') need_norm = true;
    }
    if(otkuda==4) {cout<<"-m - generirovat` matricu po formule (nujen razmer, po umolchaniu - 2)" << endl
        <<"-f - otkrit` iz faila (nujno nazvanie faila)" <<endl
        <<"-n - esli nujno schitat' normu" <<endl; return -1;}
    //cout << "Otkuda " << otkuda<<endl;

//Каждый поток напишет, сколько всего потоков. Так сразу визуально продиагностируем, что все потоки работают
    cout << num << endl;

//Если генерируем по формуле, то считываем размер генерируемой матрицы. Функцией atoi по условию было пользоваться
//запрещено
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

//Если считываем из файла, то пытаемся открыть его
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

//Итак, равномерно распределяем строки между потоками
    a = new double*[n/num];
    for(i=0;i<n/num;i++) a[i] = new double[n];
    x = new double*[n/num];
    for(i=0;i<n/num;i++) x[i] = new double[n];

//Каждый поток заполняет свою часть присоединенной матрицы
    for (i = 0; i < n/num; i++)
        for (j = 0; j < n; j++)
            if(i+me*n/num == j) x[i][j] = 1.0; else x[i][j] = 0.0;
            //x[i][j] = (double)(i+me*n/num == j);

//Проверяем, чтобы всем хватило памяти
    if (!(a && x))
    {
        cout << "Not enough memory!" << endl;

        if (a)
            {for(i=0;i<n/num;i++) delete[]a[i];
            delete[]a;}
        if (x)
            {for(i=0;i<n/num;i++) delete[]x[i];
            delete[]x;}

        return -3;
    }

//Заполняем оригинальную матрицу
    Input(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);

//Выводим матрицу на экран
    if(me==0)
    {
    cout <<endl <<"Matrica:" << endl<<endl;
    }
    Output(n,a);

//Синхронизируемся и засекаем время
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0){startwtime = MPI_Wtime();}

//Запускаем поиск обратной матрицы. Поочередно каждый столбец приводим к диагональному виду
    for(i=0;i<n;i++)
    {
        j = Inv(n,a,x,i);
        if(j==-1)
        {
            cout << "Net obratnoi" << endl;
            for(i=0;i<n/num;i++) delete[]a[i];
            delete[]a;
            for(i=0;i<n/num;i++) delete[]x[i];
            delete[]x;
            return -4;
        }
    }

//Останавливаем таймер
    if (me==0)
    {
    endwtime = MPI_Wtime();

//Выводим на экран обратную матрицу и время работы
    cout <<endl<< "Obratnaja matrica:" << endl<<endl;}
    MPI_Barrier(MPI_COMM_WORLD);
    Output(n, x);
    MPI_Barrier(MPI_COMM_WORLD);
    if(me==0)
    {
    cout << endl << "Vrem`a raboti: " << (double)(endwtime-startwtime)<<" s" <<endl;

    //cout << endl << "Norma: " << nn << endl;
    }

//Если требуется, вычисляем норму
    if(need_norm)
    {
    //Снова считываем оригинальную матрицу для подсчета нормы
        Input(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);

    //Вычисляем норму
        double nn = Norma(a,x,n);

    //Выводим на экран получившееся значение
        if(me==0)
        {
        //cout << endl << "Vrem`a raboti: " << (double)(endwtime-startwtime)<<" s" <<endl;

        cout << endl << "Norma: " << nn << endl;
        }
    }

//Чистим память
    for(i=0;i<n/num;i++) delete[]a[i];
    delete[]a;
    for(i=0;i<n/num;i++) delete[]x[i];
    delete[]x;

//Завершаем работу
    MPI_Finalize();

    return 0;
}

