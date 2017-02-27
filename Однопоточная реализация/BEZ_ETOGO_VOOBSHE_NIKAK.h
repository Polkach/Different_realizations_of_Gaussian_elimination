#include <iostream> //ввод и вывод
//#include <stdexcept>//поиск логических, динамических и временных ошибок
#include <fstream>//чтение и запись в файл
#include <stdlib.h>//выделение памяти, перевод классов друг в друга и рандом
#include <math.h>//математические функции
#include <string.h>//размерные и безразмерные строки
//#include <cstdio>//ввод и вывод по сишному
#include <iomanip>//setw
#include <time.h>//работа со временем
using namespace std;

//Генератор случайной матрицы
void generator(unsigned int n)
{
    ofstream f("gen.txt");
    f << n << endl;
    unsigned int i,j;
    srand(time(NULL));
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++) f << setw(10) << rand()%1000/100.0;
        f << endl;
    }
    f.close();
}
//Генератор матрицы по формуле
void formula(unsigned int n)
{
    ofstream f("gen.txt");
    f << n << endl;
    unsigned int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++){ f << setw(10);
        f<<1.0/(i+j+1.0);}
        f << endl;
    }
    f.close();
}
//Как и ранее мы сохраняем матрицу в файл gen.txt, используя для этого промежуточный файл
void visa(const char * a)
{
    ifstream infile(a);
    ofstream outfile("RAMAMBA_HARARAM_BURUM_AAA.txt");
    outfile << infile.rdbuf();
    infile.close();
    outfile.close();

    ifstream infil("RAMAMBA_HARARAM_BURUM_AAA.txt");
    ofstream outfil("gen.txt");
    outfil << infil.rdbuf();
    infil.close();
    outfil.close();

    remove("RAMAMBA_HARARAM_BURUM_AAA.txt");
}
//Сведение матрицы к треугольному виду
int treug(double** a,double** b,int n)
{
    double *c;
    c = new double[n];

    int i,j,k,t,l,m;
    double p;

    //Проверяем на каждом шаге, что матрица сводится к треугольному виду.
    //Иначе ранк не максимален и нет обратной матрицы.
    for(i=0;i<n;i++)
    {
        t = -1;
        for(j=i;j<n;j++)
            if(a[j][i]>1e-18 || a[j][i]<-1e-18) t = 1;
        if(t==-1)
        {
            cout << "Net obratnoi" << endl;
            return -1;
        }

        p = a[i][i];
        t = i;

    //Метод Гаусса с выбором главного элемента по строке
        for(j=i;j<n;j++)
        {
            if(fabs(a[j][i])>fabs(p))
            {
                p = a[j][i];
                t = j;
            }
        }

        for(j=0;j<n;j++) c[j] = a[t][j];
        for(j=0;j<n;j++) a[t][j] = a[i][j];
        for(j=0;j<n;j++) a[i][j] = c[j];

        for(j=0;j<n;j++) c[j] = b[t][j];
        for(j=0;j<n;j++) b[t][j] = b[i][j];
        for(j=0;j<n;j++) b[i][j] = c[j];

        p = a[i][i];

        for(j=0;j<n;j++)
        {
            a[i][j] = a[i][j]/p;
            b[i][j] = b[i][j]/p;
        }

        for(j=i+1;j<n;j++)
        {
            p = a[j][i];
            for(k=0;k<n;k++)
            {
                a[j][k] = a[j][k] - a[i][k]*p;
                b[j][k] = b[j][k] - b[i][k]*p;
            }
        }

    }
     for(l=0;l<min(n,6);l++)
    {
        for(m=0;m<min(n,6);m++)
        cout << setw(12) << a[l][m];
        cout << endl;
    }
    cout << endl;
    delete[]c;
    return 1;
}
//Сведение треугольной матрицы к диагональному виду
void diag(double** a,double** b,int n)
{
    int i,j,k;
    double p;
    for(i=n-1;i>0;i--)
    {
        for(j=0;j<i;j++)
        {
            p = a[j][i];
            for(k=0;k<n;k++)
            {
                a[j][k] = a[j][k] - a[i][k]*p;
                b[j][k] = b[j][k] - b[i][k]*p;
            }
        }
    }
}
