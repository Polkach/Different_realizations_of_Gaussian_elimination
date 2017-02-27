//Нахождение обратной матрциы с помощью присоединенной матрицы
#include "BEZ_ETOGO_VOOBSHE_NIKAK.h"
using namespace std;

int main(int argc, char * argv[])
{
    int n,i,t,k,j;
    double p;
    unsigned int skolko=2,otkuda=4;

    for(i=1;i<argc;i++)
    {
        if(argv[i][0]=='-' && argv[i][1]=='g' && argv[i][2]=='\0') otkuda = 0;
        if(argv[i][0]=='-' && argv[i][1]=='m' && argv[i][2]=='\0') otkuda = 1;
        if(argv[i][0]=='-' && argv[i][1]=='f' && argv[i][2]=='\0') otkuda = 2;
    }
    if(otkuda==4) {cout << "-g - generirovat matricu sluchaino (nujen razmer, po umolchaniu - 2)" <<endl
        <<"-m - generirovat matricu po formule (nujen razmer, po umolchaniu - 2)" << endl
        <<"-f - otktir iz faila (nujno nazvanie faila)" <<endl; return -1;}

//В случае генерирования матрицы по формуле или случайно считываем размер матрицы. По условию задания
//пользоваться библиотечной функцией atoi было запрещено. В переменную skolko будет занесен размер
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
        generator(skolko);
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
        formula(skolko);
    }

//В случае считывания из файла, мы пытаемся открыть указанный файл
    if(otkuda==2)
    {
        t = 0;
        for(i=1;i<argc;i++)
        {
            if(ifstream(argv[i]))
            {
                visa(argv[i]);
                t = 1;
            }
        }
        if(t==0){cout << "Ne ukazan fail" << endl; return -1;}
    }

    ifstream f("gen.txt");
    f >> n;

    double **a;
    a = new double*[n];
    for(i=0;i<n;i++) a[i] = new double[n];
    double **b;
    b = new double*[n];
    for(i=0;i<n;i++) b[i] = new double[n];

//Считываем нашу матрицу и создаем присоединенную матрицу
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        f >> a[i][j];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        b[i][j] = 0;
    for(i=0;i<n;i++)
        b[i][i] = 1;
    f.close();

//Засечем время работы алгоритма
    clock_t t1 = clock();

//Попытаемся свести оригинальную матрциу к треугольному виду
    i = treug(a,b,n);
    if(i==-1)
    {
        cout << "Net obratnoi" << endl;
        return -1;
    }

//Теперь приводим оригинальную матрицу к диагонаьному виду
    diag(a,b,n);

    clock_t t2 = clock();

//Вычислим квадратичную норму разности единичной матрицы и матрицы, полученной
//умножением оригинальной матрицы на полученную нами
    double norma = 0;
    ifstream g("gen.txt");
    g >> n;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        g >> a[i][j];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            p = 0;
            for(k=0;k<n;k++)
                p = p + a[i][k]*b[k][j];
            if(i==j) norma = norma + (p-1)*(p-1);
            else norma = norma + p*p;
        }
    }
    g.close();

//Вывод результатов
    cout << "Norma: " << sqrt(norma) << endl;
    cout << "Vrem`a: " << (double)(t2-t1)/ (double)CLOCKS_PER_SEC  << " s" << endl;

    ofstream h("rezult.txt");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        h << setw(12) << b[i][j];
        h << endl;
    }
    for(i=0;i<min(n,6);i++)
    {
        for(j=0;j<min(n,6);j++)
        cout << setw(12) << b[i][j];
        cout << endl;
    }
    h << endl <<  "Norma: " << sqrt(norma);

    for(i=0;i<n;i++) delete[]a[i];
    delete[]a;
    for(i=0;i<n;i++) delete[]b[i];
    delete[]b;
    h.close();
    return 0;
}
