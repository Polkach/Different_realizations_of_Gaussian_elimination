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
int me;
int num;

double ravno(double x)
{
    if(fabs(x)<1e-7){
        return 0.0;}
    else {return x;}
}

int visa(const char * a)
{
int n=2;
ifstream ff(a);
ff >> n;
ff.close();
return n;
}

void Input(int skolko, double **a, int otkuda, const char * name)
{
    int i,j,n;
    double p;

    if (otkuda == 1)
    {
        n = skolko;
        for (i = 0; i < n/num; i++)
            for (j = 0; j < n; j++)
                a[i][j] = min(n-j,n-i-me*n/num);//1.0/(1.0+i+me*n/num+j);
    }
    if(otkuda==2)
    {
        ifstream f(name);
        f >> n;
        for(i=0;i<me*n/num;i++)
            for(j=0;j<n;j++)
            f>>p;
        for (i = 0; i < n/num; i++)
            for (j = 0; j < n; j++)
                f >> a[i][j];
        f.close();
    }
}

void Output(int n, double **a)
{
    int i,j,k = 0;

    for(i=0;i<min(n,6);i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(me==k)
        {
            for(j=0;j<min(n,6);j++)
                cout<<setw(12)<< a[i-me*n/num][j];
            cout << endl;
        }
        if(i-k*n/num==n/num-1) k++;
    }
}

int Inv(int n, double **a,double **x,int k)
{
    double *ms;
    ms = new double[n];
    double *mf;
    mf = new double[n];
    int i;

    for(i=0;i<n;i++) ms[i] = 0.0;
    for(i=me*n/num;i<(me+1)*n/num;i++) ms[i] = ravno(a[i-me*n/num][k]);

    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0;i<num;i++)
    MPI_Reduce(ms, mf, n, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   // if(me==0) for(i=0;i<n;i++)
    //    cout << mf[i] << ' ';
    //if(me==0) cout<<endl;

    int M = k;
    for(i=k+1;i<n;i++)
        if(fabs(mf[i]) > fabs(mf[M]))
        M = i;

    if(fabs(mf[M])<1e-7)
    {
        delete[]ms;
        delete[]mf;
        return -1;
    }


    delete[]ms;
    delete[]mf;

    int t1,t2;
    t1 = k/(n/num);
    t2 = M/(n/num);

    double *s1;
    s1 = new double[4*n];
    double *s2;
    s2 = new double[4*n];

    for(i=0;i<4*n;i++)
    {
        s1[i]=0.0;
    }

    if(me==t1)
        for(i=n;i<2*n;i++)
        {
            s1[i] = ravno(a[k-me*n/num][i-n]);
            s1[i+2*n] = ravno(x[k-me*n/num][i-n]);
        }

    if(me==t2)
        for(i=0;i<n;i++)
        {
            s1[i] = ravno(a[M-me*n/num][i]);
            s1[i+2*n] = ravno(x[M-me*n/num][i]);
        }

    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0;i<num;i++)
    MPI_Reduce(s1, s2, 4*n, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //if(me==0) for(i=0;i<4*n;i++)
    //    cout << s2[i] << ' ';
    //if(me==0) cout<<endl;
    double c = s2[k];
    for(i=0;i<n;i++)
    {
        s2[i] = ravno(s2[i]/c);
        s2[i+2*n] = ravno(s2[i+2*n]/c);
    }
   // if(me==0) for(i=0;i<4*n;i++)
    //    cout << s2[i] << ' ';
    //if(me==0) cout<<endl;
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

    delete[]s1;
    delete[]s2;

    return 1;
}

double Norma(double **a,double **x,int n)
{
    double * s;
    double * f;
    s = new double[2*n*n];
    f = new double[2*n*n];
    int i,j;
    for(i=0;i<2*n*n;i++) s[i] = 0;
    for(i=0;i<n;i++)
        for(j=0;j<n/num;j++) s[n*(me*(n/num)+j)+i] = a[j][i];
    for(i=0;i<n;i++)
        for(j=0;j<n/num;j++) s[n*(me*(n/num)+j)+i+n*n] = x[j][i];

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

    for(i=1;i<argc;i++)
    {
        if(argv[i][0]=='-' && argv[i][1]=='m' && argv[i][2]=='\0') otkuda = 1;
        if(argv[i][0]=='-' && argv[i][1]=='f' && argv[i][2]=='\0') otkuda = 2;
    }
    if(otkuda==4) {cout<<"-m - generirovat` matricu po formule (nujen razmer, po umolchaniu - 2)" << endl
        <<"-f - otkrit` iz faila (nujno nazvanie faila)" <<endl; return -1;}
    //cout << "Otkuda " << otkuda<<endl;

    cout << num << endl;

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

    a = new double*[n/num];
    for(i=0;i<n/num;i++) a[i] = new double[n];

    x = new double*[n/num];
    for(i=0;i<n/num;i++) x[i] = new double[n];
    for (i = 0; i < n/num; i++)
        for (j = 0; j < n; j++)
            if(i+me*n/num == j) x[i][j] = 1.0; else x[i][j] = 0.0;
            //x[i][j] = (double)(i+me*n/num == j);

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

    Input(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);

    if(me==0)
    {
    cout <<endl <<"Matrica:" << endl<<endl;
    }

    Output(n,a);
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0){startwtime = MPI_Wtime();}

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

    if (me==0)
    {
    endwtime = MPI_Wtime();

    cout <<endl<< "Obratnaja matrica:" << endl<<endl;}
    MPI_Barrier(MPI_COMM_WORLD);
    Output(n, x);
    MPI_Barrier(MPI_COMM_WORLD);
    if(me==0)
    {
    cout << endl << "Vrem`a raboti: " << (double)(endwtime-startwtime)<<" s" <<endl;

    //cout << endl << "Norma: " << nn << endl;
    }
    Input(skolko, a, otkuda, argv[min(skolko,100*(otkuda-1)*otkuda)]);
    double nn = Norma(a,x,n);
    if(me==0)
    {
    //cout << endl << "Vrem`a raboti: " << (double)(endwtime-startwtime)<<" s" <<endl;

    cout << endl << "Norma: " << nn << endl;
    }
    for(i=0;i<n/num;i++) delete[]a[i];
    delete[]a;
    for(i=0;i<n/num;i++) delete[]x[i];
    delete[]x;


    MPI_Finalize();

    return 0;
}

