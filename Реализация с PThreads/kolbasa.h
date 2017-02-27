#ifndef __KOLBASA_H_INCLUDED__
#define __KOLBASA_H_INCLUDED__
#include <pthread.h>//библиотека для работы с потоками
#include <sys/resource.h>
#include <sys/time.h>//засекаем время
#include <iostream> //ввод и вывод
#include <stdexcept>//поиск логических, динамических и временных ошибок
#include <fstream>//чтение и запись в файл
#include <stdlib.h>//выделение памяти, перевод классов друг в друга
#include <math.h>//математические функции
#include <string.h>//размерные и безразмерные строки
#include <cstdio>//ввод и вывод по сишному
#include <iomanip>//setw
#include <time.h>//работа со временем
using namespace std;

int visa(const char * a);

void Input1(int skolko, double **a, int otkuda, const char * name);

void Output1(int n, double **a);

int Input2(int skolko, double **a, int otkuda, const char * name);

int InvMatrix(int n, double **a, double **x, int my_rank, int total_threads);

double Norma(int n, double **a, double **x, int my_rank, int total_threads);

void time_start();

long time_stop();

#endif /* __KOLBASA_H_INCLUDED__ */
