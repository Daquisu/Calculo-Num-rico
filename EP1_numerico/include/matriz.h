#ifndef MATRIZ_H
#define MATRIZ_H
#include <vector>
#include <iostream>
using namespace std;

class Matriz
{
    public:
        Matriz();
        Matriz( vector<vector<double> > lista_de_listas);
        Matriz(int n_linhas, int n_colunas);
        void show(int digitos_total, int digitos_depois_virgula);
        void show();
        virtual ~Matriz();
        int get_n_colunas();
        int get_n_linhas();
        Matriz transposta();
        Matriz pegarcoluna(int colunas);
        Matriz operator *(Matriz &outra);
        Matriz operator =(Matriz const&outra);
        Matriz operator -(Matriz const&outra);
        void mudar_elemento(double a_ij, int i, int j);
        void mudar_a(vector<vector<double> > a);
        double get_elemento(int i, int j);
        vector<vector<double> > get_a();
    protected:

    private:
        vector<vector<double> > a;
        int n_colunas;
        int n_linhas;
};

#endif // MATRIZ_H
