#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector> 
#include <time.h>
using namespace std;

int l = 2;
clock_t starttime, endtime;
double cpu_time_used1;
double cpu_time_used2;
double cpu_time_used3;
//loat tempC[l][l];

double fronorm(vector<vector<float>> &arrA, vector<vector<float>> &arrB, int size1, int size2)
{
    double result = 0.0;
    for(int i = 0; i < size1; ++i)
    {
        for(int j = 0; j < size2; ++j)
        {
            double value = arrA[i][j] - arrB[i][j];
            result += value * value;
        }
    }
    return sqrt(result);
}

int matmul(int si,int ei,int sj,int ej,int sk,int ek, vector<vector<float>> &arrA, vector<vector<float>> &arrB, vector<vector<float>> &arrC)
{
    int rows = ei - si +1;
    int cols = ej - sj +1;
    float tempC[rows][cols];

    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            tempC[i][j]=0;

    for (int i = si; i <= ei; ++i)
        for (int j = sj; j <= ej; ++j)
            for (int k = sk; k <= ek; ++k)
                tempC[i][j] += arrA[i][k] * arrB[k][j];

    for(int i = si; i <= ei; ++i)
        for(int j = sj; j <= ej; ++j)
            arrC[i][j]+=tempC[i][j];

    return 0;
}

int matmul_NP(int A , int B, int C, vector<vector<float>> &arrA, vector<vector<float>> &arrB, vector<vector<float>> &arrC)
{
    for(int i = 0; i < A; ++i)
        for(int j = 0; j < C; ++j)
            arrC[i][j]=0;

    for (int i = 0; i < A; ++i)
        for (int j = 0; j < C; ++j)
            for (int k = 0; k < B; ++k)
                arrC[i][j] += arrA[i][k] * arrB[k][j];

    return 0;
}

/*
int matmul_recur(int lm, int hm, int ln, int hn, int lr, int hr, vector<vector<float>> &arrA, vector<vector<float>> &arrB, vector<vector<float>> &arrC)
{
    int midm = floor((hm - lm + 1)/2);
    int midn = floor((hn - ln + 1)/2);
    int midr = floor((hr - lr + 1)/2);

    c111 = matmul_recur(lm,midm,ln,midn,lr,midr,arrA,arrB);
    c112 = matmul_recur(lm,midm,1+midn,hn,lr,midr,arrA,arrB);
    c121 = matmul_recur(lm,midm,ln,midn,midr+1,hr,arrA,arrB);
    c122 = matmul_recur(lm,midm,1+midn,hn,1+midr,hr,arrA,arrB);
    c211 = matmul_recur(1+midm,hm,ln,midn,lr,midr,arrA,arrB);
    c212 = matmul_recur(1+midm,hm,1+midn,hn,lr,midr,arrA,arrB);
    c221 = matmul_recur(1+midm,hm,ln,midn,midr+1,hr,arrA,arrB);
    c222 = matmul_recur(1+midm,hm,1+midn,hn,1+midr,hr,arrA,arrB);

    c111 +c112;
    c121 + c122;
    c211 + c212;
    c221 + c222;

    return 0;
}
*/

int matmul_recur(int ln, int hn, int lm, int hm, int lp, int hp, vector<vector<float>> &arrA, vector<vector<float>> &arrB, vector<vector<float>> &arrC)
{
    int n, m, p;
    n = hn - ln + 1;
    m = hm - lm + 1;
    p = hp - lp + 1;
    //cout<<"staty"<<n<<" "<<m<<" "<<p<<endl;
    int max = ( n >= m ) ? n : m;
    max = ( max >= p ) ? max : p;
    
    if (max > l)
    {
        if (max == n)
        {
            cout<<"hern";
            int midn = ln + floor(float(hn - ln + 1)/2);
            cout<<midn<<endl;
            matmul_recur(ln,midn,lm,hm,lp,hp,arrA,arrB,arrC);
            matmul_recur(1+midn,hn,lm,hm,lp,hp,arrA,arrB,arrC);
        }
        else if (max == p)
        {
            cout<<"herp";
            int midp = lp + floor((float)(hp - lp + 1)/2);
            matmul_recur(ln,hn,lm,hm,lp,midp,arrA,arrB,arrC);
            matmul_recur(ln,hn,lm,hm,1+midp,hp,arrA,arrB,arrC);
        }
        else
        {
            int midm = lm + floor((float)(hm - lm + 1)/2);
            matmul_recur(ln,hn,lm,midm,lp,hp,arrA,arrB,arrC);
            matmul_recur(ln,hn,1+midm,hm,lp,hp,arrA,arrB,arrC);
        }
    }
    else
    {
        matmul(ln,hn,lp,hp,lm,hm,arrA,arrB,arrC);
    }
    return 0;
}

int matmul_1_level_part(int A , int B, int C, vector<vector<float>> &arrA, vector<vector<float>> &arrB, vector<vector<float>> &arrC)
{
    
    int x = ceil((float)A/l);
    int y = ceil((float)C/l);
    int z = ceil((float)B/l);
    //cout<<x<<" "<<y<<" "<<z<<endl;
    for(int a=1; a<=x ; a++)
    {
        int starti = (a-1)*l;
        int endi;
        if (a*l < A)
            endi = a*l-1;
        else
            endi = A-1;
         
        for(int b = 1;b <= y; b++)    
        {
            int startj = (b-1)*l;
            int endj;
            if (b*l < C)
                endj = b*l-1;
            else
                endj = C-1;

            for (int e = starti; e <= endi; ++e)
                for (int f = startj; f <= endj; ++f)
                    arrC[e][f] = 0;

            for(int c = 1;c <= z;c++)
            {
                int startk = (c-1)*l;
                int endk;
                if (c*l < B)
                    endk = c*l-1;
                else
                    endk = B-1;

                matmul(starti,endi,startj,endj,startk,endk,arrA,arrB,arrC);
            }
        }
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    int P,Q,R;
    vector<vector<float>> arrA;
    vector<vector<float>> arrB;
    vector<vector<float>> arrC;
    vector<vector<float>> arrC_RP;
    vector<vector<float>> arrC_NP;

    ifstream infile;
    infile.open ("matA.txt");
    if (infile.is_open())
    {
        infile >> P >> Q;
        arrA.resize(P);
        for(int i = 0 ; i < P ; ++i)
            arrA[i].resize(Q);

        for(int i = 0 ; i < P ; ++i)
            for(int j = 0 ; j < Q ; ++j)
                infile >> arrA[i][j];
    }
    infile.close();
    
    
    infile.open ("matB.txt");
    if (infile.is_open())
    {
        infile >> Q >> R;
        arrB.resize(Q);
        for(int i = 0 ; i < Q ; ++i)
            arrB[i].resize(R);

        for(int i = 0 ; i < Q ; ++i)
            for(int j = 0 ; j < R ; ++j)
                infile >> arrB[i][j];
    }
    infile.close();

    arrC.resize(P);
    for(int i = 0 ; i < R ; ++i)
        arrC[i].resize(R);

    arrC_NP.resize(P);
    for(int i = 0 ; i < R ; ++i)
        arrC_NP[i].resize(R);

    arrC_RP.resize(P);
    for(int i = 0 ; i < R ; ++i)
        arrC_RP[i].resize(R);

    for (int e = 0; e < P; ++e)
        for (int f = 0; f < R; ++f)
            arrC_RP[e][f] = 0;


    starttime = clock();
    matmul_NP(P,Q,R,arrA,arrB,arrC_NP);
    endtime = clock();
    cpu_time_used1 = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;

    starttime = clock();
    matmul_1_level_part(P,Q,R,arrA,arrB,arrC);
    endtime = clock();
    cpu_time_used2 = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;

    starttime = clock();
    matmul_recur(0,P-1,0,Q-1,0,R-1,arrA,arrB,arrC_RP);
    endtime = clock();
    cpu_time_used3 = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;

    /*matmul_1_level_part(P,Q,R,arrA,arrB,arrC);
    matmul_recur(0,P-1,0,Q-1,0,R-1,arrA,arrB,arrC_RP);
    matmul_NP(P,Q,R,arrA,arrB,arrC_NP);
/*
        for (int e = 0; e < P; ++e){
        for (int f = 0; f < R; ++f)
            cout << arrC_NP[e][f] << " ";
        cout<<"\n";
    }
*/
    ofstream outfile;
    outfile.open ("matC_P.txt");
    for(int i = 0 ; i < P ; ++i)
    {
        for(int j = 0 ; j < R ; ++j)
            outfile << arrC[i][j] << " ";
        outfile<<"\n";
    }
    outfile.close();
    
    outfile.open ("matC_N.txt");
    for(int i = 0 ; i < P ; ++i)
    {
        for(int j = 0 ; j < R ; ++j)
            outfile << arrC_NP[i][j] << " ";
        outfile<<"\n";
    }
    outfile.close();

    outfile.open ("matC_RP.txt");
    for(int i = 0 ; i < P ; ++i)
    {
        for(int j = 0 ; j < R ; ++j)
            outfile << arrC_RP[i][j] << " ";
        outfile<<"\n";
    }
    outfile.close();

    outfile.open ("output.txt");
    outfile << cpu_time_used1 << endl<< cpu_time_used2<<endl<< cpu_time_used3<<endl<<fronorm(arrC_NP,arrC_RP,P,R);
    outfile<<"\n";
    outfile.close();

    return 0;
}