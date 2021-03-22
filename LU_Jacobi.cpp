#include <iostream>
#include <vector>

using namespace std;




void mineurs(vector<double> a, vector<double> &delt){
int N = a.size();
int n = (N+2)/3;
int i = 2;
delt[0] = 1;
delt[1] = a[n-1];


for(int i = 2;i<=n;i++){
    delt[i] = a[i-2+n - 1]*delt[i-1]- a[i-2]*a[i-2+2*n-1]*delt[i-2];
}

}

void decLU(vector<double> a, vector<double> delt, float L[4][4], float U[4][4]){
//la fonction renvoie les matrice L, U pour une matrice tridiagonale
//represente par un vector a = (a_1,a_2,..,a_{n-1},b_1,..,b_n,c_1,..,c_{n-1})
int N = a.size();
int n = (N+2)/3;

vector<double> alpha(n+1);
int k =1;


//on calcule les alpha_k
// le vector alpha commence a 1
for(int k =1; k<=n;k++){
    alpha[k] = delt[k]/delt[k-1];
}

//on va contruire les deux de la decomposition LU


int i = 0;
int j = 0;
int idx_a = 0;

for(int i = 0;i<n;i++){
    for(int j = 0;j<n;j++){
        if(i==j){
            L[i][j] = 1;
        }else if(i-j == 1){
            L[i][j] = a[idx_a]/alpha[idx_a+1];
            idx_a++;
        }else{
        L[i][j] = 0;}

    }
}

int idx_c = 2*n-1;

for(int i = 0;i<n;i++){
    for(int j = 0;j<n;j++){
        if(i==j){
            U[i][j] = alpha[i+1];
            idx_a++;
        }else if(j-i == 1){
            U[i][j] = a[idx_c] ;
            idx_c++;
        }else{
        U[i][j] = 0;}

    }
}

}


void mjacobi(vector<double> a, vector<double> B, vector<double> XINIT, double tol, int it_max, vector<double> &SOL, int &niter, int &info){

SOL = XINIT;
int n = SOL.size();
int N = a.size();
vector<double> XTEMP(n);
double norm = 100;
niter = 0;

while(norm>tol && niter<it_max){
XTEMP[0] = (a[2*n-1]*SOL[1]+B[0])/a[n-1];
XTEMP[n-1] = (a[n-2]*SOL[n-1]+B[n-1])/a[2*n-2];
norm=0;
norm += abs(XTEMP[0]-SOL[0]);
norm+= abs(XTEMP[n-1] - SOL[n-1]);



for(int i=1;i<n-1;i++){
    XTEMP[i]=(a[i-1]*SOL[i-1]+a[(i-1)+2*n-1]*SOL[i+1]+B[i])/a[(i-1)+n-1];
    norm +=abs(XTEMP[i]-SOL[i]);
}
SOL = XTEMP;


niter++;
}



if(norm<= tol){info = 1;}

}

int main()
{

    vector<double> a (10);
    a = {-1,-1,-1,2,2,2,2,-1,-1,-1};
    vector<double> delt(4);
    mineurs(a,delt);

int N = a.size();
int n = (N+2)/3;


int i1 =0;
int i2=n-1;
int i3=2*n-1;

int i,j;

cout << "La matrice A est \n";
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            {if (i == j){
            cout<<a[i2]<<' ';
            i2++;
        }else if(i-j == 1){
        cout<<a[i1]<<' ';
        i1++;
        }else if(j-i == 1){
        cout<<a[i3]<<' ';
        i3++;
        }else{
        cout<<0<<' ';}


            }
        cout << "\n";
    }




cout<<"\n"<<"Les mineurs sont"<<"\n";

for(int i = 0; i<=n; i++){

    cout<<"Delt_"<<i<<" = "<<delt[i]<<"\n";

}
cout<<endl;


float L[4][4];
float U[4][4];

decLU(a,delt,L,U);



cout<<"A = LU";
cout<<endl;
cout<<endl;
    cout << "La matrice L est \n";
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            cout << L[i][j] << " ";
        cout << "\n";}
cout<<"\n";

        cout << "La matrice U est \n";
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            cout << U[i][j] << " ";
        cout << "\n";
    }

    cout<<"\n";


vector<double> B(3);
vector<double> XINIT(4);
vector<double> SOL(4);
XINIT = {0,0,0,0};
B = {13,12,-23,10};
double tol= 0.001;
int it_max = 1000;
int niter = 0; int info = 0;


mjacobi(a, B, XINIT, tol, it_max, SOL, niter, info);


cout<<endl<<endl<<"JACOBI"<<endl<<endl;
cout<<endl<<"Le vector B est"<<"("<<B[0]<<", "<<B[1]<<", "<<B[2]<<", "<<B[3]<<")"<<endl;;
cout<<"La solution converge (0= non, 1 = oui): "<<info<<endl;
cout<<"Le resultat de l'algo est ";
cout<<"("<<SOL[0]<<", "<<SOL[1]<<", "<<SOL[2]<<", "<<SOL[3]<<")";


    return 0;
}
