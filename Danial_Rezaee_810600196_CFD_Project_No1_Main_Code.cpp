//Danial Rezaee-810600196-March 2022-Computational Fluid Dynamics Project No. 1-Dr. Nejat
#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<array>
#include<cstring>
#include<string>
#include<ctime>
#include<algorithm>

using namespace std;



double f(double );
double BC_left(int , int , int );
double BC_right(int, int , int );
double BC_south(int , int , int );
double BC_north(int , int , int );
double exact(double ,double );
void initialize(vector<vector<double>> &);
void Gauss(vector<vector<double>> &,int );
double L=1.0;
int main()
{

int i,j;
int n;


int m=5000;
int WS,MS,MSG;
int Meshes[]={10,20,40};
int MeshG[]={10,20,40};
long double s1,s2,L1,L2,Linf,EE,w;
long double W[]={0,1,1.1,1.3,1.5};




MS=(sizeof(Meshes)/ sizeof(*Meshes));
WS=(sizeof(W)/ sizeof(*W));
MSG=(sizeof(MeshG)/ sizeof(*MeshG));

sort(Meshes,Meshes+MS);
sort(MeshG,MeshG+MSG);
sort(W,W+WS);


int cc=MS*(WS)*6;
int dd=MS*WS;
vector<vector<double>>T(n,vector<double>(n,0));
vector<vector<double>>E(n,vector<double>(n,0));
vector<vector<double>>R(m,vector<double>(cc,0));
vector<vector<double>>SOR(n,vector<double>(n,0));
vector<vector<double>>J(n,vector<double>(n,0));
vector<vector<double>>TAB(dd,vector<double>(7,0));


ofstream Length;
Length.open("L.csv");
Length<<L<<endl;

ofstream NoM;
NoM.open("NoM.csv");
for(int i=0;i<MS;i++)
NoM<<Meshes[i]<<endl;

ofstream NoMG;
NoMG.open("NoMG.csv");
for(int i=0;i<MS;i++)
NoMG<<MeshG[i]<<endl;

ofstream NoW;
NoW.open("NoW.csv");
for(int i=0;i<WS;i++)
NoW<<W[i]<<endl;


ofstream FResult;
FResult.open("FResult.csv");


ofstream DResult;
DResult.open("DResult.csv");

ofstream Error;
Error.open("Error.csv");

ofstream Tempreture;
Tempreture.open("Tempreture.csv");

ofstream Table;
Table.open("Table.csv");


long double Del=10.0;
long double tt=pow(10,-10);
int k=0;
int ci=0;

ofstream Tolerance;
Tolerance.open("Tolerance.txt");
Tolerance<<tt;

ofstream N;
N.open("N.txt");
N<<Meshes[MS-1]*Meshes[MS-1];

for(int ww=0;ww<MS;ww++){
n=Meshes[ww];
E.clear();
T.clear();
SOR.clear();
J.clear();
E.resize(n,vector<double>(n,0));
T.resize(n,vector<double>(n,0));
SOR.resize(n,vector<double>(n,0));
J.resize(n,vector<double>(n,0));





for(int vv=0;vv<WS;vv++){
k=0;	
w=W[vv];
Del=10;
initialize(T);
initialize(E);
initialize(SOR);
initialize(J);

clock_t c_start = clock();
while(Del>=tt){
k=k+1;	
Linf=0;	
s2=0;
s1=0;
for(i=0;i<n;i++){
for(j=0;j<n;j++){
	

if(i==0 && j==0 && w==0 ){
T[0][0]=1/6.0*(J[1][0]+J[0][1]+2*BC_left(i,j,n)+2.0*BC_south(i,j,n));

}else if(i==0 && j !=0 && j<n-1 && w ==0){
T[0][j]=1/5.0*(J[1][j]+J[0][j+1]+J[0][j-1]+2.0*BC_left(i,j,n));}

else if(i==0 && j==n-1 && w ==0){
T[0][n-1]=1/6.0*(J[1][n-1]+J[0][n-2]+2.0*BC_north(i,j,n)+2.0*BC_left(i,j,n));}

else if(j==n-1 && i !=0 && i<n-1 && w ==0){
T[i][n-1]=(1/5.0)*(J[i+1][n-1]+J[i-1][n-1]+2.0*BC_north(i,j,n)+J[i][n-2]);}

else if(j==n-1 && i==n-1 && w ==0){
T[n-1][n-1]=(1/6.0)*(J[n-2][n-1]+J[n-1][n-2]+2.0*BC_north(i,j,n)+2.0*BC_right(i,j,n));}

else if(i==n-1 && j<n-1 && j !=0 && w ==0){
T[n-1][j]=1/5.0*(J[n-2][j]+J[n-1][j+1]+J[n-1][j-1]+2.0*BC_right(i,j,n));}

else if(i==n-1 && j==0 && w ==0){
T[n-1][0]=1/6.0*(J[n-2][j]+J[n-1][j+1]+2.0*BC_right(i,j,n)+2.0*BC_south(i,j,n));}

else if(j==0 && i !=0 && i < n-1 && w ==0){
T[i][0]=1/5.0*(J[i][1]+J[i+1][0]+J[i-1][0]+2.0*BC_south(i,j,n));}

else if(w ==0) {
T[i][j]=1/4.0*(J[i+1][j]+J[i-1][j]+J[i][j+1]+J[i][j-1]);}	

	
if(i==0 && j==0 && w !=0){
T[0][0]=1/6.0*(T[1][0]+T[0][1]+2*BC_left(i,j,n)+2.0*BC_south(i,j,n));

}else if(i==0 && j !=0 && j<n-1 && w !=0){
T[0][j]=1/5.0*(T[1][j]+T[0][j+1]+T[0][j-1]+2.0*BC_left(i,j,n));}

else if(i==0 && j==n-1 && w !=0){
T[0][n-1]=1/6.0*(T[1][n-1]+T[0][n-2]+2.0*BC_north(i,j,n)+2.0*BC_left(i,j,n));}

else if(j==n-1 && i !=0 && i<n-1 && w !=0){
T[i][n-1]=(1/5.0)*(T[i+1][n-1]+T[i-1][n-1]+2.0*BC_north(i,j,n)+T[i][n-2]);}

else if(j==n-1 && i==n-1 && w !=0){
T[n-1][n-1]=(1/6.0)*(T[n-2][n-1]+T[n-1][n-2]+2.0*BC_north(i,j,n)+2.0*BC_right(i,j,n));}

else if(i==n-1 && j<n-1 && j !=0 && w !=0){
T[n-1][j]=1/5.0*(T[n-2][j]+T[n-1][j+1]+T[n-1][j-1]+2.0*BC_right(i,j,n));}

else if(i==n-1 && j==0 && w !=0){
T[n-1][0]=1/6.0*(T[n-2][j]+T[n-1][j+1]+2.0*BC_right(i,j,n)+2.0*BC_south(i,j,n));}

else if(j==0 && i !=0 && i < n-1 && w !=0){
T[i][0]=1/5.0*(T[i][1]+T[i+1][0]+T[i-1][0]+2.0*BC_south(i,j,n));}

else if(w !=0) {
T[i][j]=1/4.0*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]);}


E[i][j]=abs(T[i][j]-(n*n/(L*L))*((-cos(M_PI*L*i/n)+cos(M_PI*L*(i+1)/n))*(cosh(M_PI*L*j/n)-cosh(M_PI*L*(j+1)/n)))/(M_PI*M_PI*sinh(M_PI)));
s2=s2+(L*L/(n*n))*E[i][j]*E[i][j];
s1=s1+(L*L/(n*n))*E[i][j];
if(E[i][j]>Linf)
Linf=E[i][j];

if(w !=0){
T[i][j]=SOR[i][j]+w*(T[i][j]-SOR[i][j]);
SOR[i][j]=T[i][j];}



}
}
if(w==0)
J=T;
L1=s1/(L*L);
L2=sqrt(s2/(L*L));
R[k-1][ci*6]=n;
R[k-1][ci*6+1]=k;
R[k-1][ci*6+2]=w;
R[k-1][ci*6+3]=L1;
R[k-1][ci*6+4]=L2;
R[k-1][ci*6+5]=Linf;

if(k>=2)
Del=abs(R[k-1][ci*6+4]-R[k-2][ci*6+4]);

}
clock_t c_end=clock();
long double time_elapsed_ms=1000.0*(c_end-c_start)/CLOCKS_PER_SEC;
if(w==0){
cout<<"Jacobi"<<"\t"<<"N="<<n*n<<"\t"<<"i="<<k<<"\t"<<"L1="<<L1<<"\t"<<"L2="<<L2<<"\t"<<"Linf="<<Linf<<"\t"<<"CPU time="<<time_elapsed_ms<<endl;
FResult<<"Jacobi"<<","<<n*n<<","<<w<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
Table<<"Jacobi"<<","<<n*n<<","<<"-"<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
}
else if(w==1){
cout<<"Gauss-Sidel"<<"\t"<<"N="<<n*n<<"\t"<<"i="<<k<<"\t"<<"L1="<<L1<<"\t"<<"L2="<<L2<<"\t"<<"Linf="<<Linf<<"\t"<<"CPU time="<<time_elapsed_ms<<endl;
FResult<<"Gauss-Sidel"<<","<<n*n<<","<<w<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
Table<<"Gauss-Sidel"<<","<<n*n<<","<<"-"<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;

}
else{
cout<<"SOR"<<"\t"<<"N="<<n*n<<"\t"<<"w="<<w<<"\t"<<"i="<<k<<"\t"<<"L1="<<L1<<"\t"<<"L2="<<L2<<"\t"<<"Linf="<<Linf<<"\t"<<"CPU time="<<time_elapsed_ms<<endl;
FResult<<"SOR"<<","<<n*n<<","<<w<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
Table<<"SOR"<<","<<n*n<<","<<w<<","<<k<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
}	
ci=ci+1;
}
}

for(i=0;i<m;i++){
for(j=0;j<cc;j++){
DResult<<R[i][j]<<",";
}
DResult<<endl;
}

for(i=0;i<n;i++){
for(j=0;j<n;j++){
Error<<E[i][j]<<",";
Tempreture<<T[i][j]<<",";
}
Error<<endl;
Tempreture<<endl;
}


vector<vector<double>>G(n*n,vector<double>(n*n+1,0));

for(int ww=0;ww<MSG;ww++){
n=MeshG[ww];
G.clear();
G.resize(n*n,vector<double>(n*n+1,0));
E.clear();
E.resize(n,vector<double>(n,0));
clock_t c_start = clock();
for(i=0;i<n;i++){
for(j=0;j<n;j++){
	
	
if(i==0 && j==0 ){
G[n*i+j][n*i+j]=6;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*i+j+1]=-1;

G[n*i+j][n*n]=2*BC_left(i,j,n)+2.0*BC_south(i,j,n);

}else if(i==0 && j !=0 && j<n-1){
G[n*i+j][n*i+j]=5;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*i+j+1]=-1;	
G[n*i+j][n*i+j-1]=-1;		
G[n*i+j][n*n]=2.0*BC_left(i,j,n);}

else if(i==0 && j==n-1 ){

G[n*i+j][n*i+j]=6;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*i+j-1]=-1;

G[n*i+j][n*n]=2.0*BC_north(i,j,n)+2.0*BC_left(i,j,n);}

else if(j==n-1 && i !=0 && i<n-1 ){

G[n*i+j][n*i+j]=5;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*i+j-1]=-1;
G[n*i+j][n*(i-1)+j]=-1;		
G[n*i+j][n*n]=2.0*BC_north(i,j,n);

}

else if(j==n-1 && i==n-1 ){
T[n-1][n-1]=(1/6.0)*(J[n-2][n-1]+J[n-1][n-2]+2.0*BC_north(i,j,n)+2.0*BC_right(i,j,n));


G[n*i+j][n*i+j]=6;
G[n*i+j][n*(i-1)+j]=-1;
G[n*i+j][n*i+j-1]=-1;
G[n*i+j][n*n]=2.0*BC_north(i,j,n)+2.0*BC_right(i,j,n);

}

else if(i==n-1 && j<n-1 && j !=0 ){

G[n*i+j][n*i+j]=5;
G[n*i+j][n*(i-1)+j]=-1;
G[n*i+j][n*i+j-1]=-1;
G[n*i+j][n*i+j+1]=-1;		
G[n*i+j][n*n]=2.0*BC_right(i,j,n);

}

else if(i==n-1 && j==0 ){


G[n*i+j][n*i+j]=6;
G[n*i+j][n*(i-1)+j]=-1;
G[n*i+j][n*i+j+1]=-1;
G[n*i+j][n*n]=2.0*BC_right(i,j,n)+2.0*BC_south(i,j,n);

}

else if(j==0 && i !=0 && i < n-1 ){

G[n*i+j][n*i+j]=5;
G[n*i+j][n*(i-1)+j]=-1;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*i+j+1]=-1;		
G[n*i+j][n*n]=2.0*BC_south(i,j,n);


}

else {
G[n*i+j][n*i+j]=4;
G[n*i+j][n*(i+1)+j]=-1;
G[n*i+j][n*(i-1)+j]=-1;
G[n*i+j][n*i+j+1]=-1;
G[n*i+j][n*i+j-1]=-1;



}	
}
}

Linf=0;
s1=0;
s2=0;
Gauss(G,n);

for(i=0;i<n;i++){
for(j=0;j<n;j++){

E[i][j]=abs(G[n*i+j][n*n]/G[n*i+j][n*i+j]-(n*n/(L*L))*((-cos(M_PI*L*i/n)+cos(M_PI*L*(i+1)/n))*(cosh(M_PI*L*j/n)-cosh(M_PI*L*(j+1)/n)))/(M_PI*M_PI*sinh(M_PI)));
s2=s2+(L*L/(n*n))*E[i][j]*E[i][j];
s1=s1+(L*L/(n*n))*E[i][j];
if(E[i][j]>Linf)
Linf=E[i][j];
}
}
L1=s1/(L*L);
L2=sqrt(s2/(L*L));
clock_t c_end=clock();
long double time_elapsed_ms=1000.0*(c_end-c_start)/CLOCKS_PER_SEC;
cout<<"Gauss Elimination"<<"\t"<<"N="<<n*n<<"\t"<<"L1="<<L1<<"\t"<<"L2="<<L2<<"\t"<<"Linf="<<Linf<<"\t"<<"CPU time="<<time_elapsed_ms<<endl;
FResult<<"Gauss Elimination"<<","<<n*n<<","<<"-"<<","<<"-"<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
Table<<"Gauss Elimination"<<","<<n*n<<","<<"-"<<","<<"-"<<","<<L1<<","<<L2<<","<<Linf<<","<<time_elapsed_ms<<endl;
}
ofstream Gauss;
Gauss.open("G.csv");


ofstream TG;
TG.open("TG.csv");


for(i=0;i<n*n;i++){
for(j=0;j<(n*n+1);j++){

Gauss<<G[i][j]<<",";
}
Gauss<<endl;
}

for(i=0;i<n;i++){
for(j=0;j<n;j++){
TG<<G[n*i+j][n*n]/G[n*i+j][n*i+j]<<",";
}
TG<<endl;
}




system("pause")	;
return 0;
}
double f(double x)
{
return sin(M_PI*x);
}
double BC_left(int i, int j, int n)
{
return 0;
}
double BC_right(int i, int j, int n)
{
return 0;
}
double BC_south(int i, int j, int n)
{
return 0;
}
double BC_north(int i, int j, int n)
{
return f(L*(2*i+1)/(2*n));
 
}
double exact(double x,double y)
{
return sin(M_PI*x)*sinh(M_PI*y)/sinh(M_PI);
}
void initialize(vector<vector<double>> &A)
{
for(int i=0;i !=A.size();i++){
for(int j=0;j !=A.size();j++){
A[i][j]=0;
}
}
}
void Gauss(vector<vector<double>> &G,int n)
{
int i,j,k;
double co;
	
for(i=0;i<n*n-1;i++){
for(j=i+1;j<n*n;j++){
co=-G[j][i]/G[i][i];
for(k=0;k<n*n+1;k++){
G[j][k]=co*G[i][k]+G[j][k];

}
}
}

i=0;
j=0;
k=0;
for(i=n*n-1;i>0;i--){
for(j=i-1;j>=0;j--){
co=-G[j][i]/G[i][i];
for(k=n*n;k>=0;k--){
G[j][k]=co*G[i][k]+G[j][k];
}	
}
}

}

	




















