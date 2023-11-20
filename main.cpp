#include <iostream>
 #include <cmath>
 #include <ctime>
 
 #define N 3
 
 using namespace std;
 //矩陣乘法
 double * mul(double A[N*N],double B[N*N])
 {
     double *C=new double[N*N]{};
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             for(int k=0;k<N;k++)
             {
                 C[i*N+j] += A[i*N+k]*B[k*N+j];
             }
         }
     }
 
     //若絕對值小於10^-10,則置為0（這是我自己定的）
     for(int i=0;i<N*N;i++)
     {
         if(abs(C[i])<pow(10,-10))
         {
             C[i]=0;
         }
     }
 
     return C;
 }
 
 //LUP分解
 void LUP_Descomposition(double A[N*N],double L[N*N],double U[N*N],int P[N])
 {
     int row=0;
     for(int i=0;i<N;i++)
     {
         P[i]=i;
     }
     for(int i=0;i<N-1;i++)
     {
         double p=0.0d;
         for(int j=i;j<N;j++)
         {
             if(abs(A[j*N+i])>p)
             {
                 p=abs(A[j*N+i]);
                 row=j;
             }
         }
         if(0==p)
         {
             cout<< "矩陣奇異，無法計算逆" <<endl;
             return ;
         }
 
         //交換P[i]和P[row]
         int tmp=P[i];
         P[i]=P[row];
         P[row]=tmp;
 
         double tmp2=0.0d;
         for(int j=0;j<N;j++)
         {
             //交換A[i][j]和 A[row][j]
             tmp2=A[i*N+j];
             A[i*N+j]=A[row*N+j];
             A[row*N+j]=tmp2;
         }
 
         //以下同LU分解
         double u=A[i*N+i],l=0.0d;
         for(int j=i+1;j<N;j++)
         {
             l=A[j*N+i]/u;
             A[j*N+i]=l;
             for(int k=i+1;k<N;k++)
             {
                 A[j*N+k]=A[j*N+k]-A[i*N+k]*l;
             }
         }
 
     }
 
     //構造L和U
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<=i;j++)
         {
             if(i!=j)
             {
                 L[i*N+j]=A[i*N+j];
             }
             else
             {
                 L[i*N+j]=1;
             }
         }
         for(int k=i;k<N;k++)
         {
             U[i*N+k]=A[i*N+k];
         }
     }
 
 }
 
 //LUP求解方程
 double * LUP_Solve(double L[N*N],double U[N*N],int P[N],double b[N])
 {
     double *x=new double[N]();
     double *y=new double[N]();
 
     //正向替換
     for(int i = 0;i < N;i++)
     {
         y[i] = b[P[i]];
         for(int j = 0;j < i;j++)
         {
             y[i] = y[i] - L[i*N+j]*y[j];
         }
     }
     //反向替換
     for(int i = N-1;i >= 0; i--)
     {
         x[i]=y[i];
         for(int j = N-1;j > i;j--)
         {
             x[i] = x[i] - U[i*N+j]*x[j];
         }
         x[i] /= U[i*N+i];
     }
     return x;
 }
 
 /*****************矩陣原地轉置BEGIN********************/
 
 /* 後繼 */
 int getNext(int i, int m, int n)
 {
   return (i%n)*m + i/n;
 }
 
 /* 前驅 */
 int getPre(int i, int m, int n)
 {
   return (i%m)*n + i/m;
 }
 
 /* 處理以下標i為起點的環 */
 void movedata(double *mtx, int i, int m, int n)
 {
   double temp = mtx[i]; // 暫存
   int cur = i;    // 當前下標
   int pre = getPre(cur, m, n);
   while(pre != i)
   {
     mtx[cur] = mtx[pre];
     cur = pre;
     pre = getPre(cur, m, n);
   }
   mtx[cur] = temp;
 }
 
 /* 轉置，即迴圈處理所有環 */
 void transpose(double *mtx, int m, int n)
 {
   for(int i=0; i<m*n; ++i)
   {
     int next = getNext(i, m, n);
     while(next > i) // 若存在後繼小於i說明重複,就不進行下去了（只有不重複時進入while迴圈）
       next = getNext(next, m, n);
     if(next == i)  // 處理當前環
       movedata(mtx, i, m, n);
   }
 }
 /*****************矩陣原地轉置END********************/
 
 //LUP求逆(將每列b求出的各列x進行組裝)
 double * LUP_solve_inverse(double A[N*N])
 {
     //創建矩陣A的副本，注意不能直接用A計算，因為LUP分解演算法已將其改變
     double *A_mirror = new double[N*N]();
     double *inv_A=new double[N*N]();//最終的逆矩陣（還需要轉置）
     double *inv_A_each=new double[N]();//矩陣逆的各列
     //double *B    =new double[N*N]();
     double *b    =new double[N]();//b陣為B陣的列矩陣分量
 
     for(int i=0;i<N;i++)
     {
         double *L=new double[N*N]();
         double *U=new double[N*N]();
         int *P=new int[N]();
 
         //構造單位陣的每一列
         for(int i=0;i<N;i++)
         {
             b[i]=0;
         }
         b[i]=1;
 
         //每次都需要重新將A複製一份
         for(int i=0;i<N*N;i++)
         {
             A_mirror[i]=A[i];
         }
 
         LUP_Descomposition(A_mirror,L,U,P);
 
         inv_A_each=LUP_Solve (L,U,P,b);
         memcpy(inv_A+i*N,inv_A_each,N*sizeof(double));//將各列拼接起來
     }
     transpose(inv_A,N,N);//由於現在根據每列b算出的x按行存儲，因此需轉置
 
     return inv_A;
 }
 
 int main()
 {
     double *A = new double[N*N]();
 
    /*
     srand((unsigned)time(0));
     for(int i=0; i<N ;i++)
     {
         for(int j=0; j<N;j++)
         {
             A[i*N+j]=rand()%100 *0.01;
         }
     }
     */
     
     A[0] = 3;
     A[1] = 471;
     A[2] = 73949;
     A[3] = 471;
     A[4] = 73949;
     A[5] = 11610621;
     A[6] = 73949;
     A[7] = 11610621;
     A[8] = 1823015393;
 
 
     double *E_test1 = new double[N*N]();
     double *E_test2 = new double[N*N]();
     double *invOfA = new double[N*N]();
     invOfA=LUP_solve_inverse(A);
 
     E_test1=mul(A,invOfA);    //驗證精確度
     E_test2=mul(invOfA,A);    //驗證精確度
 
     cout<< "Matrix A:" <<endl;
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             cout<< A[i*N+j]<< " " ;
         }
         cout<<endl;
     }
     
     cout<<endl;
 
     cout<< "inv_A:" <<endl;
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             cout<< invOfA[i*N+j]<< " " ;
         }
         cout<<endl;
     }
     
     cout<<endl;
 
     cout<< "E_test1:" <<endl;    
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             cout<< E_test1[i*N+j]<< " " ;
         }
         cout<<endl;
     }
     
     cout<<endl;

     cout<< "E_test2:" <<endl;    
     for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             cout<< E_test2[i*N+j]<< " " ;
         }
         cout<<endl;
     }
 
     return 0;
 }
