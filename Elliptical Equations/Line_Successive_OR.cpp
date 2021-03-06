///Program for Solution of Steady heat equation using Line Successive Over Relaxation Method.
//Header files for  Implicit Method
#include <iostream>
#include<cmath>
#include<fstream>

//Defining Constants for entire Program
#define l 1          //width of plate in x-direction
#define h 2           //height of plate in y-direction
#define dx 0.05
#define dy 0.05
#define w_opt 1.265

const double beta=dx/dy;

using namespace std;

int main()
{

    int rows,cols;
    rows=1+(h/dy);
    cols=1+(l/dx);

    cout<<"The number of rows are (JM):"<<rows<<endl;
    cout<<"The number of columns are (IM):"<<cols<<endl;
    cout<<"The value of beta is:"<<beta<<endl;

    double errormax;
    double T1,T2,T3,T4;
    T1=100;
    T2=0;
    T3=0;
    T4=0;

    double T[rows][cols],Tstar[rows][cols];
    double error[rows][cols];

    ///initializing the T matrix.

      for(int i=0;i<cols;i++)    ///zeroth row all columns
           T[0][i]=T1;

      for(int j=1;j<rows;j++)    ///zeroth column all rows
           T[j][0]=T2;

      for(int i=1;i<cols;i++)    ///last row all columns
           T[rows-1][i]=T3;

      for(int j=1;j<rows;j++)    ///last column all rows
           T[j][cols-1]=T4;


      //interior nodes set to zero

  for(int j=1;j<rows-1;j++)   ///For the rows
  {
      for(int i=1;i<cols-1;i++)  ///For the columns
      {
          T[j][i]=0;
      }
  }

      ///For outputting initial conditions to a file

        ofstream griddata("initial values.dat");
        if(griddata.is_open())
        {
             cout << "File Opened successfully!!!. Writing data from array to file" << endl;
             for(int j=rows-1;j>=0;j--)   ///for the rows
              {
                 for(int i=0;i<cols;i++)   ///for the columns
                 {
                      //if((i==0)||(i==5)||(i==10)||(i==15)||(i==20)||(i==25)||(i==30)||(i==35))
                      griddata<<T[j][i]<<"\t";
                 }
                      griddata<<"\n";
              }
        }
        else                  //file could not be opened
	    {
		cout << "File could not be opened." << endl;
	    }

double A1=1;
double A2=2*(1+pow(beta,2));
double A=A1/A2;

double a[rows],b[rows],c[rows],D[rows];
double H[rows],G[rows];

for(int j=0;j<rows;j++)
{
    a[j]=pow(beta,2);
    b[j]=-2*(1+pow(beta,2));
    c[j]=pow(beta,2);
}

H[0]=0;
for(int j=1;j<rows;j++)
{
    H[j]=c[j]/(b[j]-(a[j]*H[j-1]));
}


int k=0; ///counter for convergence

do
{
      for(int i=0;i<cols;i++)
      {
        for(int j=0;j<rows;j++)
        {
          Tstar[j][i]=T[j][i];
        }
      }

      for(int i=1;i<cols-1;i++)
      {
          for(int j=1;j<rows-1;j++)
          {
              D[j]=-Tstar[j][i+1]-T[j][i-1];

              G[0]=T1;
              G[j]=(D[j]-(a[j]*G[j-1]))/(b[j]-(a[j]*H[j-1]));
          }

          for(int j=rows-2;j>0;j--)
          {
              T[j][i]=(-H[j]*T[j+1][i])+G[j];
              T[j][i]=(1-w_opt)*Tstar[j][i]+w_opt*T[j][i];
          }
      }

//      for(int i=1;i<cols-1;i++)
//        for(int j=1;j<rows-1;j++)
//           T[j][i]=(1-w_opt)*Tstar[j][i]+w_opt*T[j][i];


   //}
      for(int i=0;i<cols;i++)
      {
          for(int j=0;j<rows;j++)
          {
              error[j][i]=abs(T[j][i]-Tstar[j][i]);
          }
      }

      errormax=error[0][0];
      for(int i=0;i<cols;i++)
      {
          for(int j=0;j<rows;j++)
          {
              if(error[j][i]>errormax)
                errormax=error[j][i];
          }
      }
   k++;
}while(errormax>1e-4);


    cout<<"The solution converges in"<<" "<<k<<" "<<"iterations"<<endl;

      ///For outputting initial conditions to a file

        ofstream griddata1("final values.dat");
        if(griddata.is_open())
        {
             cout << "File Opened successfully!!!. Writing data from array to file" << endl;
             for(int j=rows-1;j>=0;j--)   ///for the rows
              {
                 for(int i=0;i<cols;i++)   ///for the columns
                 {
                      //if((i==0)||(i==5)||(i==10)||(i==15)||(i==20)||(i==25)||(i==30)||(i==35))
                      griddata1<<T[j][i]<<"\t";
                 }
                      griddata1<<"\n";
              }
        }
        else                  //file could not be opened
	    {
		cout << "File could not be opened." << endl;
	    }




return 0;
}
