//Header files for  Implicit Method
#include <iostream>
#include<cmath>
#include<fstream>

//Defining Constants for entire Program
#define alpha 0.645     //Speed of bottom plate
#define b 3.5    //Distance between plates
#define h 3.5
#define time 0.4

const int nx1=36;   ///Nodes in x-direction = columns of the matrix
const int ny1=36;   ///Nodes in Y-direction = rows in the matrix

using namespace std;

/****************************************************************************************************
                                Function Prototypes Declared here
*****************************************************************************************************/

 void initialvalueout(double (&val)[ny1][nx1]);
 void TDMA_2D(double (&val)[ny1][nx1], double, double);
 void finalvalueout(double (&val)[ny1][nx1]);
 //void xsweep(double (&val)[ny1][nx1], double, double);
 //void ysweep(double (&val)[ny1][nx1], double, double);
//void resultout(double [], int);

/****************************************************************************************************
                                        Main Function here
*****************************************************************************************************/

int main()
{

    /// Declarations of Array Variables and their Definitions

    int nx, ny;    //Number of Grid Points in X and Y directions respectively.
    double dx,dy;  //Step Sizes in x and y directions.
    dx=0.1;        //Step Size in X-direction.
    dy=0.1;        //Step Size in Y-direction.
    nx=1+(b/dx);   //Number of grid points in X-direction.
    ny=1+(h/dy);   //Number of grid points in Y-direction.

    cout<<"The number of grid points in X-Direction with dx="<<dx<<" is: "<<nx<<endl;
    cout<<"The number of grid points in Y-Direction with dy="<<dy<<" is: "<<ny<<endl;

    double dt=0.01; //time step in hour.
    int nt;         //Number of time steps
    nt=1+(time/dt);

    cout<<"The number of time steps with dt="<<dt<<" are: "<<nt<<endl;

    double T[ny1][nx1];
    double d1,d2;   //Refer to notebook for these.

    d1=(0.5*alpha*dt)/pow(dx,2);
    d2=(0.5*alpha*dt)/pow(dy,2);

    cout<<"The value of d1 is: "<<d1<<endl;
    cout<<"The value of d2 is: "<<d2<<endl;


    ///Initialization of Variable Arrays

    for(int i=0;i<ny1;i++)  /// i is for rows
    {
        for(int j=0;j<nx1;j++)  /// j is for columns
        {
            if(j==0)
            {
                T[i][j]=200;
            }
            else if(j==ny1)
            {
                T[i][j]=0;
            }
            else if(i==0)
            {
                T[i][j]=200;
            }
            else if(i==nx1)
            {
                T[i][j]=0;
            }
            else
                T[i][j]=0;
        }
    }


    initialvalueout(T);   /// Function to generate the initial values at the grid
    TDMA_2D(T,d1,d2);        /// Function for TDMA
    finalvalueout(T);

    for(int i=0;i<nx1;i++)
    {
        for(int j=0;j<ny1;j++)
        {
            cout<<T[i][j];
        }
        cout<<endl;
    }

    return 0;

}

/****************************************************************************************************
                  Function for generating initial values at grid point to file
*****************************************************************************************************/

void initialvalueout(double (&val)[ny1][nx1])
{
        ofstream griddata("initial values.txt");
        if(griddata.is_open())
        {
             cout << "File Opened successfully!!!. Writing data from array to file" << endl;
             for(int j=nx1-1;j>=0;j--)
              {
                 for(int i=0;i<ny1;i++)
                 {
                      griddata<<val[i][j]<<"\t";
                 }
                      griddata<<"\n\n";
              }
        }
        else                  //file could not be opened
	    {
		cout << "File could not be opened." << endl;
	    }
}

/****************************************************************************************************
                                       TDMA 2D Subroutine
*****************************************************************************************************/

void TDMA_2D(double (&val)[ny1][nx1],double diff1,double  diff2)
{

    double T_old_x[ny1][nx1], T_old_y[ny1][nx1];    ///One old
    double a_x[nx1],b_x[nx1],c_x[nx1],D_x[nx1],G_x[nx1],H_x[nx1];   ///Constants for x-sweep
    double a_y[ny1],b_y[ny1],c_y[ny1],D_y[ny1],G_y[ny1],H_y[ny1];   ///Constants for y-sweep

    for(int i=0;i<nx1;i++)   ///definition of constants for X-Sweep.
    {
        a_x[i]=-diff1;
        b_x[i]=((2*diff1)+1);
        c_x[i]=-diff1;
    }

    for(int i=0;i<ny1;i++)  ///definition of constants for Y-Sweep.
    {
        a_y[i]=-diff2;
        b_y[i]=((2*diff2)+1);
        c_y[i]=-diff2;
    }

    for(int t=1;t<42;t++)  ///Big Time Loop
    {

        ///For X-Sweep

        for(int i=0;i<ny1;i++)
        {
            for(int j=0;j<nx1;j++)
            {
                T_old_x[i][j]=val[i][j];

                if(i==1)
                {
                    D_x[j]=(diff2*T_old_x[i][j+1])+((1-(2*diff2))*T_old_x[i][j])+(diff2*T_old_x[i][j-1])+(diff1*val[i-1][j]);
                }
                else if(i==34)
                {
                    D_x[j]=(diff2*T_old_x[i][j+1])+((1-(2*diff2))*T_old_x[i][j])+(diff2*T_old_x[i][j-1])+(diff1*val[i+1][j]);
                }
                else if(i!=1 && i!=34)
                {
                    D_x[j]=(diff2*T_old_x[i][j+1])+((1-(2*diff2))*T_old_x[i][j])+(diff2*T_old_x[i][j-1]);
                }

                if(j==0)
                {
                    G_x[j]=val[i][j];
                    H_x[j]=0;
                }
                else if(j!=0)
                {
                    G_x[j]=c_x[j]/(b_x[j]-(a_x[j]*H_x[j-1]));
                    H_x[j]=(D_x[j]-(a_x[j]*G_x[j-1]))/(b_x[j]-(a_x[j]*H_x[j-1]));
                }
            }
        }

        for(int i=1;i<ny1-1;i++)
        {
            for(int j=nx1-2;j>0;j--)
            {
                val[i][j]=(-H_x[j]*val[i][j+1])+G_x[j];
            }
        }

        ///For Y-Sweep

        for(int j=1;j<ny1;j++)
        {
            for(int i=1;i<nx1;i++)
            {
                T_old_y[i][j]=val[i][j];

                if(j==1)
                {
                   D_y[i]=(diff1*T_old_y[i+1][j])+((1-(2*diff1))*T_old_y[i][j])+(diff1*T_old_y[i-1][j])+(diff2*val[i][j-1]);
                }
                else if(j==34)
                {
                   D_y[i]=(diff1*T_old_y[i+1][j])+((1-(2*diff1))*T_old_y[i][j])+(diff1*T_old_y[i-1][j])+(diff2*val[i][j+1]);
                }
                else if(j!=1 && j!=34)
                {
                   D_y[i]=(diff1*T_old_y[i+1][j])+((1-(2*diff1))*T_old_y[i][j])+(diff1*T_old_y[i-1][j]);
                }

                if(i==0)
                {
                    G_y[i]=val[i][j];
                    H_y[i]=0;
                }
                else if(i!=0)
                {
                    G_y[i]=c_x[i]/(b_y[i]-(a_y[i]*H_y[i-1]));
                    H_y[i]=(D_y[i]-(a_y[i]*G_y[i-1]))/(b_y[i]-(a_y[i]*H_y[i-1]));
                }

            }
        }

        for(int j=1;j<nx1-1;j++)
        {
            for(int i=ny1-2;i>0;i--)
            {
                val[i][j]=(-H_y[i]*val[i+1][j])+G_y[i];
            }
        }


    }
}

/****************************************************************************************************
                                       X-Sweep Subroutine
*****************************************************************************************************/

//void xsweep(double (&val1)[nx1][ny1], double diff1, double diff2)
//{
//
//    double T_old[nx1][ny1];
//    double a_x[ny1],b_x[ny1],c_x[ny1],D_x[ny1],G_x[ny1],H_x[ny1];
//
//    for(int i=0;i<ny1;i++)
//    {
//        a_x[i]=-diff1;
//        b_x[i]=((2*diff1)+1);
//        c_x[i]=-diff1;
//    }
//
//
//    for(int i=0;i<nx1;i++)
//    {
//        for(int j=0;j<ny1;j++)
//        {
//            T_old[i][j]=val1[i][j];
//        }
//    }
//
//    for(int j=0;j<nx1;j++)
//    {
//        for(int i=0;i<ny1;i++)
//        {
//
//        if(i==1)
//        {
//            D_x[i]=diff2*T_old[i][j+1]+(1-2*diff2)*T_old[i][j]+diff2*T_old[i][j-1]+diff1*val1[i-1][j];
//        }
//
//        else if(i==34)
//        {
//            D_x[i]=diff2*T_old[i][j+1]+(1-2*diff2)*T_old[i][j]+diff2*T_old[i][j-1]+diff1*val1[i+1][j];
//        }
//
//        else
//        {
//            D_x[i]=diff2*T_old[i][j+1]+(1-2*diff2)*T_old[i][j]+diff2*T_old[i][j-1];
//        }
//
//        }
//    }
//
//    for(int j=0;j<nx1;j++)
//    {
//        for(int i=0;i<ny1;i++)
//        {
//            if(i==0)
//            {
//            G_x[i]=val1[i][0];
//            H_x[i]=0;
//            }
//            else
//            {
//            G_x[i]=c_x[i]/(b_x[i]-(a_x[i]*H_x[i-1]));
//            H_x[i]=(D_x[i]-(a_x[i]*G_x[i-1]))/(b_x[i]-(a_x[i]*H_x[i-1]));
//            }
//        }
//    }
//
//    for(int j=1;j<ny1-2;j++)
//    {
//        for(int i=nx1-2;i>0;i--)
//        {
//            val1[i][j]=(-H_x[i]*val1[i+1][j])+G_x[i];
//        }
//    }
//
//
//
//
//}
//
///****************************************************************************************************
//                                       Y-Sweep Subroutine
//*****************************************************************************************************/
//
////void ysweep(double (&val)[nx1][ny1], double diff1, double diff2)
////{
////
////}

/****************************************************************************************************
                  Function for generating final values at grid point to file
*****************************************************************************************************/

void finalvalueout(double (&val)[ny1][nx1])
{
        ofstream griddata("final values.txt");
        if(griddata.is_open())
        {
             cout << "File Opened successfully!!!. Writing data from array to file" << endl;
             for(int j=nx1-1;j>=0;j--)
              {
                 for(int i=0;i<ny1;i++)
                 {
                      griddata<<val[i][j]<<"\t";
                 }
                      griddata<<"\n\n";
              }
        }
        else                  //file could not be opened
	    {
		cout << "File could not be opened." << endl;
	    }
}
