#include <iostream>
#include<cmath>
#include<fstream>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>     // usleep

# define l 400 //Defining Constants could have been done with static constants too
# define dx 5
#define dt1 0.02
#define dt2 0.01
#define dt3 0.005
#define a 250
#define pi 3.142857143
#define time1 0.50


const double c1=(a*dt1)/dx;  ///Courant Number for case 1
const double nt1=1+(time1/dt1);

const double c2=(a*dt2)/dx;  ///Courant Number for case 2
const double nt2=1+(time1/dt2);

const double c3=(a*dt3)/dx;  ///Courant Number for case 3
const double nt3=1+(time1/dt3);

using namespace std;

int main()
{
    cout<<"***********Program for solving the first order wave equation***********"<<endl;
    cout<<"\n\n";
    cout<<"The speed of the wave is:"<<a<<endl;
    cout<<"The spatial step size is:"<<dx<<endl;
    cout<<"The time step for case 1 is:"<<dt1<<endl;
    cout<<"The courant number for case 1 is:"<<c1<<endl;
    cout<<"The time step for case 2 is:"<<dt2<<endl;
    cout<<"The courant number for case 2 is:"<<c2<<endl;
    cout<<"The time step for case 3 is:"<<dt3<<endl;
    cout<<"The courant number for case 3 is:"<<c3<<endl;
    cout<<"The length of the tube is:"<<l<<endl;
    cout<<"The value of PI is:"<<pi<<endl;
    cout<<"The total number of time steps for case 1 are:"<<nt1<<endl;
    cout<<"The total number of time steps for case 2 are:"<<nt2<<endl;
    cout<<"The total number of time steps for case 3 are:"<<nt3<<endl;
//    double d=sin();
//
//    cout<<d<<endl;

    int grid=1+l/dx;

    cout<<"The size of matrix is:"<<grid<<endl;

    double x[grid];   ///Grid matrix
    double u1[grid];   ///Velocity Matrix
    double u2[grid];
    double u3[grid];
    /// Code for gnuplot pipe
    FILE *pipe1 = popen("gnuplot -persist", "w");
    FILE *pipe2 = popen("gnuplot -persist", "w");
    FILE *pipe3 = popen("gnuplot -persist", "w");

    /// set axis ranges
    fprintf(pipe1,"set xrange [0:400]\n");
    fprintf(pipe1,"set yrange [0:150]\n");
    fprintf(pipe2,"set xrange [0:400]\n");
    fprintf(pipe2,"set yrange [0:150]\n");
    fprintf(pipe3,"set xrange [0:400]\n");
    fprintf(pipe3,"set yrange [0:150]\n");
    //fprintf(pipe, "set xlabel aa");
    ///Initializing the grid matrices

    x[0]=0;
    for(int i=1;i<grid;i++)
    {
        x[i]=x[i-1]+dx;
    }

    ///Outputting the grid matrix

    ofstream griddata("gridpoints.dat");
    if(griddata.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
        griddata<<"Point"<<"\t"<<"location"<<endl;
    for(int i=0;i<grid;i++)
    {
        griddata<<i<<"\t\t"<<x[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}




    ///Initializing the velocity matrix with the given IC's and BC's

    u1[0]=0;
    u1[grid-1]=0;

   for(int i=1;i<grid-1;i++)
    {
        if(i>=1 && i<10)
            u1[i]=0;
        else if(i>=10 && i<=22)
            {
                double f=i;
                double g=f*5;
                double h=(g-50)/60;
                double j=pi*h;
                double k=sin(j);
                u1[i]=100*k;
            }
        else if(i>22)
               u1[i]=0;

        if(u1[i]<0)
                u1[i]=0;
    }

    for(int i=0;i<grid;i++)
    {
        u2[i]=u1[i];
        u3[i]=u1[i];
    }

   ///Outputting the velocity matrix of u1

    ofstream griddata1("initial data_u1.dat");
    if(griddata1.is_open())
    {
        griddata1<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata1<<i<<"\t\t"<<x[i]<<"\t\t"<<u1[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

	///Outputting the velocity matrix of u2

    ofstream griddata2("initial data_u2.dat");
    if(griddata2.is_open())
    {
        griddata2<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata2<<i<<"\t\t"<<x[i]<<"\t\t"<<u2[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

	///Outputting the velocity matrix of u3

    ofstream griddata3("initial data_u3.dat");
    if(griddata3.is_open())
    {
        griddata3<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata3<<i<<"\t\t"<<x[i]<<"\t\t"<<u3[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}



    double u_old1[grid];
    double u_old2[grid];
    double u_old3[grid];
    ///Time loop 1----Solving the equation----for solving case 1 with c1

    for(int t=1;t<27;t++)
    {

        fprintf(pipe1,"plot '-' using 1:2 with linespoints \n");
        for (int m=0;m<grid;m++)  // 10 datapoints per plot
        {
        fprintf(pipe1, "%lf %lf\n",x[m],u1[m]);  // passing x,y data pairs one at a time to gnuplot
        }


        for(int i=0;i<grid;i++)
            u_old1[i]=u1[i];

        for(int i=1;i<grid-1;i++)
        {
            u1[i]=u_old1[i]-c1*(u_old1[i]-u_old1[i-1]);
        }
        fprintf(pipe1, "set terminal gif animate\n");
        fprintf(pipe1, "set output 'c1.gif'\n");
        fprintf(pipe1, "replot\n");
        fprintf(pipe1,"e \n");    // finally, e
        fflush(pipe1);   // flush the pipe to update the plot
        usleep(100000);// wait a second before updating again
        //fprintf(pipe, "replot\n");
    }

    ///Time loop 2----Solving the equation----for solving case 2 with c2
    for(int t=1;t<52;t++)
    {

        fprintf(pipe2,"plot '-' using 1:2 with linespoints \n");
        for (int m=0;m<grid;m++)  // 10 datapoints per plot
        {
        fprintf(pipe2, "%lf %lf\n",x[m],u2[m]);  // passing x,y data pairs one at a time to gnuplot
        }


        for(int i=0;i<grid;i++)
            u_old2[i]=u2[i];

        for(int i=1;i<grid-1;i++)
        {
            u2[i]=u_old2[i]-c2*(u_old2[i]-u_old2[i-1]);
        }
        fprintf(pipe2, "set terminal gif animate\n");
        fprintf(pipe2, "set output 'c2.gif'\n");
        fprintf(pipe2, "replot\n");
        fprintf(pipe2,"e \n");    // finally, e
        fflush(pipe2);   // flush the pipe to update the plot
        usleep(100000);// wait a second before updating again
    }

    ///Time loop 3----Solving the equation----for solving case 3 with c3
    for(int t=1;t<102;t++)
    {

        fprintf(pipe3,"plot '-' using 1:2 with linespoints \n");
        for (int m=0;m<grid;m++)  // 10 datapoints per plot
        {
        fprintf(pipe3, "%lf %lf\n",x[m],u3[m]);  // passing x,y data pairs one at a time to gnuplot
        }


        for(int i=0;i<grid;i++)
            u_old3[i]=u3[i];

        for(int i=1;i<grid-1;i++)
        {
            u3[i]=u_old3[i]-c3*(u_old3[i]-u_old3[i-1]);
        }
        fprintf(pipe3, "set terminal gif animate\n");
        fprintf(pipe3, "set output 'c3.gif'\n");
       // fprintf(pipe, "replot\n");
        fprintf(pipe3,"e \n");    // finally, e
        fflush(pipe3);   // flush the pipe to update the plot
        usleep(100000);// wait a second before updating again
        //fprintf(pipe, "replot\n");
    }


    fclose(pipe1);
    fclose(pipe2);
    fclose(pipe3);

    ofstream griddata4("final data_u1.dat");
    if(griddata4.is_open())
    {
        griddata4<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata4<<i<<"\t\t"<<x[i]<<"\t\t"<<u1[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

	ofstream griddata5("final data_u2.dat");
    if(griddata5.is_open())
    {
        griddata5<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata5<<i<<"\t\t"<<x[i]<<"\t\t"<<u2[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

	ofstream griddata6("final data_u3.dat");
    if(griddata6.is_open())
    {
        griddata4<<"grid"<<"\t"<<"point"<<"\t"<<"Velocity"<<endl;
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(int i=0;i<grid;i++)
    {
        griddata6<<i<<"\t\t"<<x[i]<<"\t\t"<<u3[i]<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

    return 0;
}
