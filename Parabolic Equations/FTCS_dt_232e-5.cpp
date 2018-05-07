#include <iostream>
#include<cmath>
#include<fstream>
# define h 0.04 //Defining Constants could have been done with static constants too
# define u_0 40

using namespace std;

void FTCS(double *,int,double);
void velocityout(int, double *);
void gridout(double *);

int main()
{
  //Main Constants
    double nu=0.000217;  //Viscosity of fluid

 //Problem Dependent Constants
    double dt=0.00232; //Time Step
    double dy=0.001; //Space Step
    double d=((nu*dt))/(pow(dy,2)); //Diffusion Number
    cout<<"The Value of Diffusion Number is:"<<d<<endl;  //You will find that diffusion number is greater than 0.5 and hence solution is unstable.
    float error_time=1.08;  //Error to be checked at this time
    double eta1=h/(2*sqrt(nu*error_time));
    cout<<eta1<<endl;


    double y[41];    //Array Declaration for Grid points
    double u[41];    //Array Declaration for Velocity

    int i;

    y[0]=0;                //Defining Grid point Array
    for(i=1;i<41;i++)
    {
        y[i]=y[i-1]+dy;
    }
    gridout(y); //Function output for grid point output in .txt file.


    for(i=1;i<40;i++)     //Initializing the Velocity array
    {
        u[i]=0;
    }
    u[40]=0; //These are constant Values
    u[0]=u_0; //This is a constant


    int time_steps=541; //Total Time Steps

    FTCS(u,time_steps,d);

    for(i=0;i<41;i++)
    {
        cout<<u[i]<<endl;
    }

    return 0;
}


//Function For Solving the difference equation using FTCS Method

void FTCS(double *pt,int t,double alpha)
{
    int i;
    double *pj;
    for(t=1;t<542;t++)
    {
        double velocity1[41];
        pj=velocity1;
        for(i=0;i<41;i++)
        {
           *(pj+i)=*(pt+i);
        }
        for(i=1;i<40;i++)
        {
            *(pt+i)=*(pj+i)+(alpha*(*(pj+(i+1))-(2*(*(pj+i)))+(*(pj+(i-1)))));
        }
        if(t==91||t==181||t==271||t==361||t==451||t==541)
            velocityout(t,pt);
    }
}


//Function for giving a Text Printout of the Grid Points

void gridout(double *pt)
{
    int i;
    ofstream griddata("gridpoints.txt");
    if(griddata.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        griddata<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}

}


//Function for giving a Text printout of a Velocity Field at Different times
void velocityout(int time, double *pt)
{
    int i;

    if(time==91)
    {
    ofstream fout("0.18sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }

    else if(time==181)
    {
    ofstream fout("0.36sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }

    else if(time==271)
    {
    ofstream fout("0.54sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }

    else if(time==361)
    {
    ofstream fout("0.72sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }

    else if(time==451)
    {
    ofstream fout("0.90sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }

    else if(time==541)
    {
    ofstream fout("1.08sec.txt");
    if(fout.is_open())
    {
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;
    for(i=0;i<41;i++)
    {
        fout<<*(pt+i)<<endl;
    }
    }
    else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
    }
}

