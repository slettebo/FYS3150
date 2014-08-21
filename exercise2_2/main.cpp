#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

int main()
{
    int n = 1;
    int i;
    unsigned int N;
    double difference = 0;

    //double N;
    //cout << "Sum(1/n) from n=1 to n=N & from n=N to n=1. Enter N:";
    //cin >> N;
    //N = int (N);

    float sum_up = 0;
    float sum_down = 0;

    // preparing file for writing
    ofstream myfile;
    myfile.open ("difference.dat");

    for (i=1; i<=10; i++){
        N = pow(10,i);
        while (n < N+1){
            sum_up += 1.0/n;
            sum_down += 1.0/(N-n+1);
            n++;
        }
        difference = abs((sum_up - sum_down)/sum_down);
        myfile << i << "\t" << difference <<"\n";
        cout << "N=10^" << i << ". Difference=" << difference << "\n";
    }
    //cout << setprecision(16) << "Sum up:" << sum_up << "\nSum down:" << sum_down << endl;
    myfile.close();
    return 0;
}
