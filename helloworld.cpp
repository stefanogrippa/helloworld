#include <iostream>
#include <vector>
#include <string>

// CMakeProject1.cpp: definisce il punto di ingresso dell'applicazione.
//
#define _USE_MATH_DEFINES

#include <thread>

// https://stackoverflow.com/questions/66734911/monte-carlo-integration-for-area-under-curve
#include <math.h>
//#include <cmath>
#include <ctime>
//#include <xmemory>
using namespace std;
#ifdef USA_MERTON
// https://www.quantstart.com/articles/Jump-Diffusion-Models-for-European-Options-Pricing-in-C/
// Standard normal probability density function
double norm_pdf(const double x) {
    return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}

// An approximation to the cumulative distribution function
// for the standard normal distribution
// Note: This is a recursive function
double norm_cdf(const double x) {
    double k = 1.0 / (1.0 + 0.2316419 * x);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return (1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
    }
    else {
        return 1.0 - norm_cdf(-x);
    }
}

// This calculates d_j, for j in {1,2}. This term appears in the closed
// form solution for the European call or put price
double d_j(const int j, const double S, const double K, const double r, const double v, const double T) {
    return (log(S / K) + (r + (pow(-1, j - 1)) * 0.5 * v * v) * T) / (v * (pow(T, 0.5)));
}

// Calculate the European vanilla call price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double bs_call_price(const double S, const double K, const double r,
    const double sigma, const double T) {
    return S * norm_cdf(d_j(1, S, K, r, sigma, T)) - K * exp(-r * T) *
        norm_cdf(d_j(2, S, K, r, sigma, T));
}

// Calculate the Merton jump-diffusion price based on 
// a finite sum approximation to the infinite series
// solution, making use of the BS call price.
double bs_jd_call_price(const double S, const double K, const double r,
    const double sigma, const double T, const int N, const double m,
    const double lambda, const double nu) {
    double price = 0.0;  // Stores the final call price
    double factorial = 1.0;

    // Pre-calculate as much as possible
    double lambda_p = lambda * m;
    double lambda_p_T = lambda_p * T;

    // Calculate the finite sum over N terms
    for (int n = 0; n < N; n++) {
        double sigma_n = sqrt(sigma * sigma + n * nu * nu / T);
        double r_n = r - lambda * (m - 1) + n * log(m) / T;

        // Calculate n!
        if (n == 0) {
            factorial *= 1;
        }
        else {
            factorial *= n;
        }

        // Refine the jump price over the loop
        price += ((exp(-lambda_p_T) * pow(lambda_p_T, n)) / factorial) *
            bs_call_price(S, K, r_n, sigma_n, T);
    }
    std::cout << "Call Price under JD:      " << price << "\n" << std::endl;
    return price;
}

/////////////////////////////////////////////////////////////////////////////////////////

double pow(double x, int n) {
    double r = 1.0;
    for (int i = 0; i < n; i++) {
        r = r * x;
    }
    return x;
}

bool belongTo_quadrifolium(double x, double y) {
    return pow(pow(x, 2) + pow(y, 2), 3) - 4 * (pow(x, 2) * pow(y, 2)) < 0;
}



double montecarlo(double accuracy) {
    unsigned int n = 0, d = 0;
    double x, y;
    do
    {
        x = 1.0 * rand() / (RAND_MAX - 1);
        y = 1.0 * rand() / (RAND_MAX - 1);
        if (belongTo_quadrifolium(x, y)) {
            d++;
            n++;
        }
        else {
            n++;
        }
    //} while (d < 3100 || (1.0 / n > accuracy));
    } while (d < 3100 && (1.0 / n > accuracy));
    double mypi = 4.0 * d / n;
    cout << "pi from thread =" << mypi << endl;
    return mypi;
}
double pi_(double accuracy)
{
    int n = 0, d = 0;
    double x, y, latest_pi = 0;
    double Origin_dist = 0;
    do
    {
        x = 0;
        y = 0;
        x = rand() % 100;
        y = rand() % 100;
        Origin_dist = sqrt(x * x + y * y);
        if (Origin_dist < 100.0)
        {
            d++;
            n++;
        }
        else
        {
            n++;
        }
        latest_pi = 4.0 * (d + 1.0) / (n + 1.0);
    } while ((d < 3100) || (4.0 / (n + 1.0) < accuracy));
    cout << "pi from thread =" << latest_pi << endl;
    return latest_pi;
}
#endif
/////////////////////////////////////////////////////////////////////////////////////////



// The function we want to execute on the new thread.
static void mytaskAction(std::string msg)
{
    cout << "task says: " << msg << "\n";
}
class Bar
{
public:
    void operator()(int a)
    {
        std::cout << a << '\n';
    }
};

int main()
{
    
    double const accuracy = 0;
    srand((int)time(0));
    //cout << "Enter the accuracy: \n";
    //cin >> accuracy;
    //cout << pi_(accuracy) << endl;

    // Constructs the new thread and runs it. Does not block execution.
    std::thread t1(mytaskAction, "Hello from task 1");
    std::thread t2(mytaskAction, "Hello from task 2");
    Bar bar;
    std::thread t3(bar, 10); // Pass 10 to functor object
    //thread tPI(montecarlo, accuracy);
    // Do other things...

    // First we create the parameter list
    double S = 100.0;     // Option price
    double K = 100.0;     // Strike price
    double r = 0.01;      // Risk-free rate (5%)
    double v = 0.2;       // Volatility of the underlying (20%)
    double T = 1.0;       // One year until expiry
    int N = 50;           // Terms in the finite sum approximation
    double m = 1.083287;  // Scale factor for J
    double lambda = 1.0;  // Intensity of jumps
    double nu = 0.4;      // Stdev of lognormal jump process

    // Then we calculate the call jump-diffusion value
    #ifdef USA_MERTON
    double call_jd = bs_jd_call_price(S, K, r, v, T, N, m, lambda, nu);
    std::cout << "Call Price under JD:      " << call_jd << "\n" << std::endl;
    #endif
    // Makes the main thread wait for the new thread to finish execution, therefore blocks its own execution.
    //thread tMerton(bs_jd_call_price, S, K, r, v, T, N, m, lambda, nu);
    t1.join();
    t2.join();
    t3.join();
    //tMerton.join();
    //tPI.join();
    

    cout << "Hello CMake." << endl;
    return 0;
}