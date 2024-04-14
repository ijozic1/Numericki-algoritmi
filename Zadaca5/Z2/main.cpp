#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>

template <typename FunType>
void BracketMin(FunType f, double x0, double &a, double &b, double &c, double hinit, double hMax, double lambda){
    a=x0-hinit;
    c=x0;
    b=x0+hinit;

    double fa=f(a), fb=f(b), fc=f(c);

    if(fb<fc || (hinit*=-1, fa<fc)){
        while(std::fabs(hinit)<hMax){
            a=c-hinit;
            fa=f(a);

            b=c+hinit;
            fb=f(b);

            while(!std::isfinite(fa)){
                hinit/=2*(1+lambda);
                a=c-hinit;
                fa=f(a);
            }

            while(!std::isfinite(fb)){
                hinit/=2*(1+lambda);
                b=c+hinit;
                fb=f(b);
            }

            if(fc<fa && fc<fb) break;

            hinit*=lambda;
            c=b;
            fc=fb;
        }
        if(std::fabs(hinit)>hMax || std::fabs(std::fabs(hinit)-hMax)<std::numeric_limits<double>::epsilon())
            throw std::logic_error("Minimum has not found");
    }
    if(a>b) std::swap(a,b);
}

template <typename FunType>
double FindMinimum(FunType f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    if(eps<0 || hinit<0 || hmax<0 || lambda<0)
        throw std::domain_error("Invalid parameters");
    
    double a=0.0, b=0.0, c=0.0, d=0.0;
    try{
        BracketMin(f, x0, a, b, c, hinit, hmax, lambda);
    }
    catch(std::logic_error &e){
        throw e;
    }

    const double zlatniPresjek=(1+std::sqrt(5))/2;

    if(std::fabs(c-a)<std::fabs(b-c)) d=b-(b-c)/zlatniPresjek;
    else{
        d=c;
        c=a+(c-a)/zlatniPresjek;
    }

    double u=f(c), v=f(d);
    
    while (std::fabs(b-a)>eps) {
        if(u<v){
            b=d; 
            d=c;
            c=a+(c-a)/zlatniPresjek;
            v=u;
            u=f(c);
        }
        else{
            a=c;
            c=d;
            d=b-(b-d)/zlatniPresjek;
            u=v;
            v=f(d);
        }
    }
    return (a + b) / 2.0;
}

int main(){
    try{
        //eps<0
        double r1=FindMinimum([](double x){return std::sin(x);}, 0.1, -0.5);
        std::cout<<"Rjesenje: "<<r1<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska FindMinimum: "<<e.what()<<std::endl;
    }
    try{
        //hinit<0
        double r1=FindMinimum([](double x){return std::sin(x);}, 0.1, 1e-8, -0.5);
        std::cout<<"Rjesenje: "<<r1<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska FindMinimum: "<<e.what()<<std::endl;
    }
    try{
        //hmax<0
        double r1=FindMinimum([](double x){return std::sin(x);}, 0.1, 1e-8, 1e-5, -0.5);
        std::cout<<"Rjesenje: "<<r1<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska FindMinimum: "<<e.what()<<std::endl;
    }
    try{
        //lambda<0
        double r1=FindMinimum([](double x){return std::sin(x);}, 0.1, 1e-8, 1e-5, 1e10, -0.5);
        std::cout<<"Rjesenje: "<<r1<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska FindMinimum: "<<e.what()<<std::endl;
    }

    try {
        double rez = FindMinimum([](double x) { return x * x; }, 0.0);
        std::cout << "Minimum pronadjen u x = " << rez << std::endl;
    } 
    catch (std::exception &e) {
        std::cout <<"Greska FindMinimum: "<< e.what() << std::endl;
    }
    try {
        //Minimum nije pronađen (nije ogradjen)
        double rez = FindMinimum([](double x) { return x * x; }, 1.0, 1e-8, 1e-5, 1e10, 1.4);
        std::cout << "Minimum pronadjen u x = " << rez << std::endl;
    } 
    catch (std::logic_error &e) {
        std::cout << "Greska FindMinimum: " << e.what() << std::endl;
    }

    try {
        // Nema minimuma
        double rez = FindMinimum([](double x) { return x*x*x;}, 1.0);
        std::cout << "Minimum pronadjen u x = " << rez << std::endl;
    } 
    catch (std::logic_error &e) {
        std::cout << "Greska FindMinimum: " << e.what() << std::endl;
    }
    return 0;
}