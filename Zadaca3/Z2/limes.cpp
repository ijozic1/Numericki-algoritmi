#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <functional>

template<typename FunType>
std::pair<double, bool> Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20){
    if (eps <std::numeric_limits<double>::epsilon() || nmax < 3 || nmax > 30)
        throw std::domain_error("Invalid parameters");

    bool inf=false;
    if(std::isinf(x0) && x0<0){
        inf=true;
        x0=0;
        x0=-x0;
    }
    else if(std::isinf(x0)){
        inf=true;
        x0=0;
    }

    if (std::fabs(h)<std::numeric_limits<double>::epsilon()){
       1.-std::fabs(x0)>1e-5 ? h=0.001 : h=0.001*std::fabs(x0);
    }
    
    double yOld = std::numeric_limits<double>::infinity();
    std::vector<double> y(nmax,0);
    for (int i = 0; i < nmax; i++) {
        if(inf) y[i]=f(1./x0+h);
        else y[i]=f(x0+h);
        double p=2;
        for(int k=i-1; k>=0; k--){
            y[k]=(p*y[k+1]-y[k])/(p-1);
            p=2*p;
        }
        if(std::fabs(y[0]-yOld)<eps) return std::make_pair(y[0],true);
        yOld=y[0];
        h/=2;
    }
    return std::make_pair(y[0],false);
}

int main(){
    try {
        // Testiranje limesa funkcije 1/x kad x teži ka 0
        auto rez = Limit([](double x) {return 1/x;}, 0);
        std::cout << "\nLimes 1/x kad x tezi 0: " << rez.first;
        if (rez.second)
            std::cout << " (postoji)\n";
        else 
            std::cout << " (ne postoji)\n";
        
        // Testiranje limesa funkcije 1/x kad x teži ka beskonacnosti
        auto rez1 = Limit([](double x) {return 1/x;}, std::numeric_limits<double>::infinity());
        std::cout << "Limes 1/x kad x tezi beskonacnosti: " << rez.first;
        if (rez1.second)
            std::cout << " (postoji)\n";
        else 
            std::cout << " (ne postoji)\n";
        
        // Testiranje limesa funkcije sin(x)/x kad x teži ka 0
        auto rez2 = Limit([](double x) {return (std::sin(x)/x);}, 0);
        std::cout << "Limes sin(x)/x kad x tezi 0: " << rez2.first;
        if (rez2.second)
            std::cout << " (postoji)\n";
        else 
            std::cout << " (ne postoji)\n";

        // Testiranje limesa funkcije e^(-x) kad x teži ka beskonačnosti
        auto rez3 = Limit([](double x) {return std::exp(-x);}, std::numeric_limits<double>::infinity());
        std::cout << "Limes e^(-x) kad x tezi ka beskonacnosti: " << rez3.first;
        if (rez3.second)
            std::cout << " (postoji)\n";
        else 
            std::cout << " (ne postoji)\n";
    } 
    catch (const std::exception &e) {
        std::cout<< "Greska: " << e.what() << std::endl;
    }
    return 0;
}