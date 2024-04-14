#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <functional>

template <typename FunType>
std::pair<double, bool> RombergIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 50){
    if(eps<0 || nmin<0 || nmax<0 || nmax<nmin) 
        throw std::domain_error("Bad parameter");
    int okrenuteGranice=1;
    if(a>b){ 
        okrenuteGranice*=-1;
        std::swap(a,b);
    }

    int N=2; 
    double h=(b-a)/N, s=(f(a)+f(b))/2;
    double Iold=s;
    std::vector<double> I(30);
    for(int i=0; i<30 && N<nmax; i++){
        for(int j=1; j<=N/2; j++) s+=f(a+(2*j-1)*h);
        I[i]=h*s;
        double p=4;
        for(int k=i-1; k>=0; k--){
            I[k]=(p*I[k+1]-I[k])/(p-1);
            p*=4;
        }
        if(N>=nmin && std::fabs(I[0]-Iold)<=eps) return{I[0]*okrenuteGranice, true};

        Iold=I[0];
        h/=2;
        N*=2;
    }
    //throw std::domain_error("Tacnost nije postignuta");
    return {Iold*okrenuteGranice, false};
}

template <typename FunType>
std::pair<double, bool> TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 20, double range = 3.5){
    if(eps<0 || nmin<0 || nmax<0 || nmax<nmin || range<0) 
        throw std::domain_error("Bad parameter");
    int okrenuteGranice=1;
    if(a>b){ 
        okrenuteGranice*=-1;
        std::swap(a,b);
    }

    int N=2; 
    double h=2*range/N, p=(b+a)/2, q=(b-a)/2, s=0;
    double Iold=s, I=0, pi=4*std::atan(1);
    long double t=0, u=0, v=0;
    while(N<nmax){
        for(int i=1; i<=N/2; i++){
            t=-range+(2*i-1)*h;
            u=pi*std::sinh(t)/2;
            v=f(p+q*tanh(u));
            if(std::isfinite(v)) 
                s+=q*pi*std::cosh(t)*v/(2*std::cosh(u)*std::cosh(u));
        }
        I=h*s;
        if(N>=nmin && std::fabs(I-Iold)<eps) 
            return{I*okrenuteGranice, true};
        Iold=I;
        N*=2;
        h/=2;
    }
    //throw std::domain_error("Tacnost nije postignuta");
    return {I*okrenuteGranice, false};
}

template <typename FunType>
std::pair<double, bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2, double f3, int maxdepth){
    if(!std::isfinite(f1)) f1 = 0;
    if(!std::isfinite(f2)) f2 = 0;
    if(!std::isfinite(f3)) f3 = 0;
    
    double c = (a + b) / 2;
    long double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;
    long double f4 = f((a + c) / 2);
    
    if(!std::isfinite(f4)) 
        f4 = 0;
    long double f5 = f((c + b) / 2);
    if(!std::isfinite(f5)) 
        f5 = 0;
    long double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;

    if(std::abs(I1 - I2) < eps)
        return {I2, true};
    if(maxdepth <= 0)
        return {I2, false};

    auto rez1 = AdaptiveAux(f, a, c, eps, f1, f3, f4, maxdepth - 1);
    auto rez2 = AdaptiveAux(f, c, b, eps, f3, f2, f5, maxdepth - 1);

    return {rez1.first + rez2.first, rez1.second && rez2.second};
}

template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1){
    if(eps < std::numeric_limits<double>::epsilon() || nmin <= 0 || maxdepth <= 0)
        throw std::domain_error("Bad parameter");

    if(a - b > std::numeric_limits<double>::epsilon()){
        auto rez = AdaptiveIntegration(f, b, a, eps, maxdepth, nmin);
        return {-rez.first, rez.second};
    }
    
    long double s = 0;
    long double h = (b - a) / nmin;
    bool rez = true;

    for(int i = 1; i <= nmin; i++){
        auto aux = AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f(a + h / 2), maxdepth);
        s += aux.first;
        if(!aux.second) rez = false;
        a += h;
    }

    return {s, rez};
}

int main(){
    const double pi=4*std::atan(1);
    //izuzeci Romberg
    try{
        std::pair<double, double> romberg= RombergIntegration([](double x){return std::sin(x);}, 0, pi, -1e-10);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Romberg: "<<e.what()<<std::endl;
    }
    //nmax
    try{
        std::pair<double, double> romberg= RombergIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, -1000000);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Romberg: "<<e.what()<<std::endl;
    }
    //nmin
    try{
        std::pair<double, double> romberg= RombergIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, 1000000, -50);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Romberg: "<<e.what()<<std::endl;
    }
    //nmin>nmax
    try{
        std::pair<double, double> romberg= RombergIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, 50, 1000000);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Romberg: "<<e.what()<<std::endl;
    }
    //funkcija sin(x)
    try{
        std::pair<double, double> romberg= RombergIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, 1000000, 50);
        std::cout<<"Integracija Rombergom zavrsena.";
        if(romberg.second) std::cout<<" Vrijednost je: "<<romberg.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Romberg: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl;
    //izuzeci tanh-sinh
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, -1e-10);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //nmax
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, -1000000);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //nmin
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 1000000, -50);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinhn: "<<e.what()<<std::endl;
    }
    //nmin>nmax
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 50, 1000000);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //R<0
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 1000000, 50, -3.5);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //funkcija sin(x)/x
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 1000000, 50);
        std::cout<<"Integracija TanhSinh zavrsena.";
        if(tanhSinh.second) std::cout<<" Vrijednost je: "<<tanhSinh.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //funckija 1/sqrt(x)
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return 1/std::sqrt(x);}, 0, 1);
        std::cout<<"Integracija TanhSinh zavrsena.";
        if(tanhSinh.second) std::cout<<" Vrijednost je: "<<tanhSinh.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }
    //funkcija sin(x)
    try{
        std::pair<double, double> tanhSinh= TanhSinhIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, 1000000, 50);
        std::cout<<"Integracija TanhSinh zavrsena.";
        if(tanhSinh.second) std::cout<<" Vrijednost je: "<<tanhSinh.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska TanhSinh: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl;
    //izuzeci Adaptivna eps<0
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return std::sin(x)/x;}, 0, 1, -1e-10);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    //maxdepth
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, -30);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    //nmin
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 30, -1);
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    //funkcija sin(x)/x
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return std::sin(x)/x;}, 0, 1, 1e-8, 1000000, 50);
        std::cout<<"Integracija Adaptivna zavrsena.";
        if(adaptive.second) std::cout<<" Vrijednost je: "<<adaptive.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    //funckija 1/sqrt(x)
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return 1/std::sqrt(x);}, 0, 1);
        std::cout<<"Integracija Adaptivna zavrsena.";
        if(adaptive.second) std::cout<<" Vrijednost je: "<<adaptive.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    //funkcija sin(x)
    try{
        std::pair<double, double> adaptive= AdaptiveIntegration([](double x){return std::sin(x);}, 0, pi, 1e-8, 1000000, 50);
        std::cout<<"Integracija Adaptivna zavrsena.";
        if(adaptive.second) std::cout<<" Vrijednost je: "<<adaptive.first<<", a tacnost je postignuta\n";
        else std::cout<<"Tacnost nije postignuta\n";
    }
    catch(std::domain_error &e){
        std::cout<<"Greska Adaptivna: "<<e.what()<<std::endl;
    }
    return 0;
}