#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <ctime>

template <typename FunType>
bool OneWayBracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    a=x0;
    double fa=f(a);

    while(std::fabs(hinit)<hmax){
        b=a+hinit;
        double fb=f(b);

        while(!std::isfinite(fb)){
            hinit/=2*(1+lambda);
            if(std::fabs(hinit)<std::fabs(a)*std::numeric_limits<double>::epsilon() || 
                std::fabs(std::fabs(hinit)-std::fabs(a)*std::numeric_limits<double>::epsilon())<std::numeric_limits<double>::epsilon())
                    return false;
            b=a+hinit;
            fb=f(b);
        }

        if(fa*fb<0 || std::fabs(fa*fb)<std::numeric_limits<double>::epsilon())
            return true;
        
        hinit*=lambda;
        a=b;
        fa=fb;
    }

    return false;
}

template <typename FunType>
bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    if(hinit<0 || hmax<0 || lambda <0) 
        throw std::domain_error("Invalid parameters");
    if(!OneWayBracketRoot(f, x0, a, b, hinit, hmax, lambda)){
        if(!OneWayBracketRoot(f, x0, a, b, -hinit, hmax, lambda)) return false;
        std::swap(a,b);
    }
    return true;
}

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename FunType>
double UnmodifiedRF(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    double fa=f(a), fb=f(b), c=a, cOld=b;

    while(--maxiter && std::fabs(c-cOld)>eps){
        cOld=c;
        c=(a*fb-b*fa)/(fb-fa);

        double fc=f(c);
        if(std::fabs(fc)<std::numeric_limits<double>::epsilon()) return c;
        if(fa*fc<0){
            b=c;
            fb=fc;
        }
        else{
            a=c;
            fa=fc;
        }
    }

    if(maxiter==0)
        throw std::logic_error("Given accuracy has not achieved");
    return c;
}

template <typename FunType>
double IllinoisRF(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    double fa=f(a), fb=f(b), c=a, cOld=b;

    while(--maxiter && std::fabs(c-cOld)>eps){
        cOld=c;
        c=(a*fb-b*fa)/(fb-fa);

        double fc=f(c);
        if(std::fabs(fc)<std::numeric_limits<double>::epsilon()) return c;
        if(fa*fc<0){
            b=a;
            fb=fa;
        }
        else fb/=2;

        a=c;
        fa=fc;
    }
    
    if(maxiter==0)
        throw std::logic_error("Given accuracy has not achieved");
    return c;
}

template <typename FunType>
double SlavicRF(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    try{
        return UnmodifiedRF([f](double x){
            double y=f(x);
            return y/(1+std::fabs(y));
        }, a, b, eps, maxiter);
    }
    catch(std::logic_error &e){
        throw e;
    }
}

template <typename FunType>
double IllinoisSlavicRF(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    try{
        return IllinoisRF([f](double x){
            double y=f(x);
            return y/(1+std::fabs(y));
        }, a, b, eps, maxiter);
    }
    catch(std::logic_error &e){
        throw e;
    }
}

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100){
    if(f(a)*f(b)>0)
        throw std::range_error("Root must be bracketed");
    if(eps<0 || maxiter<0) 
        throw std::domain_error("Invalid parameters");

    if(a>b) std::swap(a,b);

    switch (mode) {
        case Unmodified:
            try{
                return UnmodifiedRF(f, a, b, eps, maxiter);
            }
            catch(std::logic_error &e){
                throw e;
            }
        case Illinois:
            try{
                return IllinoisRF(f, a, b, eps, maxiter);
            }
            catch(std::logic_error &e){
                throw e;
            }
        case Slavic:
            try{
                return SlavicRF(f, a, b, eps, maxiter);
            }
            catch(std::logic_error &e){
                throw e;
            }
        case IllinoisSlavic:
            try{
                return IllinoisSlavicRF(f, a, b, eps, maxiter);
            }
            catch(std::logic_error &e){
                throw e;
            }
    }
    throw std::logic_error("Error...");
}

int znak(double n){
    if(n<0) return -1;
    if(n>0) return 1;
    return 0;
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    if(f(a)*f(b)>0)
        throw std::range_error("Root must be bracketed");
    if(eps<0 || maxiter<0) 
        throw std::domain_error("Invalid parameters");
    
    if(a>b) std::swap(a,b);
    
    double fa = f(a),fb = f(b);
    while(--maxiter && std::fabs(b-a)>eps){
        double c=(a+b)/2, fc=f(c);
        if(std::fabs(fc)<std::numeric_limits<double>::epsilon()) return c;
        
        double d=c+fc*(c-a)*znak((fa-fb))/std::sqrt(fc*fc-fa*fb);
        double fd=f(d);
        if(std::fabs(fd)<std::numeric_limits<double>::epsilon()) return d;
        
        if(fc*fd<0 || std::fabs(fc*fd)<std::numeric_limits<double>::epsilon()){
            a=c;
            fa=fc;
            b=d;
            fb=fd;
        }
        else if(fa*fd<0 || std::fabs(fa*fd)<std::numeric_limits<double>::epsilon()){
            b=d;
            fb=fd;
        }
        else{
            a=d;
            fa=fd;
        }
    }
    
    if(maxiter==0)
        throw std::logic_error("Given accuracy has not achieved");
    return (a + b) / 2.0;
}

template <typename FunType1, typename FunType2>
double NewtonRaphsonSolve(FunType1 f, FunType2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100){
    if(eps<0 || maxiter<0) 
        throw std::domain_error("Invalid parameters");
    if(damping <0 || damping>1 || std::fabs(damping - 1)<std::numeric_limits<double>::epsilon())
        throw std::domain_error("Invalid parameters");
    
    double deltaX=std::numeric_limits<double>::infinity();
    double v=f(x0), d=fprim(x0);
    bool noDamping = (std::fabs(damping)<std::numeric_limits<double>::epsilon());

    while(--maxiter && std::fabs(deltaX)>eps){
        if(std::fabs(v)<eps || std::fabs(v-eps)<std::numeric_limits<double>::epsilon())
            return x0;
        
        deltaX=v/d;
        double w=v;
        v=f(x0-deltaX);
        d=fprim(x0-deltaX);

        if(noDamping){
            if(!std::isfinite(v) || std::fabs(d)<std::numeric_limits<double>::epsilon())
                throw std::logic_error("Convergence has not achieved");
        }
        else{
            while(std::fabs(v)>std::fabs(w) || !std::isfinite(v) || std::fabs(d)<std::numeric_limits<double>::epsilon()){
                deltaX*=damping;
                v=f(x0-deltaX);
                d=fprim(x0-deltaX);
            }
        }
        x0-=deltaX;
    }
    if(maxiter==0)
        throw std::logic_error("Convergence has not achieved");
    return x0;
}

std::complex<double> operator*(double n, std::complex<double>c){
    return{n*c.real(), n*c.imag()};
}

double Random(){
    //vraca random broj [-10, 10]
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    double random_broj = static_cast<double>(std::rand()) / RAND_MAX * 20.0 - 10.0;

    return random_broj;
}

template<typename Type>
std::pair<std::complex<double>,bool> LaguerreRoot(const std::vector<Type>& coeffs, int order, std::complex<double>x0, double eps, int maxiter){
    std::complex<double> deltaX=std::numeric_limits<double>::infinity();
    int iter=1;
    while(iter<maxiter && std::fabs(deltaX) > eps){
        std::complex<double>f=coeffs[order], d=0, s=0;

        for(int i=order-1; i>=0; i--){
            s=s*x0 + 2*d;
            d=d*x0 + f;
            f=f*x0 + coeffs[i];
        }

        if(std::fabs(f)<eps || std::fabs(std::fabs(f)-eps)<std::numeric_limits<double>::epsilon())
            return{x0, true};
        
        std::complex<double> r=std::sqrt((order-1)*((order-1)*d*d-order*f*s));
        if(std::fabs(d+r)>std::fabs(d-r)) deltaX=order*f/(d+r);
        else deltaX=order*f/(d-r);

        x0-=deltaX;
        iter++;
    }
    if((std::fabs(deltaX)<eps || std::fabs(std::fabs(deltaX)-eps)<std::numeric_limits<double>::epsilon()))
        return {x0,true};
    
    return {x0, false};
}

std::vector<std::complex<double>> PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps < 0 || maxiters <= 0 || maxtrials <= 0)
        throw std::domain_error("Invalid parameters");
    
    std::vector<std::complex<double>> koef=coefficients, nule;
    int vel=coefficients.size();
    std::complex<double> tacnaNula(0.0);

    for(int i=vel-1; i>0; i--){
        int t=1;
        bool c=false;

        std::pair<std::complex<double>,bool> p(0.0, false);
        for(int l=1; !c && l<maxtrials; l++)
            p=LaguerreRoot(koef, i, {Random(), Random()}, eps, maxiters);
        
        if(!p.second) 
            throw std::logic_error("Convergence has not achieved");
        

        while(!c && t<maxtrials){
            tacnaNula.real(Random());
            tacnaNula.imag(Random());
            auto par=LaguerreRoot(coefficients, i, tacnaNula, eps, maxiters);
            tacnaNula=par.first;
            c=par.second;

            //poliranje nula 
            par=LaguerreRoot(koef, vel-1, tacnaNula, eps, maxiters);
            tacnaNula=par.first;
            c=par.second;

            t++;
        }
        if(!c) 
            throw std::logic_error("Convergence has not achieved");
        if(std::fabs(tacnaNula.imag())<eps)
            tacnaNula.imag(0.0);
        
        nule.push_back(tacnaNula);
        std::complex<double> v=coefficients[i];

        for(int j=i-1; j>=0; j--){
            std::complex<double>w=coefficients[j];
            coefficients[j]=v;
            v=w + tacnaNula*v;
        }
    }

    std::sort(nule.begin(), nule.end(), [](const std::complex<double> &c1, const std::complex<double> &c2){
        return c1.real()<c2.real() || (std::fabs(c1.real()-c2.real())<std::numeric_limits<double>::epsilon() && c1.imag()<c2.imag());
    });
    return nule;
}

std::vector<std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps < 0 || maxiters <= 0 || maxtrials <= 0)
        throw std::domain_error("Invalid parameters");
    
    std::vector<double> koef=coefficients;
    std::vector<std::complex<double>>nule;
    int vel=coefficients.size();
    std::complex<double> tacnaNula(0.0);

    int i=vel-1;

    while(i>0){
        int t=1;
        bool c=false;

        std::pair<std::complex<double>,bool> p(0.0, false);
        for(int l=1; !c && l<maxtrials; l++)//
            p=LaguerreRoot(koef, i, {Random(), Random()}, eps, maxiters);
        
        if(!p.second) 
            throw std::logic_error("Convergence has not achieved");

        while(!c && t<maxtrials){
            tacnaNula.real(Random());
            tacnaNula.imag(Random());
            auto par=LaguerreRoot(coefficients, i, tacnaNula, eps, maxiters);
            tacnaNula=par.first;
            c=par.second;

            //poliranje nula 
            par=LaguerreRoot(koef, vel-1, tacnaNula, eps, maxiters);
            tacnaNula=par.first;
            c=par.second;

            t++;
        }
        if(!c)
            throw std::logic_error("Convergence has not achieved");
        if(std::fabs(tacnaNula.imag())<eps){
            tacnaNula.imag(0.0);
            nule.push_back(tacnaNula);
            double v=coefficients[i];

            for(int j=i-1; j>=0; j--){
                double w=coefficients[j];
                coefficients[j]=v;
                v=w + tacnaNula.real()*v;
            }
            i--;
        }
        else{
            nule.push_back(std::complex<double>(tacnaNula.real(),-tacnaNula.imag()));
            nule.push_back(tacnaNula);
            double alpha=2*tacnaNula.real(),
                   beta=std::fabs(tacnaNula)*std::fabs(tacnaNula);
            double u=coefficients[i], v=coefficients[i-1]+alpha*u;

            for(int j=i-2; j>=0; j--){
                double w=coefficients[j];
                coefficients[j]=u;
                u=v;
                v=w + alpha*v - beta*coefficients[j];
            }
            i-=2;
        }
    }

    std::sort(nule.begin(), nule.end(), [](const std::complex<double> &c1, const std::complex<double> &c2){
        return c1.real()<c2.real() || (std::fabs(c1.real()-c2.real())<std::numeric_limits<double>::epsilon() && c1.imag()<c2.imag());
    });
    return nule;
}

int main(){
    double pi=4*std::atan(1);
    std::cout<<"BracketRoot\n";
    //BracketRoot
    try{
        //hinit<0
        double rjesenjeA=-1., rjesenjeB=-1.;
        bool izuzetak1=BracketRoot([](double x){return std::sin(x);}, 0.1, rjesenjeA, rjesenjeB, -0.5);
        if(izuzetak1) std::cout<<"Rjesenje: a="<<rjesenjeA<<", b="<<rjesenjeB<<std::endl;
        else std::cout<<"Pretraga nije uspjela"<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska BracketRoot: "<<e.what()<<std::endl;
    }
    try{
        //hmax<0
        double rjesenjeA=-1., rjesenjeB=-1.;
        bool izuzetak2=BracketRoot([](double x){return std::sin(x);}, 0.1, rjesenjeA, rjesenjeB, 1e-5, -0.5);
        if(izuzetak2) std::cout<<"Rjesenje: a="<<rjesenjeA<<", b="<<rjesenjeB<<std::endl;
        else std::cout<<"Pretraga nije uspjela"<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska BracketRoot: "<<e.what()<<std::endl;
    }
    try{
        //lambda<0
        double rjesenjeA=-1., rjesenjeB=-1.;
        bool izuzetak3=BracketRoot([](double x){return std::sin(x);}, 0.1, rjesenjeA, rjesenjeB, 1e-5, 1e10, -0.5);
        if(izuzetak3) std::cout<<"Rjesenje: a="<<rjesenjeA<<", b="<<rjesenjeB<<std::endl;
        else std::cout<<"Pretraga nije uspjela"<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska BracketRoot: "<<e.what()<<std::endl;
    }
    try{
        //ok
        double a1=0, b1=0;
        bool rez1 = BracketRoot([](double x) { return x * x - 4; }, 0, a1, b1);
        std::cout << "BracketRoot rezultat (x^2-4): " << rez1 << ", a: " << a1 << ", b: " << b1 << std::endl;

        double a2=0, b2=0;
        bool rez2 = BracketRoot([](double x) { return std::cos(x);}, 0, a2, b2);
        std::cout << "BracketRoot rezultat (cos(x)): " << rez2 << ", a: " << a2 << ", b: " << b2 << std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska BracketRoot: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl<<"RegulaFalsiSolve\n";
    //RegulaFalsiSolve
    try{
        //f(a) i f(b) istog znaka - oba pozitivna
        double nula=RegulaFalsiSolve([](double x){return std::sin(x);}, 0.1, pi-0.1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::range_error &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }
    try{
        //f(a) i f(b) istog znaka - oba negativna
        double nula=RegulaFalsiSolve([](double x){return std::sin(x);}, pi+0.1, 2*pi-0.1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::range_error &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }
    try{
        //eps<0
        double nula=RegulaFalsiSolve([](double x){return std::sin(x);}, 0.1, pi+0.1, Slavic, -0.5);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }
    try{
        //maxiter<0
        double nula=RegulaFalsiSolve([](double x){return std::sin(x);}, 0.1, pi+0.1, Slavic, 1e-10, -1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }
    try{
        //ok
        double nula = RegulaFalsiSolve([](double x) { return x * x - 4; }, 1, 3);
        std::cout << "Rjesenje (x^2-4): nula=" << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }
    try{
        //ok
        double nula = RegulaFalsiSolve([](double x) { return std::cos(x); }, 0, 3);
        std::cout << "Rjesenje (cos(x)): nula=" << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska RegulaFalsiSolve: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl<<"RiddersSolve\n";
    //RiddersSolve
    try{
        //f(a) i f(b) istog znaka - oba pozitivna
        double nula=RiddersSolve([](double x){return std::sin(x);}, 0.1, pi-0.1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::range_error &e){
        std::cout<<"Greska RiddersSolve: "<<e.what()<<std::endl;
    }
    try{
        //f(a) i f(b) istog znaka - oba negativna
        double nula=RiddersSolve([](double x){return std::sin(x);}, pi+0.1, 2*pi-0.1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::range_error &e){
        std::cout<<"Greska RiddersSolve: "<<e.what()<<std::endl;
    }
    try{
        //eps<0
        double nula=RiddersSolve([](double x){return std::sin(x);}, 0.1, pi+0.1, -0.5);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska RiddersSolve: "<<e.what()<<std::endl;
    }
    try{
        //maxiter<0
        double nula=RiddersSolve([](double x){return std::sin(x);}, 0.1, pi+0.1, 1e-10, -1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska RiddersSolve: "<<e.what()<<std::endl;
    } 
    try{
        //ok
        double nula = RiddersSolve([](double x) { return x * x - 4; }, 1, 3);
        std::cout << "Rjesenje (x^2-4): nula=" << nula << std::endl;

        nula = RiddersSolve([](double x) { return std::cos(x); }, 0, 3);
        std::cout << "Rjesenje (cos(x)): nula=" << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska RiddersSolve: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl<<"NewtonRaphsonSolve\n";
    //NewtonRaphsonSolve
    try{
        //eps<0
        double nula=NewtonRaphsonSolve([](double x){return std::sin(x);}, [](double x){return std::cos(x);}, 0.1, -0.5);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    }
    try{
        //maxiter<0
        double nula=NewtonRaphsonSolve([](double x){return std::sin(x);}, [](double x){return std::cos(x);}, 0.1, 1e-10, 0, -1);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    } 
    try{
        //damping<0
        double nula=NewtonRaphsonSolve([](double x){return std::sin(x);}, [](double x){return std::cos(x);}, 0.1, 1e-10, -0.5, 100);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    }
    try{
        //damping>1
        double nula=NewtonRaphsonSolve([](double x){return std::sin(x);}, [](double x){return std::cos(x);}, 0.1, 1e-10, 1.5, 100);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    } 
    try{
        //damping=1
        double nula=NewtonRaphsonSolve([](double x){return std::sin(x);}, [](double x){return std::cos(x);}, 0.1, 1e-10, 1, 100);
        std::cout<<"Rjesenje: nula="<<nula<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    } 
    try{
        // ok
        double nula = NewtonRaphsonSolve([](double x) { return x * x - 4; }, [](double x) { return 2 * x; }, 3);
        std::cout << "Rjesenje (x^2-4): nula=" << nula << std::endl;

        nula= NewtonRaphsonSolve([](double x) { return std::cos(x); }, [](double x) { return -std::sin(x); }, 1);
        std::cout << "Rjesenje (cos(x)): nula=" << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    }
    try {
        //Newton ne konvergira za x<0
        double nula = NewtonRaphsonSolve([](double x) { return std::sqrt(x); }, [](double x) { return 0.5 / std::sqrt(x); }, -1);
        std::cout << "Rjesenje (sqrt(x)): nula=" << nula << std::endl;
    } 
    catch(std::logic_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    }
    try {
        //Newton ne konvergira za ovaj broj iteracija
        double nula = NewtonRaphsonSolve([](double x) { return x * x - 2; }, [](double x) { return 2 * x; }, 0, 1e-10, 10);
        std::cout << "Rjesenje (x^2-2): nula=" << nula << std::endl;
    } 
    catch(std::logic_error &e){
        std::cout<<"Greska NewtonRaphsonSolve: "<<e.what()<<std::endl;
    }

    std::cout<<std::endl<<"PolyRoots\n";
    //PolyRoots
    try{
        //eps<0
        std::vector<std::complex<double>> a=PolyRoots(std::vector<std::complex<double>>{1.2, 0.1, 2.3}, -0.5, 100, 10);
        std::cout<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska PolyRoots: "<<e.what()<<std::endl;
    }
    try{
        //maxiter<0
        std::vector<std::complex<double>> a=PolyRoots(std::vector<double>{1.2, 0.1, 2.3}, 1e-10, -1, 10);
        std::cout<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska PolyRoots: "<<e.what()<<std::endl;
    } 
    try{
        //maxtrials<0
        std::vector<std::complex<double>> a=PolyRoots(std::vector<std::complex<double>>{1.2, 0.1, 2.3}, 1e-10, 100, -1);
        std::cout<<std::endl;
    }
    catch(std::domain_error &e){
        std::cout<<"Greska PolyRoots: "<<e.what()<<std::endl;
    }
    try{
        //vektor double
        std::vector<double> koef = {1, 0, -4};
        auto nule = PolyRoots(koef);
        for (const auto &nula : nule)
            std::cout << "PolyRoots nule: " << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska PolyRoots: "<<e.what()<<std::endl;
    }
    try{
        //vektor complex
        std::vector<std::complex<double>> koef = {1, 0, -4}; 
        auto nule = PolyRoots(koef);
        for (const auto &nula : nule) 
            std::cout << "PolyRoots nule: " << nula << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska PolyRoots: "<<e.what()<<std::endl;
    }
    
    std::cout<<std::endl<<std::endl;
    return 0;
}