#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <functional>

class ChebyshevApproximation{
    int red;
    double Xmin, Xmax;
    std::vector<double> koeficijenti;
    template <typename FunType>
    void izracunajKoef(FunType f);
   
    ChebyshevApproximation(const std::vector<double> koef, double xmin, double xmax)
    : Xmin(xmin), Xmax(xmax), red(koef.size()-1) {
        koeficijenti=koef;
    }

 public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n);
    void set_m(int m);
    void trunc(double eps);
    double operator()(double x) const;
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

//implementacija ChebyshevApproximation
template <typename FunType>
ChebyshevApproximation::ChebyshevApproximation(FunType f, double xmin, double xmax, int n){
    if(n<1 || xmin>xmax|| std::fabs(xmin-xmax)<std::numeric_limits<double>::epsilon())
        throw std::domain_error("Bad parameters");
    ChebyshevApproximation::red=n;
    Xmin=xmin; Xmax=xmax;
    izracunajKoef(f);
}

template <typename FunType>
void ChebyshevApproximation::izracunajKoef(FunType f){
    std::vector<double> v, w;
    double pi=4*std::atan(1);
    for(int i=0; i<=red+1; i++) w.push_back(std::cos(pi*i/(2*red+2)));
    for(int i=0; i<=red/2; i++) v.push_back(f((Xmax+Xmin+(Xmax-Xmin)*w[2*i+1])/2));
    for(int i=(red/2)+1; i<=red; i++) v.push_back(f((Xmax+Xmin+(Xmax-Xmin)*w[2*red+1-2*i])/2));
    
    for(int k=0; k<=red; k++){
        long double s=0;
        for(int i=0; i<=red; i++){
            int p=(k*(2*i+1))%(4*red+4);
            if(p>2*red+2) p=4*red+4-p;
            if(p>red+1) s-=v[i]*w[2*red+2-p];
            else s+=v[i]*w[p];
        }
        koeficijenti.push_back(2*s/(red+1));
    }
}


void ChebyshevApproximation::set_m(int m){
    if(m<=1 || m>red) throw std::domain_error("Bad order");
    red=m;
}
   
void ChebyshevApproximation::trunc(double eps){
    if(eps<0) throw std::domain_error("Bad tolerance");

    while (red > 1 && std::abs(koeficijenti[red]) < eps) red--;

    if (red == 1 && std::abs(koeficijenti[0]) < eps)
        throw std::domain_error("Bad tolerance");
}

double ChebyshevApproximation::operator()(double x) const{
    if(x<Xmin || x>Xmax) throw std::domain_error("Bad argument");
    double t=(2*x-Xmin-Xmax)/(Xmax-Xmin);
    double p=1, q=t, s=koeficijenti[0]/2+koeficijenti[1]*t;
    for(int k=2; k<=red; k++){
        double r=2*t*q-p;
        s+=koeficijenti[k]*r;
        p=q;
        q=r;
    }
    return s;
}
    

double ChebyshevApproximation::derivative(double x) const{
    if(x<Xmin || x>Xmax) throw std::domain_error("Bad argument");
    double t=(2*x-Xmin-Xmax)/(Xmax-Xmin);
    double p=1.0, q=4*t, s=koeficijenti[1]+q*koeficijenti[2];
    for(int k=3; k<=red; k++){
        double r=k*(2*t*q/(k-1)-p/(k-2));
        s+=koeficijenti[k]*r;
        p=q;
        q=r;
    }
    return 2*s/(Xmax-Xmin);
    //vjv treba formula sa str 231 ali sta je s u njoj??
}
    
ChebyshevApproximation ChebyshevApproximation::derivative() const{
    std::vector<double> koef(red);
    long double mi=4/(Xmax-Xmin);
    koef[red-1]=mi*red*koeficijenti[red];
    koef[red-2]=mi*(red-1)*koeficijenti[red-1];
    for(int k=red-3; k>=0; k--) koef[k]=koef[k+2]+mi*(k+1)*koeficijenti[k+1];

    return ChebyshevApproximation(koef, Xmin, Xmax);
}
    
ChebyshevApproximation ChebyshevApproximation::antiderivative() const{
    std::vector<double> koef(red+2);
    for(int k=1; k<=red-1; k++) koef[k]=(Xmax-Xmin)/(4*k)*(koeficijenti[k-1]-koeficijenti[k+1]);
    koef[0]=0; 
    koef[red]=(Xmax-Xmin)/(4*red)*koeficijenti[red-1]; 
    koef[red+1]=(Xmax-Xmin)/(4*red+4)*koeficijenti[red];
    return ChebyshevApproximation (koef, Xmin, Xmax);
}

double ChebyshevApproximation::integrate(double a, double b) const{
    //if (a > b) return -integrate(b, a);

    if (a < Xmin || b > Xmax)
        throw std::domain_error("Bad interval");

    ChebyshevApproximation antideriv = this->antiderivative();
    return antideriv(b) - antideriv(a);
}

double ChebyshevApproximation::integrate() const{
    const double pi = 4*std::atan(1);
    ChebyshevApproximation integral=this->antiderivative();
    /*long double rez=(Xmax-Xmin)/2, suma=0;
    for(int i=1; i<=(red+1)/2; i++) suma+=(koeficijenti[2*i-2]-koeficijenti[2*i])/(2*i-1);
    return rez*suma;*/
    return integral(Xmax)-integral(Xmin);
}


int main(){
    double pi=4*std::atan(1);
    try{
        //konstruktor baca izuzetak xmax<Xmin
        try{
            ChebyshevApproximation t1([](double x){return std::sin(x);},pi,0,7);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska konstruktor: "<<e.what()<<std::endl;
        }
        //konstruktor baca izuzetak n<1
        try{
            ChebyshevApproximation t1([](double x){return std::sin(x);},0,pi,0);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska konstruktor: "<<e.what()<<std::endl;
        }

        //kontruktor je ok
        ChebyshevApproximation t1([](double x){return std::sin(x);},0,pi,7);
        std::cout<<"Objekat je uspjesno kreiran\n";

        //izuzetak iz set m
        try{
            t1.set_m(1);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska set_m: "<<e.what()<<std::endl;
        }
        try{
            t1.set_m(9);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska set_m: "<<e.what()<<std::endl;
        }

        t1.set_m(5);
        std::cout<<"Red je promijenjen sa 7 na 5\n";

        //izuzetak trunc m<1
        try{
            t1.trunc(2.0);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska trunc: "<<e.what()<<std::endl;
        }

        //izuzetak trunc eps<0
        ChebyshevApproximation t([](double x){return std::sin(x);},0,pi,7);
        try{
            t.trunc(-1);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska trunc: "<<e.what()<<std::endl;
        }

        t.trunc(0.2);
        std::cout<<"Red je promijenjen metodom trunc\n";

        //izuzetak iz operatora x=2pi
        try{
            t1(2*pi);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska u metodi operator: "<<e.what()<<std::endl;
        }

        std::cout<<"Aproksimirana vrijednost funkcije sin(x) u tacki 0.0 je: "<<t1(0.0)<<std::endl;

        //izuzetak iz derivative x=2pi
        try{
            t.derivative()(2*pi);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska u metodi derivative: "<<e.what()<<std::endl;
        }
        
        std::cout<<"Aproksimirana vrijednost prvog izvoda funkcije sin(x) u tacki 0.0 je: "<<t1.derivative(0.0)<<std::endl;

        t1=ChebyshevApproximation([](double x){return std::sin(x);},0,pi,7);
        ChebyshevApproximation t2=t1.derivative();
        std::cout<<"Uspjesno kreiran objekat koji predstavlja aproksimaciju prvog izvoda funkcije sin(x)\n";
        std::cout<<"Vrijednost prvog izvoda u 0.0 je: "<<t2(0.0)<<std::endl;

        ChebyshevApproximation t3=t1.antiderivative();
        std::cout<<"Uspjesno kreiran objekat koji predstavlja aproksimaciju antiderivacije funkcije sin(x)\n";
        std::cout<<"Njena vrijednost u 0.0 je: "<<t3(0.0)<<std::endl;
        
        //izuzetak integrate - granice izvan [xmin, xmax]
        try{
            t2.integrate(0,2*pi);
        }
        catch(std::domain_error &e){
        std::cout<<"Greska u metodi integrate sa paramerima: "<<e.what()<<std::endl;
        }
        //granice se mogu zamijeniti
        try{
            std::cout<<"Ako zamijenimo granice, integral cos(x) na [pi,0]: "<<t2.integrate(pi,0)<<std::endl;
        }
        catch(std::domain_error &e){
        std::cout<<"Greska u metodi integrate sa paramerima: "<<e.what()<<std::endl;
        }
        //integral t2 koja je cos(x)
        std::cout<<"Integral cos(x) na [0, pi/2]: "<<t2.integrate(0,pi/2)<<std::endl;
        //integral t3 koja je -cos(x)
        std::cout<<"Integral -cos(x) na [0, pi/2]: "<<t3.integrate(0,pi/2)<<std::endl;
        //integral na cijelom intervalu
        std::cout<<"Integral cos(x) na [0, pi]: "<<t2.integrate()<<std::endl;
        std::cout<<"Integral -cos(x) na [0, pi]: "<<t3.integrate()<<std::endl;

        //Ismarov -cos(pi/4)
        t1=ChebyshevApproximation ([](double x){return std::sin(x);},0,pi,7);
        std::cout<<"-cos(pi/4)="<<t1.antiderivative()(pi/4);
        auto temp=t1.antiderivative();
        std::cout<<"\n-cos(pi/4)="<<temp(pi/4);

        //Ismarov polinom 3x^4-2x^3+5x^2-7x+1
        temp=ChebyshevApproximation([](double x){return 3*x*x*x*x-2*x*x*x+5*x*x-7*x+1;},-1,1,10);
        std::cout<<"\nVrijednost integrala funkcije 3x^4-2x^3+5x^2-7x+1 u 0.8 je "<<temp.antiderivative()(0.8);
        
    }
    catch(std::domain_error &e){
        std::cout<<e.what()<<std::endl;
    }
    return 0;
}