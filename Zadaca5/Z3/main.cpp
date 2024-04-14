#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

template <typename FunType>
double RK4Step(FunType f, double x, double y, double h) {
    long double k1 = f(x, y);
    long double k2 = f(x + h / 2, y + h * k1 / 2);
    long double k3 = f(x + h / 2, y + h * k2 / 2);
    long double k4 = f(x + h, y + h * k3);

    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

template <typename FunType>
std::vector<std::pair<double, double>> RK4Integrator(FunType f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false) {
    std::vector<std::pair<double, double>> rezultat;
    bool negativanH=(h<-std::numeric_limits<double>::epsilon());

    if(!adaptive){
        while((!negativanH && ((xmax+h/1000)-x0 > std::numeric_limits<double>::epsilon() || 
                std::fabs((xmax+h/1000)-x0) < std::numeric_limits<double>::epsilon())) ||
            (negativanH && (x0-(xmax+h/1000) > std::numeric_limits<double>::epsilon() ||
                std::fabs((xmax+h/1000)-x0) < std::numeric_limits<double>::epsilon()))){
            
            rezultat.push_back({x0,y0});
            y0=RK4Step(f, x0, y0, h);
            x0+=h;
        }
    }
    else{
        while((!negativanH && (xmax-x0 > std::numeric_limits<double>::epsilon() || 
                std::fabs(xmax-x0) < std::numeric_limits<double>::epsilon())) ||
            (negativanH && (x0-xmax > std::numeric_limits<double>::epsilon() ||
                std::fabs(xmax-x0) < std::numeric_limits<double>::epsilon()))){

            double u=RK4Step(f,x0,y0,h/2);
            double v=RK4Step(f,x0+h/2,u, h/2);
            double w=RK4Step(f,x0,y0,h);

            double delta=std::fabs(w-v)/std::fabs(h);

            if(delta<eps || std::fabs(delta-eps)<std::numeric_limits<double>::epsilon()){
                if((!negativanH && (x0+h)-xmax > std::numeric_limits<double>::epsilon()) || 
                        (negativanH && xmax-(x0+h) > std::numeric_limits<double>::epsilon()))
                            break;
                x0+=h;
                y0=v;
                rezultat.push_back({x0,y0});
            }
            
            long double t=0.9*std::sqrt(std::sqrt(eps/delta));
            if(t-5>std::numeric_limits<double>::epsilon()) t=5;
            h*=t;
        }

        h=xmax-x0;
        double u=RK4Step(f,x0,y0,h/2);
        double v=RK4Step(f,x0+h/2,u, h/2);
        double w=RK4Step(f,x0,y0,h);

        double delta=std::abs(w-v)/std::fabs(h);

        if(delta<eps || std::fabs(delta-eps)<std::numeric_limits<double>::epsilon())
            rezultat.push_back({xmax,v});
    }
    return rezultat;
}


template <typename FunType>
std::vector<std::pair<double, std::vector<double>>> RK4SystemIntegrator(FunType f, double x0, std::vector<double> y0, double xmax, double h) {
    std::vector<std::pair<double, std::vector<double>>> rezultat;
    bool negativanH=(h<-std::numeric_limits<double>::epsilon());

    if((!negativanH && x0-xmax > std::numeric_limits<double>::epsilon()) ||
            (negativanH && xmax - x0 > std::numeric_limits<double>::epsilon())){
        rezultat.push_back({x0,y0});
        return rezultat;
    }
    
    int vel=y0.size();

    while((!negativanH && ((xmax+h/1000)-x0 > std::numeric_limits<double>::epsilon() || 
                std::fabs((xmax+h/1000)-x0) < std::numeric_limits<double>::epsilon())) ||
            (negativanH && (x0-(xmax+h/1000) > std::numeric_limits<double>::epsilon() ||
                std::fabs((xmax+h/1000)-x0) < std::numeric_limits<double>::epsilon()))){
        
        rezultat.push_back({x0,y0});
        
        std::vector<double> k1=f(x0,y0),u;
        if(k1.size()!= y0.size())
            throw std::range_error("Incompatible formats");
        for(int k=0; k<vel; k++) u.push_back(y0[k] + h*k1[k]/2);
        x0+=h/2;
        
        std::vector<double> k2=f(x0,u),v;
        if(k2.size()!= y0.size())
            throw std::range_error("Incompatible formats");
        for(int k=0; k<vel; k++) v.push_back(y0[k] + h*k2[k]/2);

        std::vector<double> k3=f(x0,v);
        if(k3.size()!= y0.size())
            throw std::range_error("Incompatible formats");
        for(int k=0; k<vel; k++) u.push_back(y0[k] + h*k3[k]);
        x0+=h/2;
        
        std::vector<double> k4=f(x0,u);
        if(k4.size()!= y0.size())
            throw std::range_error("Incompatible formats");
        for(int k=0; k<vel; k++) y0[k]+=h*(k1[k] + 2*k2[k] + 2*k3[k]+ k4[k]) /6;

    }
    return rezultat;
}

int main() {
    try{
        //y' = 3x + 2y
        auto rjesenje1 = RK4Integrator([](double x, double y) { return 3 * x + 2 * y; }, 0, 1, 1, 0.1);

        std::cout << "Rjesenje za y' = 3x + 2y:\n";
        for (const auto& t : rjesenje1)
            std::cout << "(" << t.first << ", " << t.second << ")\n";
    }
    catch(std::exception &e){
        std::cout<<"Greska: "<<e.what()<<std::endl;
    }
    try{
        // sistem dif. jednacina
        std::vector<double> y0 = {1, 0};
        auto rjesenje2 = RK4SystemIntegrator([](double x, const std::vector<double>& y) {
            return std::vector<double>{2 * y[0] + 3 * y[1] + x, y[0] - 2 * y[1] + 1};}, 0, y0, 1, 0.1);

        std::cout << "\nRjesenje sistema:\n";
        for (const auto& t : rjesenje2) {
            std::cout << "(" << t.first << ", ";
            for (const auto& v : t.second)
                std::cout << v << " ";
            std::cout << ")\n";
        }
    }
    catch(std::exception &e){
        std::cout<<"Greska: "<<e.what()<<std::endl;
    }
    try {
        //y' = -2y, h<0
        auto rjesenje3 = RK4Integrator([](double x, double y) { return -2 * y; }, 2, 1, 0, -0.1);

        std::cout << "\nRjesenje za y' = -2y:\n";
        for (const auto& t : rjesenje3)
            std::cout << "(" << t.first << ", " << t.second << ")\n";
    } 
    catch (std::exception& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    
    return 0;
}