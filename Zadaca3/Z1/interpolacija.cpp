#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <functional>

class AbstractInterpolator{
    mutable int kesiraniIndeks; 
 protected:
    std::vector<std::pair<double, double>> tacke;
    int Locate(double x) const;

 public:
    AbstractInterpolator(const std::vector<std::pair<double, double>> &data);
    virtual double operator()(double x) const = 0;
};


class LinearInterpolator : public AbstractInterpolator{
  public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){}
    double operator()(double x) const override;
};

class PolynomialInterpolator : public AbstractInterpolator{
    std::vector<double> NewtonoviKoef;
    std::vector<double>y;

    void IzracunajKoeficijente();

 public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){
        IzracunajKoeficijente();
    }
    double operator()(double x) const override;
    void AddPoint(const std::pair<double, double> &p);
    std::vector<double> GetCoefficients() const;
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator{
    int stepen;
  public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data, int order):
        AbstractInterpolator(data){
            if(order<1 || order>=tacke.size()) 
                throw std::domain_error("Invalid order");
            stepen=order;
    }
    double operator()(double x) const override;
};

class SplineInterpolator : public AbstractInterpolator{
    std::vector<double> koeficijenti;
    std::vector<double> alpha;

    void buildSpline();

  public:
    SplineInterpolator(const std::vector<std::pair<double, double>> &data):
        AbstractInterpolator(data){ buildSpline(); }
    
    double operator()(double x) const override;
};

class BarycentricInterpolator : public AbstractInterpolator{
    std::vector<double>tezinskiFaktori;
    int red;

    void IzracunajTezinskeFaktore(int red);

  public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order): AbstractInterpolator(data){
        if (order < 0 || order >= tacke.size()) 
            throw std::domain_error("Invalid order");
        red=order;

        IzracunajTezinskeFaktore(red);
    }
    double operator()(double x) const override;
    std::vector<double> GetWeights() const{
        return tezinskiFaktori;
    }
};

class TrigonometricInterpolator : public AbstractInterpolator{ 
  public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){
        int vel=tacke.size();
        if (std::fabs(tacke[0].second - tacke[vel-1].second)>std::numeric_limits<double>::epsilon()){
            tacke.clear();
            throw std::domain_error("Function is not periodic.");
        }
    }

    double operator()(double x) const override;
};

//AbstractInterpolator implementacija
int AbstractInterpolator::Locate(double x) const{
    //prije prvog
    if (x <= tacke[0].first) {
        kesiraniIndeks = 0;
        return 0;
    }
    //nakon posljednjeg
    if (x > tacke[tacke.size()-1].first) {
        kesiraniIndeks = tacke.size();
        return kesiraniIndeks;
    }
    //u trenutnom
    if(kesiraniIndeks!=-1 && kesiraniIndeks>0 && x>tacke[kesiraniIndeks-1].first && x<=tacke[kesiraniIndeks].first)
        return kesiraniIndeks;
    //u prethodnom
    if(kesiraniIndeks!=-1 && kesiraniIndeks>1 && x>tacke[kesiraniIndeks-2].first && x<=tacke[kesiraniIndeks-1].first) 
        return --kesiraniIndeks;
    //u sljedecem
    if(kesiraniIndeks!=-1 && kesiraniIndeks<tacke.size()-1 && x>tacke[kesiraniIndeks].first && x<=tacke[kesiraniIndeks+1].first) 
        return ++kesiraniIndeks;

    /*auto it = std::lower_bound(tacke.begin(), tacke.end(), std::make_pair(x, 0.0),
                                [](const std::pair<double, double>& a, const std::pair<double, double>& b) {return a.first < b.first;});*/
    auto it=std::lower_bound(tacke.begin(), tacke.end(),x,
                            [](const std::pair<double, double> &p, double vrijednost){
                                return p.first<vrijednost;});

    kesiraniIndeks = std::distance(tacke.begin(), it);
    return kesiraniIndeks;
}

AbstractInterpolator::AbstractInterpolator(const std::vector<std::pair<double, double>> &data){
    /*for (int i = 0; i < data.size(); i++) {
        for(int j=i+1; j<data.size(); j++)
            if (std::fabs(data[i].first -data[j].first)<std::numeric_limits<double>::epsilon())
                throw std::domain_error("Invalid data set");
    }*/
    tacke=data;
    kesiraniIndeks=-1;
    std::sort(tacke.begin(), tacke.end(), [](const std::pair<double, double>a, const std::pair<double, double> &b){
        if(std::fabs(a.first-b.first)<std::numeric_limits<double>::epsilon()) 
            throw std::domain_error("Invalid data set");
        return a.first<b.first;});
}

//LinearInterpolator implementacija
double LinearInterpolator::operator()(double x) const{
    int index = Locate(x);
    double y;

    if (index == 0)
        y=((tacke[1].first-x)/(tacke[1].first-tacke[0].first))*tacke[0].second + ((x-tacke[0].first)/(tacke[1].first-tacke[0].first))*tacke[1].second;
    else if (index == tacke.size()){
        int vel=tacke.size();
        y=((tacke[vel-1].first-x)/(tacke[vel-1].first-tacke[vel-2].first))*tacke[vel-2].second + 
            ((x-tacke[vel-2].first)/(tacke[vel-1].first-tacke[vel-2].first))*tacke[vel-1].second;
    }
    else 
        y=((tacke[index].first-x)/(tacke[index].first-tacke[index-1].first))*tacke[index-1].second + 
            ((x-tacke[index-1].first)/(tacke[index].first-tacke[index-1].first))*tacke[index].second;

    return y; 
}

//PolynomialInterpolator implementacija
void PolynomialInterpolator::IzracunajKoeficijente(){
    int brTacaka=tacke.size();
    NewtonoviKoef.clear();

    for(int i=0; i<brTacaka; i++)
        NewtonoviKoef.push_back(tacke[i].second);
        
    y.push_back(NewtonoviKoef[brTacaka-1]);

    for(int j=1; j<brTacaka; j++){
        for(int i=brTacaka-1; i>=j; i--)
            NewtonoviKoef[i]=(NewtonoviKoef[i]-NewtonoviKoef[i-1])/(tacke[i].first - tacke[i-j].first);
        y.push_back(NewtonoviKoef[brTacaka-1]);
    }
}

double PolynomialInterpolator::operator()(double x) const{
    int brTacaka=tacke.size();
    double f=NewtonoviKoef[brTacaka-1];
    for(int i=brTacaka-1; i>0; i--)
        f=f*(x-tacke[i-1].first)+NewtonoviKoef[i-1];
    return f;
}

void PolynomialInterpolator::AddPoint(const std::pair<double, double> &p){
    int brTacaka=tacke.size();
    for (int i = 0; i < brTacaka; i++) {
        if (std::fabs(tacke[i].first - p.first)<std::numeric_limits<double>::epsilon())
                throw std::domain_error("Invalid point");
    }
    tacke.push_back(p);
    NewtonoviKoef.push_back(p.second);
    double koef=p.second;

    for(int i=1; i<=brTacaka; i++){
        double temp=y[i-1];
        y[i-1]=koef;
        koef=(koef-temp)/(tacke[brTacaka].first-tacke[brTacaka-i].first);
    }

    NewtonoviKoef[brTacaka]=koef;
    y.push_back(koef);
}

std::vector<double> PolynomialInterpolator::GetCoefficients() const{
    int brTacaka=tacke.size();
    std::vector<double> p(brTacaka, 0.0), w(brTacaka+1), v(brTacaka+1);

    w[0]=1;
    for(int i=1; i<=brTacaka; i++){
        w[i]=w[i-1];
        for(int j=i-1; j>0; j--)
            w[j]=w[j-1]-tacke[i-1].first*w[j];
        w[0]=-tacke[i-1].first*w[0];
    }
    for(int i=0; i<brTacaka; i++){
        double alpha=1.0;
        for(int j=0; j<brTacaka; j++)
            if(j!=i) alpha*=tacke[i].first-tacke[j].first;
        alpha=tacke[i].second/alpha;

        for(int j=0; j<=brTacaka; j++) v[j]=w[j];
        for(int j=brTacaka-1; j>=0; j--){
            v[j]+=tacke[i].first*v[j+1];
            p[j]+=alpha*v[j+1];
        }
    }
    return p;
}

//PiecewisePolynomialInterpolator implementacija
double PiecewisePolynomialInterpolator::operator()(double x) const {
    int index=Locate(x);
    int brTacaka=tacke.size();

    if (index == brTacaka || x <= tacke[0].first) index-=2;
    else if (index >0) index--;

    int pocetak=0, kraj=0;

    if(stepen%2==1 && index<(stepen-1)/2 || (stepen%2==0 && index<stepen/2)){
        pocetak=0; 
        kraj=stepen;
    }
    else if((stepen%2==1 && index+(stepen+1)/2>brTacaka-1) || (stepen%2==0 && index+stepen/2>brTacaka-1)){
        pocetak=brTacaka-stepen-1;
        kraj=brTacaka-1;
    }
    else if(stepen%2==1){
        pocetak=index-(stepen-1)/2;
        kraj=index+(stepen+1)/2;
    }
    else if(stepen%2==0){
        pocetak=index-stepen/2;
        kraj=index+stepen/2;
    }

    long double s=0;

    for(int i=pocetak; i<=kraj; i++){
        long double p=tacke[i].second;
        for(int j=pocetak; j<=kraj; j++){
            if(i!=j)
                p*=(x-tacke[j].first)/(tacke[i].first-tacke[j].first);
        }
        s+=p;
    }
    return s;
}

//SplineInterpolator implementacija
void SplineInterpolator::buildSpline() {
    int brIntervala = tacke.size();

    koeficijenti.resize(brIntervala, 0.0);
    alpha.resize(brIntervala-1, 0.0);
        

    for (int i = 1; i < brIntervala-1; i++) {
        alpha[i] = 2.0*(tacke[i+1].first-tacke[i-1].first);
        koeficijenti[i] = 3.0 * ((tacke[i + 1].second - tacke[i].second) / (tacke[i + 1].first - tacke[i].first) - 
                        (tacke[i].second - tacke[i-1].second)/(tacke[i].first - tacke[i-1].first));
    }

    for (int i = 1; i < brIntervala-2; i++) {
        long double mi=(tacke[i+1].first-tacke[i].first)/alpha[i];
        alpha[i+1] = alpha[i+1]-mi*(tacke[i + 1].first - tacke[i].first);
        koeficijenti[i+1]=koeficijenti[i+1]-mi*koeficijenti[i];
    }
    koeficijenti[brIntervala-2]/=alpha[brIntervala-2];

    for (int j = brIntervala - 3; j>0 ; j--) {
         koeficijenti[j] = (koeficijenti[j]-(tacke[j+1].first-tacke[j].first)*koeficijenti[j+1])/alpha[j];
    }
}

double SplineInterpolator::operator()(double x) const{
    int index = Locate(x);
    if(index==tacke.size()) index-=2;
    else if(index>0) index--;

    long double t, dX, s, q;
    t=x-tacke[index].first;
    dX=tacke[index+1].first-tacke[index].first;
    s=(koeficijenti[index+1]-koeficijenti[index])/(3*dX);
    q=(tacke[index+1].second-tacke[index].second)/dX-dX*(koeficijenti[index+1]+2*koeficijenti[index])/3;
        
    return tacke[index].second+t*(q+t*(koeficijenti[index]+t*s));
}

//BarycentricInterpolator implementacija
void BarycentricInterpolator::IzracunajTezinskeFaktore(int red){
    int brTacaka = tacke.size();
    tezinskiFaktori.resize(brTacaka, 0.0);

    for (int i = 0; i < brTacaka; i++) {
        tezinskiFaktori[i] = 0.0;
        double p=1;
        int max=std::max(0, i - red), min=std::min(brTacaka - red, i );
        for (int k = max; k <= min; k++) {
            p=1;
            for(int j=k; j<k+red; j++){
                if (i != j) {
                    tezinskiFaktori[i]/= (tacke[i].first - tacke[j].first);
                }
            }
            if(k%2==1) p=-p;
        }
        tezinskiFaktori[i] += p;
    }
}

double BarycentricInterpolator::operator()(double x) const{
    double p = 0.0, q = 0.0;

    for (int i = 0; i < tacke.size(); i++) {
        if(std::fabs(x-tacke[i].first) < std::numeric_limits<double>::epsilon())
            return tacke[i].second;
        double w = tezinskiFaktori[i] / (x - tacke[i].first);
        p += w * tacke[i].second;
        q += w;
    }
    return p / q;
}

//TrigonometricInterpolator implementacija
double TrigonometricInterpolator::operator()(double x) const{
    double rez = 0.0;
    int brTacaka = tacke.size();

    for (int i = 0; i < brTacaka; i++) {
        double term = tacke[i].second;

        for (int j = 1; j < brTacaka; j++) {
            if (j % 2 == 1) 
                term += 2.0 * tacke[j].second * std::cos(2.0 * M_PI * j * (x - tacke[i].first));
            else
                term += 2.0 * tacke[j].second * std::sin(2.0 * M_PI * j * (x - tacke[i].first));
        }
        rez += term;
    }
    return rez;
}


// Funkcija za generisanje podataka za spline interpolaciju
std::vector<std::pair<double, double>> generateSplineTestData() {
    std::vector<std::pair<double, double>> data;
    for (double x = 0.0; x <= 2.0 * M_PI; x += 0.2) {
        data.emplace_back(x, std::sin(x));
    }
    return data;
}

int main(){
    //provjera za klasu AbstractInterpolator
    /*try{
        std::vector<std::pair<double, double>> data1 = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        AbstractInterpolator li1(data1);
        std::cout<<"li1 kreiran\n";

        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        AbstractInterpolator interpolator(data);

        int index = interpolator.Locate(0.5);
        std::cout << "Index: " << index << std::endl;

        index = interpolator.Locate(4.0);
        std::cout << "Index: " << index << std::endl;

        index = interpolator.Locate(2.5);
        std::cout << "Index: " << index << std::endl;

        std::vector<std::pair<double, double>> data2 = {{1.0, 2.0}, {1.0, 4.0}, {2.0, 3.0}};
        AbstractInterpolator li2(data2);
        std::cout<<"li2 kreiran";
    }
    catch(std::domain_error &e){
        std::cout<<"Konstruktor bacio izuzetak: "<<e.what()<<std::endl;
    }*/

    //provjera za LinearInterpolator
    try{
        std::vector<std::pair<double, double>> data1 = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        LinearInterpolator li1(data1);
        std::cout<<"li1 kreiran\n";

        double y=li1(0.5);
        std::cout << "Vrijednost u 0.5: " << y << std::endl;

        y=li1(4.0);
        std::cout << "Vrijednost u 4.0: " << y << std::endl;

        y=li1(2.5);
        std::cout << "Vrijednost u 2.5: " << y << std::endl;
        
        //ovdje treba baciti gresku
        std::vector<std::pair<double, double>> data2 = {{1.0, 2.0}, {1.0, 4.0}, {2.0, 3.0}};
        LinearInterpolator li2(data2);
        std::cout<<"li2 kreiran";
    }
    catch(std::domain_error &e){
        std::cout<<"Konstruktor bacio izuzetak: "<<e.what()<<std::endl;
    }

    //provjera za polinomski interpolator
    try{
        //ovdje treba baciti gresku
        std::vector<std::pair<double, double>> data2 = {{1.0, 2.0}, {1.0, 4.0}, {2.0, 3.0}};
        LinearInterpolator li2(data2);
        std::cout<<"li2 kreiran";
    }
    catch(std::domain_error &e){
        std::cout<<"Konstruktor bacio izuzetak: "<<e.what()<<std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        PolynomialInterpolator interpolator(data);

        double f = interpolator(2.5);
        std::cout << "Interpolirana vrijednost: " << f << std::endl;

        interpolator.AddPoint({4.0, 5.0});

        std::vector<double> NewtonoviKoef = interpolator.GetCoefficients();
        std::cout << "Koeficijenti: ";
        for (double koef : NewtonoviKoef) {
            std::cout << koef << " ";
        }
        std::cout << std::endl;

        //addPoint treba baciti izuzetak
        try{
            interpolator.AddPoint({1.0, 5.0});
        }
        catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
        }
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }

    //provjera za PieceWisePolynomial
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        int red=0;
        PiecewisePolynomialInterpolator interpolator(data,red);
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        int red=3;
        PiecewisePolynomialInterpolator interpolator(data,red);
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        int red=2;
        PiecewisePolynomialInterpolator interpolator(data,red);
        std::cout<<"Interpolator kreiran\n";
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}, {4.0, 5.0}};
        int red = 2;

        PiecewisePolynomialInterpolator interpolator(data, red);

        double rezultat = interpolator(2.5);
        std::cout << "Interpolirana vrijednost: " << rezultat << std::endl;
    } 
    catch (const std::domain_error& e) {
        std::cerr << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        int red=2;

        PiecewisePolynomialInterpolator interpolator(data, red);

        // na početku opsega
        double rez1 = interpolator(0.5);
        std::cout << "Extrapolacija u x = 0.5: " << rez1 << std::endl;

        // na kraju opsega
        double rez2 = interpolator(6.0);
        std::cout << "Extrapolacija u x = 6.0: " << rez2 << std::endl;

        //unutar opsega
        double rez3 = interpolator(2.5);
        std::cout << "Interpolacija u x = 2.5: " << rez3 << std::endl;

    } 
    catch (const std::domain_error& e) {
        std::cerr << "Greska: " << e.what() << std::endl;
    }

    //provjera SplineInterpolatora
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        SplineInterpolator splineInterpolator(data);

        // Testiranje splajna u nekim tačkama
        for (double x = 0.0; x <= 6.0; x += 1.0) {
            double rez = splineInterpolator(x);
            std::cout << "Spline u x = " << x << ": " << rez << std::endl;
        }

    } 
    catch (const std::exception& e) {
        std::cout<< "Greska: " << e.what() << std::endl;
    }

    //provjera baricentricnog interpolatora
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        int red=-1;
        BarycentricInterpolator interpolator(data,red);
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        int red=4;
        BarycentricInterpolator interpolator(data,red);
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        int red = 1;
        BarycentricInterpolator barycentricInterpolator(data, red);

        for (double x = 0.0; x <= 6.0; x += 1.0) {
            double rez = barycentricInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        std::vector<double> weights = barycentricInterpolator.GetWeights();
        std::cout << "Tezinski faktori: ";
        for (size_t i = 0; i < weights.size(); i++) {
            std::cout << weights[i] << " ";
        }
        std::cout << std::endl;

    }  
    catch (const std::exception& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }

    //provjera trigonometrijskog interpolatora
    try {
        std::vector<std::pair<double, double>> data = {{1.0, 2.0}, {3.0, 4.0}, {2.0, 3.0}};
        TrigonometricInterpolator interpolator(data);
        //treba baciti izuzetak zbog neperiodicnosti
    } 
    catch (const std::domain_error& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }
    try {
        std::vector<std::pair<double, double>> data = {{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}};
        TrigonometricInterpolator trigInterpolator(data);

        for (double x = -1.0; x <= 3.0; x += 0.5) {
            double rez = trigInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

    } 
    catch (const std::exception& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }

    //provjera za instrukcije iz postavke
    try {
        // PolynomialInterpolator test
        std::cout << "PolynomialInterpolator Test:" << std::endl;
        PolynomialInterpolator polyInterpolator({{-1.0, 3.0}, {0.0, 1.0}, {1.0, 1.0}});
        for (double x = -1.0; x <= 1.0; x += 0.5) {
            double rez = polyInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        // LinearInterpolator test
        std::cout << "\nLinearInterpolator Test:" << std::endl;
        LinearInterpolator linearInterpolator({{-1.0, 3.0}, {0.0, 1.0}, {1.0, 1.0}});
        for (double x = -1.0; x <= 1.0; x += 0.5) {
            double rez = linearInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        // SplineInterpolator test
        std::cout << "\nSplineInterpolator Test:" << std::endl;
        std::vector<std::pair<double, double>> splineData = generateSplineTestData();
        SplineInterpolator splineInterpolator(splineData);
        for (double x = 0.0; x <= 2.0 * M_PI; x += 0.5) {
            double rez = splineInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        // PiecewisePolynomialInterpolator test
        std::cout << "\nPiecewisePolynomialInterpolator Test:" << std::endl;
        int red=2;
        PiecewisePolynomialInterpolator piecewisePolyInterpolator({{-1.0, 3.0}, {0.0, 1.0}, {1.0, 1.0}},red);
        for (double x = -1.0; x <= 1.0; x += 0.5) {
            double rez = piecewisePolyInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        // BarycentricInterpolator test
        std::cout << "\nBarycentricInterpolator Test:" << std::endl;
        BarycentricInterpolator barycentricInterpolator({{-1.0, 3.0}, {0.0, 1.0}, {1.0, 1.0}}, 2);
        for (double x = -1.0; x <= 1.0; x += 0.5) {
            double rez = barycentricInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }

        // TrigonometricInterpolator test
        std::cout << "\nTrigonometricInterpolator Test:" << std::endl;
        TrigonometricInterpolator trigInterpolator(generateSplineTestData());
        for (double x = 0.0; x <= 2.0 * M_PI; x += 0.5) {
            double rez = trigInterpolator(x);
            std::cout << "Interpolacija u x = " << x << ": " << rez << std::endl;
        }
    } 
    catch (const std::exception& e) {
        std::cout << "Greska: " << e.what() << std::endl;
    }

    return 0;
}