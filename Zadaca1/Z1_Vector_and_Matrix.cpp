#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <iomanip>

class Vector{
    std::vector<double> vektor;
  public:
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);
    int NElems() const;
    double &operator[](int i);
    double operator[](int i) const;
    double &operator()(int i);
    double operator()(int i) const;
    double Norm() const;
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const;
    void Print(char separator = '\n', double eps = -1) const;
    friend void PrintVector(const Vector &v, char separator , double eps);
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v);
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v);
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s);
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s);
};

class Matrix{
    std::vector<std::vector<double>> matrica;
  public:
    Matrix(int m, int n);
    Matrix(const Vector &v);
    Matrix(std::initializer_list<std::vector<double>> l);
    int NRows() const;
    int NCols() const;
    double *operator[](int i);
    const double *operator[](int i) const;
    double &operator()(int i, int j);
    double operator()(int i, int j) const;
    double Norm() const;
    friend double MatrixNorm(const Matrix &m);
    double GetEpsilon() const;
    void Print(int width = 10, double eps = -1) const;
    friend void PrintMatrix(const Matrix &m, int width, double eps);
    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    Matrix &operator +=(const Matrix &m);
    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    Matrix &operator -=(const Matrix &m);
    friend Matrix operator *(double s, const Matrix &m);
    friend Matrix operator *(const Matrix &m, double s);
    Matrix &operator *=(double s);
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix &operator *=(const Matrix &m);
    friend Vector operator *(const Matrix &m, const Vector &v);
    friend Matrix Transpose(const Matrix &m);
    void Transpose();
};

//implementacija metoda iz klase Vector
Vector::Vector(int n){
    if(n<=0) throw std::range_error("Bad dimension");
    vektor.resize(n,0);
}

Vector::Vector(std::initializer_list<double> l){
    if(l.size()==0) throw std::range_error("Bad dimension");
    for(double v : l) vektor.push_back(v);
}

int Vector::NElems()const{
    return vektor.size();
}

double & Vector::operator[](int i){
    return vektor[i]; //mogu li u ovoj i sljedecoj koristiti at, da program ne bi krahirao
}

double Vector::operator[](int i) const{
    return vektor[i];
}

double & Vector::operator()(int i){
    if(i>vektor.size() || i<1) throw std::range_error("Invalid index");
    return vektor[i-1];
}

double Vector::operator()(int i) const{
    if(i>vektor.size() || i<1) throw std::range_error("Invalid index");
    return vektor[i-1];
}

double Vector::Norm() const{
    //ovdje sam radi problema dodavanja malih brojeva na velike mogla uvesti pomocni vektor 
    //auto v1=this->vektor;
    //std::sort(v1.begin(), v1.end(),std::less<double>());
    //i onda bi u rasponskoj vrtila v1 umjesto v
    double sumaKomponenti=0;
    int vel=NElems();
    for(double v:this->vektor) sumaKomponenti+=(v*v);
    //for(int i=0; i<vel; i++) sumaKomponenti+=(this->vektor.at(i)*this->vektor.at(i));
    return std::sqrt(sumaKomponenti);
}

double VectorNorm(const Vector &v){
    return v.Norm();
}

double Vector::GetEpsilon() const{
    return 10*this->Norm()*std::numeric_limits<double>::epsilon();
}

void Vector::Print(char separator, double eps) const{
    if (eps<0) eps=this->GetEpsilon();
    int vel=this->NElems();
    for(int i=0; i<vel; i++){
        if(std::fabs(vektor[i])<eps) std::cout<<'0';
        else std::cout<<vektor[i];
        
        if(i!=vel-1) std::cout<<separator;
        else if(i==vel-1 && separator=='\n') std::cout<<separator;
    }
}

void PrintVector(const Vector &v, char separator = '\n', double eps= -1){
    v.Print(separator, eps);
}

Vector operator +(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3=v1;
    return v3+=v2;
    /*Vector v3(v1.NElems());
    int vel=v3.NElems();
    //for(int i=0; i<vel; i++) v3.vektor.at(i)=v1.vektor.at(i)+v2.vektor.at(i);
    for(int i=0; i<vel; i++) v3.vektor[i]=v1.vektor[i]+v2.vektor[i];
    return v3;*/
}
    
Vector & Vector::operator +=(const Vector &v){
    if(this->NElems()!=v.NElems()) throw std::domain_error("Incompatible formats");
    int vel=this->NElems();
    for(int i=0; i<vel; i++) this->vektor[i]+=v.vektor[i];
    return *this;
}

Vector operator -(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    /*Vector v3(v1.NElems());
    int vel=v3.NElems();
    for(int i=0; i<vel; i++) v3.vektor[i]=v1.vektor[i]-v2.vektor[i];
    return v3;*/
    Vector v3=v1;
    return v3-=v2;
}

Vector & Vector::operator -=(const Vector &v){
    if(this->NElems()!=v.NElems()) throw std::domain_error("Incompatible formats");
    int vel=this->NElems();
    for(int i=0; i<vel; i++) this->vektor[i]-=v.vektor[i];
    return *this;
}

Vector operator *(double s, const Vector &v){
    /*Vector novi(v.NElems());
    int vel=v.NElems();
    for(int i=0; i<vel; i++) novi.vektor[i]=v.vektor[i]*s;
    return novi;*/
    Vector v1=v;
    return v1*=s;
}

Vector operator *(const Vector &v, double s){
    return s*v;
}

Vector & Vector::operator *=(double s){
    int vel=this->NElems();
    //for(double & v : this->vektor) v*=s;
    for(int i=0; i<vel; i++) this->vektor[i]*=s;
    return *this;
}

double operator *(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    int vel=v1.NElems();
    double skalarniP=0;
    for(int i=0; i<vel; i++) skalarniP+=(v1.vektor[i]*v2.vektor[i]);
    return skalarniP;
}

Vector operator /(const Vector &v, double s){
    if(std::fabs(s-0)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    //moze if(s==0) sasvim je legitimno, isto i u sljedecoj
    /*int vel=v.NElems();
    Vector novi(vel);
    for(int i=0; i<vel; i++) novi.vektor[i]=v.vektor[i]/s;
    return novi;*/
    Vector v1=v;
    return v1/=s;
}

Vector & Vector::operator /=(double s){
    if(std::fabs(s-0)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    int vel=this->NElems();
    for(int i=0; i<vel; i++) this->vektor[i]/=s;
    return *this;
}

//implementacija metoda iz klase Matrix
Matrix::Matrix(int m, int n){
    if(m<=0 || n<=0) throw std::range_error("Bad dimension");
    matrica.resize(m);
    for(int i=0; i<m; i++) matrica[i].resize(n);
}

Matrix::Matrix(const Vector &v){
    matrica.resize(v.NElems());
    for(int i=0; i<matrica.size(); i++) matrica.at(i).push_back(v[i]);
}

Matrix::Matrix(std::initializer_list<std::vector<double>> l){
    if(l.size()==0 || l.begin()->size()==0) throw std::range_error("Bad dimension");
    int velReda=l.begin()->size();
    for(const auto &v : l){
        if(v.size()==0) throw std::range_error("Bad dimension");
        if(v.size()!=velReda) throw std::logic_error("Bad matrix");
    }
    for(const auto &v:l) matrica.emplace_back(v);
}

int Matrix::NRows() const{
    return matrica.size();
}

int Matrix::NCols() const{
    return matrica[0].size();
}

double * Matrix::operator[](int i){
    return & matrica[i][0];
    //return & matrica.at(i).at(0);
}

const double * Matrix:: operator[](int i) const{
    return & matrica[i][0];
    //return & matrica.at(i).at(0);
}

double & Matrix::operator()(int i, int j){
    if(i<1 || j<1 || i>this->matrica.size() || j>this->matrica[0].size())
        throw std::range_error("Invalid index");
    return this->matrica.at(i-1).at(j-1);
}

double Matrix::operator()(int i, int j) const{
    if(i<1 || j<1 || i>this->matrica.size() || j>this->matrica[0].size())
        throw std::range_error("Invalid index");
    return this->matrica.at(i-1).at(j-1);
}

double Matrix::Norm() const{
    double sumaKomponenti=0;
    int brRed=this->NRows(), brKol=this->NCols();
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            sumaKomponenti+=(this->matrica.at(i).at(j)*this->matrica.at(i).at(j));
    return std::sqrt(sumaKomponenti);
}

double MatrixNorm(const Matrix &m){
    return m.Norm();
}

double Matrix::GetEpsilon() const{
    return 10*this->Norm()*std::numeric_limits<double>::epsilon();
}

void Matrix::Print(int width, double eps) const{
    if(eps<0) eps=this->GetEpsilon();
    auto w=std::setw(width);
    int brRed=this->NRows(), brKol=this->NCols();
    for(int i=0; i<brRed; i++){
        for(int j=0; j<brKol; j++) {
            if(std::fabs(matrica.at(i).at(j))<eps) std::cout<<w<<'0';
            else std::cout<<w<<this->matrica.at(i).at(j);
        }
        std::cout<<std::endl;
    }
}

void PrintMatrix(const Matrix &m, int width = 10, double eps = -1){
    m.Print(width, eps);
}

Matrix operator +(const Matrix &m1, const Matrix &m2){
    if(m1.NCols()!=m2.NCols() || m1.NRows()!=m2.NRows()) 
        throw std::domain_error("Incompatible formats");
    /*int brKol=m1.NCols(), brRed=m1.NRows();
    Matrix rez(brRed, brKol);
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            rez.matrica.at(i).at(j)=m1.matrica.at(i).at(j)+m2.matrica.at(i).at(j);
            //rez(i,j)=m1(i,j)+m2(i, j);
            //ako ces ovako onda u petljama ide od 1 do brRed+1/brKol+1
    return rez;*/
    Matrix rez=m1;
    return rez+=m2;
}

Matrix & Matrix::operator +=(const Matrix &m){
    if(this->NCols()!=m.NCols() || this->NRows()!=m.NRows()) 
        throw std::domain_error("Incompatible formats");
    int brKol=this->NCols(), brRed=this->NRows();
    Matrix rez(brRed, brKol);
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            this->matrica.at(i).at(j)+=m.matrica.at(i).at(j);
            //this->operator()(i,j)+=m(i, j);
            //ako ces ovako onda u petljama ide od 1 do brRed+1/brKol+1
    return *this;
}

Matrix operator -(const Matrix &m1, const Matrix &m2){
    if(m1.NCols()!=m2.NCols() || m1.NRows()!=m2.NRows()) 
        throw std::domain_error("Incompatible formats");
    /*int brKol=m1.NCols(), brRed=m1.NRows();
    Matrix rez(brRed, brKol);
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            rez.matrica.at(i).at(j)=m1.matrica.at(i).at(j)-m2.matrica.at(i).at(j);
            //rez(i,j)=m1(i,j)-m2(i, j);
    return rez;*/
    Matrix rez=m1;
    return rez-=m2;
}

Matrix & Matrix::operator -=(const Matrix &m){
    if(this->NCols()!=m.NCols() || this->NRows()!=m.NRows()) 
        throw std::domain_error("Incompatible formats");
    int brKol=this->NCols(), brRed=this->NRows();
    Matrix rez(brRed, brKol);
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            this->matrica.at(i).at(j)-=m.matrica.at(i).at(j);
            //this->operator()(i,j)-=m(i, j);
    return *this;
}

Matrix operator *(double s, const Matrix &m){
    /*int brKol=m.NCols(), brRed=m.NRows();
    auto rez=m;
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            rez.matrica.at(i).at(j)*=s;
            //rez(i,j)*=s;
    return rez;*/
    Matrix rez=m;
    return rez*=s;
}

Matrix operator *(const Matrix &m, double s){
    return s*m;
}

Matrix & Matrix::operator *=(double s){
    int brKol=this->NCols(), brRed=this->NRows();
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++) 
            this->matrica.at(i).at(j)*=s;
            //rez(i,j)*=s; - ali pazi na granice petlji
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2){
    if(m1.NCols()!=m2.NRows()) throw std::domain_error("Incompatible formats");
    /*Matrix nova(m1.NRows(), m2.NCols());
    int brRedNove=nova.NRows(), brKolNove=nova.NCols(), ColRed=m1.NCols();
    for (int i = 0; i < brRedNove; i++) 
        for (int j = 0; j < brKolNove; j++) 
            for (int k = 0; k < ColRed; k++) 
                nova.matrica[i][j] += m1.matrica[i][k] * m2.matrica[k][j];
    return nova;*/
    Matrix rez=m1;
    return rez*=m2;
}

Matrix & Matrix::operator *=(const Matrix &m){
    if(this->NCols()!=m.NRows()) throw std::domain_error("Incompatible formats");
    Matrix nova(this->NRows(), m.NCols());
    int brRedNove=nova.NRows(), brKolNove=nova.NCols(), ColRed=this->NCols();
    for (int i = 0; i < brRedNove; i++) 
        for (int j = 0; j < brKolNove; j++) 
            for (int k = 0; k < ColRed; k++) 
                nova.matrica[i][j] += this->matrica[i][k] * m.matrica[k][j];
    *this=nova;
    return *this;
}

Vector operator *(const Matrix &m, const Vector &v){
    if(m.NCols()!=v.NElems()) throw std::domain_error("Incompatible formats");
    Vector novi(m.NRows());
    int brElemenata=m.NRows(), ColRed=v.NElems();
    for (int i = 0; i < brElemenata; i++) 
        for (int j = 0; j < ColRed; j++) 
            novi[i] += m.matrica[i][j] * v[j];
    return novi;
}

Matrix Transpose(const Matrix &m){
    int brRed=m.NRows(), brKol=m.NCols();
    Matrix transponovana(brKol, brRed);
    for(int i=0; i<brKol; i++)
        for(int j=0; j<brRed; j++) 
            transponovana.matrica[i][j] = m.matrica[j][i];
    return transponovana;
}

void Matrix::Transpose(){
    int brRed=this->NRows(), brKol=this->NCols();
    if(brRed!= brKol) {
        Matrix transponovana(brKol, brRed);
        for(int i=0; i<brKol; i++)
            for(int j=0; j<brRed; j++)
                transponovana.matrica[i][j] = this->matrica[j][i];
        *this = transponovana;
        //*this=Transpose(*this); - zasto nece
    }
    else {
        for(int i=0; i<brRed; i++){
            for(int j=i; j<brKol; j++) {
                if(i!=j) {
                    double temp;
                    temp = this->matrica[i][j];
                    this->matrica[i][j] = this->matrica[j][i];
                    this->matrica[j][i] = temp;
                }
            }
        }
    }
}

int main(){
    try{
        //test Vector
        try{ Vector v({});
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;
        }
        try{ Vector v(-1);
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;
        }
        Vector v1({1.11,2.23,3.43}), v2({4.01,5.55,6.78});
        Vector vp(3);
        vp[0]=10.11; vp[1]=11.222; vp(2)++;
        PrintVector(vp, '+'); std::cout<<std::endl;
        std::cout<<"Velicine:\nv1: "<<v1.NElems()<<"\nv2: "<<v2.NElems()<<"\n";
        std::cout<<v1[2];
        v1[2]+=3;
        std::cout<<"\n"<<v1[2]<<"\n";
        std::cout<<v2(2);
        v2(2)+=3.02; std::cout<<"\n"<<v2(2);
        std::cout<<"\nNorma v1 je "<<v1.Norm()<<std::endl;
        v2.Print('*'); std::cout<<std::endl;
        try{
            Vector v(2);
            Vector rez=v+v1;
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;
        }
        try{
            Vector v(2);
            v1-=v;
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;
        }
        try{
            Vector v(2);
            double p=v*v2;
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;
        }
        Vector v3=v1+v2;
        v3.Print('-'); std::cout<<std::endl;
        v3+=v1;
        v3.Print('-'); std::cout<<std::endl;
        v3-=v2;
        v3.Print('-'); std::cout<<std::endl;
        v3*=2.5;
        v3.Print('-'); std::cout<<std::endl;
        std::cout<<"Skalarni proizvod v1 i v2 je "<<v1*v2<<"\n";
        v1/=2;
        v1.Print('*'); std::cout<<std::endl;
        std::cout<<"Epsilon za v1: "<<v1.GetEpsilon()<<std::endl;
        try{v2/=0;}
        catch(std::exception &izuzetak){std::cout<<"Greska Vector: "<<izuzetak.what()<<std::endl;}

        //test Matrix
        try{ Matrix matT1(3,0);
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        try{ Matrix matT2({});
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        try{ Matrix matT3({{1.1,2,3.09}, {}});
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        try{ Matrix matT4({{1,2},{1,2,3}});
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }

        Matrix m1(2,3);
        std::cout<<"m1\n"; m1.Print(); std::cout<<std::endl;
        Matrix m2({{1.11,2.22,3.33},{4.33,5.22,6.11}});
        std::cout<<"m2\n"; PrintMatrix(m2, 3); std::cout<<std::endl;
        Matrix m3({10.23,11.03,12.22});
        std::cout<<"m3\n"; m3.Print(); std::cout<<std::endl;
        int m2Col=m2.NCols(), m2Row=m2.NRows();

        auto m2r0=m2[0];
        for(int i=0; i<m2Col; i++) std::cout<<m2r0[i]<<" ";
        std::cout<<std::endl;
        try{ std::cout<<m2(0,0);
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        std::cout<<"m2[1][1]: "<<m2(1,1);
        m2(1,1)=11;
        std::cout<<", a nakon promjene m2[1][1]: "<<m2(1,1)<<std::endl;
        std::cout<<"Frobeniusova norma m3 iznosi: "<<MatrixNorm(m3)<<std::endl;
        std::cout<<"Epsilon za m3: "<<m3.GetEpsilon()<<std::endl;
        
        try{ auto rez=m2+m3;
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        Matrix m4({{3,2,1},{6,5,4}});
        std::cout<<"m4\n"; PrintMatrix(m4, 3); std::cout<<std::endl;
        auto rez=m2+m4;
        std::cout<<"m2+m4\n"; rez.Print(); std::cout<<std::endl;
        rez-=m2;
        std::cout<<"rez-m2\n"; rez.Print(); std::cout<<std::endl;
        rez=rez*2;
        std::cout<<"rez*2\n"; rez.Print(); std::cout<<std::endl;
        rez=3*rez;
        std::cout<<"3*rez\n"; rez.Print(); std::cout<<std::endl;
        rez*=(1./3);
        std::cout<<"rez*1/3\n"; rez.Print(); std::cout<<std::endl;
        Matrix m5({{1.33,2},{3,4.56},{5.32,6}});
        try{
            Matrix rez1=m2*m4;
            std::cout<<"m2*m4\n"; rez1.Print(); std::cout<<std::endl;
        }
        catch(std::exception &izuzetak){
            std::cout<<"Greska Matrix: "<<izuzetak.what()<<std::endl;
        }
        Matrix rez1=m2*m5;
        std::cout<<"rez1=m2*m5\n"; rez1.Print(); std::cout<<std::endl;
        rez1=Transpose(rez1);
        std::cout<<"rez 1 transponovano\n"; rez1.Print(); std::cout<<std::endl;
        rez1*=m4;
        std::cout<<"rez1*=m4\n"; rez1.Print(); std::cout<<std::endl;
        Vector MxV=m2*v1;
        std::cout<<"Matrica*vektor\n"; MxV.Print(); std::cout<<std::endl;
        rez1.Transpose();
        std::cout<<"rez 1 transponovano\n"; rez1.Print(); std::cout<<std::endl;
    }
    catch(std::exception &izuzetak){
        std::cout<<izuzetak.what()<<std::endl;
    }
    return 0;
}