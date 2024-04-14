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
    void Chop(double eps = -1);
    bool EqualTo(const Vector &v, double eps = -1) const;
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
    void Chop(double eps = -1);
    bool EqualTo(const Matrix &m, double eps = -1) const;
    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);
    friend Matrix operator /(const Matrix &m, double s);
    Matrix &operator /=(double s);
    friend Matrix operator /(Matrix m1, Matrix m2);
    Matrix &operator /=(Matrix m);
    double Det() const;
    friend double Det(Matrix m);
    void Invert();
    friend Matrix Inverse(Matrix m);
    void ReduceToRREF();
    friend Matrix RREF(Matrix m);
    int Rank() const;
    friend int Rank(Matrix m);
};

class LUDecomposer{
  private:
    Matrix kompaktnaLU;
    Vector permutacije;
  public:
    LUDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(const Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Matrix GetCompactLU() const;
    Matrix GetL() const;
    Matrix GetU() const;
    Vector GetPermuation() const;
};

class QRDecomposer{
   private:
    Matrix R;
    Vector diagR;
   public:
    QRDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Vector MulQWith(Vector v) const;
    Matrix MulQWith(Matrix m) const;
    Vector MulQTWith(Vector v) const;
    Matrix MulQTWith(Matrix m) const;
    Matrix GetQ() const;
    Matrix GetR() const;
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
    return vektor[i];
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
    //int vel=NElems();
    for(double v:this->vektor) sumaKomponenti+=(v*v);
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
}
    
Vector & Vector::operator +=(const Vector &v){
    if(this->NElems()!=v.NElems()) throw std::domain_error("Incompatible formats");
    int vel=this->NElems();
    for(int i=0; i<vel; i++) this->vektor[i]+=v.vektor[i];
    return *this;
}

Vector operator -(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
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
    Vector v1=v;
    return v1*=s;
}

Vector operator *(const Vector &v, double s){
    return s*v;
}

Vector & Vector::operator *=(double s){
    int vel=this->NElems();
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
    Vector v1=v;
    return v1/=s;
}

Vector & Vector::operator /=(double s){
    if(std::fabs(s-0)<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    int vel=this->NElems();
    for(int i=0; i<vel; i++) this->vektor[i]/=s;
    return *this;
}

void Vector::Chop(double eps){
    if(eps<0) eps=this->GetEpsilon();
    int vel=this->NElems();
    for(int i=0; i<vel; i++){
        if(std::fabs(this->vektor[i])<eps) this->vektor[i]=0;
    }
}

bool Vector::EqualTo(const Vector &v, double eps) const{
    if(this->vektor.size()!=v.NElems()) return false;
    if(eps<0) eps=this->GetEpsilon();
    int vel=this->NElems();
    for(int i=0; i<vel; i++){
        if(std::fabs(this->vektor[i]-v[i])>eps) return false;
    }
    return true;
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
    int brRed=this->NRows(), brKol=this->NCols();
    for(int i=0; i<brRed; i++){
        for(int j=0; j<brKol; j++) {
            std::cout.width((*this)[i][j] <0 ? width+1 : width);
            std::cout<<(std::fabs((*this)[i][j]) <eps ? 0 : (*this)[i][j]);
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
    return *this;
}

Matrix operator -(const Matrix &m1, const Matrix &m2){
    if(m1.NCols()!=m2.NCols() || m1.NRows()!=m2.NRows()) 
        throw std::domain_error("Incompatible formats");
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
    return *this;
}

Matrix operator *(double s, const Matrix &m){
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
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2){
    if(m1.NCols()!=m2.NRows()) throw std::domain_error("Incompatible formats");
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

void Matrix::Chop(double eps){
    if(eps<0) eps=this->GetEpsilon();
    int velR=this->NRows(), velK=this->NCols();
    for(int i=0; i<velR; i++)
        for(int j=0; j<velK; j++)
            if(std::fabs(this->matrica[i][j])<eps) this->matrica[i][j]=0;
}

bool Matrix::EqualTo(const Matrix &m, double eps) const{
    int velR=this->NRows(), velK=this->NCols();
    if(velR!=m.NRows() || velK!=m.NCols()) return false;
    if(eps<0) eps=this->GetEpsilon();
    for(int i=0; i<velR; i++)
        for(int j=0; j<velK; j++)
            if(std::fabs(this->matrica[i][j]-m[i][j])>eps) return false;
    return true;
}

Matrix LeftDiv(Matrix m1, Matrix m2){
    if(m1.NCols()!=m1.NRows()) 
        throw std::domain_error("Divisor matrix is not square");
    if(m1.NRows()!=m2.NRows())
        throw std::domain_error("Incompatible formats");
    int velM1=m1.NRows(), brKolM2=m2.NCols(), brKolM1=m1.NCols();
    double eps=m1.GetEpsilon();

    for(int k=0; k<brKolM1; k++){
        int p=k;
        for(int i=k+1; i<velM1; i++){
            if(std::fabs(m1[i][k])>std::fabs(m1[p][k]))
                p=i;
        }
        if(std::fabs(m1[p][k])<eps)
            throw std::domain_error("Divisor matrix is singular");
        if(p!=k){
            std::swap(m1.matrica[k], m1.matrica[p]);
            std::swap(m2.matrica[k], m2.matrica[p]);
        }
        for(int i=k+1; i<brKolM1; i++){
            long double u=m1[i][k]/m1[k][k];
            for(int j=k+1; j<brKolM1; j++)
                m1[i][j]-=u*m1[k][j];
            for(int j=0; j<brKolM2; j++)
                m2[i][j]-=u*m2[k][j];
        }
    }
    //supstitucija unazad
    Matrix rez(velM1, brKolM2);
    for(int k=0; k<brKolM2; k++){
        for(int i=velM1-1; i>=0; i--){
            long double s=m2[i][k];
            for(int j=i+1; j<velM1; j++)
                s=s-m1[i][j]*rez[j][k];
            rez[i][k]=s/m1[i][i];
        }
    }
    return rez;
}

Vector LeftDiv(Matrix m, Vector v){
    if(m.NCols()!=m.NRows()) 
        throw std::domain_error("Divisor matrix is not square");
    if(m.NRows()!=v.NElems())
        throw std::domain_error("Incompatible formats");
    int velM=m.NRows();
    double eps=m.GetEpsilon();

    for(int k=0; k<velM; k++){
        int p=k;
        for(int i=k+1; i<velM; i++){
            if(std::fabs(m[i][k])>std::fabs(m[p][k]))
                p=i;
        }
        if(std::fabs(m[p][k])<eps)
            throw std::domain_error("Divisor matrix is singular");
        if(p!=k){
            std::swap(m.matrica[k], m.matrica[p]);
            std::swap(v[k], v[p]);
        }
        for(int i=k+1; i<velM; i++){
            long double u=m[i][k]/m[k][k];
            for(int j=k+1; j<velM; j++)
                m[i][j]-=u*m[k][j];
            v[i]-=u*v[k];
        }
    }
    //supstitucija unazad
    Vector rez(velM);
    for(int i=velM-1; i>=0; i--){
        long double s=v[i];
        for(int j=i+1; j<velM; j++)
            s-=m[i][j]*rez[j];
        rez[i]=s/m[i][i];
    }
    return rez;
}

Matrix operator /(const Matrix &m, double s){
    Matrix rez=m;
    try{
        return rez/=s;
    }
    catch(std::domain_error &e){
        throw e;
    }
}

Matrix & Matrix::operator /=(double s){
    if(std::fabs(s-0)<std::numeric_limits<double>::epsilon()) 
        throw std::domain_error("Division by zero");
    int brRed=this->NRows(), brKol=this->NCols();
    for(int i=0; i<brRed; i++)
        for(int j=0; j<brKol; j++)
            this->matrica[i][j]/=s;
    return *this;
}

Matrix operator /(Matrix m1, Matrix m2){
    try{
        return m1/=m2;
    }
    catch(std::domain_error &e){
        throw e;
    }
}

Matrix & Matrix::operator /=(Matrix m){
    if(m.NCols()!=m.NRows()) 
        throw std::domain_error("Divisor matrix is not square");
    if(this->NCols()!=m.NRows())
        throw std::domain_error("Incompatible formats");
    int velM1=m.NRows(), brRedM2=this->NRows(), brKolM1=m.NCols();
    double eps=this->GetEpsilon();

    for(int k=0; k<velM1; k++){
        int p=k;
        for(int i=k+1; i<velM1; i++){
            if(std::fabs(m[i][k])>std::fabs(m[p][k]))
                p=i;
        }
        if(std::fabs(m[p][k])<eps)
            throw std::domain_error("Divisor matrix is singular");
        if(p!=k){
            for (int i = 0; i < velM1; i++) {
                std::swap(m[i][k], m[i][p]);
            }
            for (int i = 0; i < brRedM2; i++) {
                std::swap((*this)[i][k], (*this)[i][p]);
            }
        }
        for(int i=k+1; i<brKolM1; i++){
            long double u=m[k][i]/m[k][k];
            for(int j=k+1; j<velM1; j++)
                m[j][i]-=u*m[j][k]; 
            for(int j=0; j<brRedM2; j++)
                (*this)[j][i]-=u*(*this)[j][k];
        }
    }
    //supstitucija unazad
    Matrix rez(brRedM2, velM1);
    for (int k = 0; k < brRedM2; k++) {
        for (int i = brKolM1 - 1; i >= 0; i--) {
            long double s=(*this)[k][i];
            for (int j = i + 1; j < velM1; j++) {
                s-=m[j][i] * rez[k][j];
            }
            rez[k][i]=s/ m[i][i];
        }
    }
    *this=rez;
    return *this;
}

double Matrix::Det() const{
    if(this->NRows()!=this->NCols())
        throw std::domain_error("Matrix is not square");
    long double det=1;
    int brRed=this->NRows();
    double eps=this->GetEpsilon();
    Matrix pomocna=(*this);

    for(int k=0; k<brRed; k++){
        int p=k;
        for(int i=k+1; i<brRed; i++){
            if(std::fabs(pomocna[i][k])>std::fabs(pomocna[p][k]))
                p=i;
        }
        if(std::fabs(pomocna[p][k])<eps)
            //throw std::domain_error("Matrix is singular");
            return 0;
        if(p!=k){
            std::swap(pomocna.matrica[k],pomocna.matrica[p]);
            det*=-1;
        }
        det=det*pomocna[k][k];
        for(int i=k+1; i<brRed; i++){
            long double u=pomocna[i][k]/pomocna[k][k];
            for(int j=k+1; j<brRed; j++)
                pomocna[i][j]-=u*pomocna[k][j];
        }
    }
    return det;
}

double Det(Matrix m){
    if(m.NRows()!=m.NCols())
        throw std::domain_error("Matrix is not square");
    try{
        return m.Det();
    }
    catch(std::domain_error &e){
        throw e;
    }
}

void Matrix::Invert(){
    if(this->NRows()!=this->NCols())
        throw std::domain_error("Matrix is not square");
    int brRed=this->NRows();
    double eps=this->GetEpsilon();
    std::vector<int> pivot(brRed);

    for(int k=0; k<brRed; k++){
        int p=k;
        for(int i=k+1; i<brRed; i++){
            if(std::fabs(matrica[i][k])>std::fabs(matrica[p][k]))
                p=i;
        }
        if(std::fabs(matrica[p][k])<eps)
            throw std::domain_error("Matrix is singular");
        if(p!=k){
            std::swap(matrica[k],matrica[p]);
        }
        pivot[k]=p;
        long double u=matrica[k][k];
        matrica[k][k]=1.0;
        for(int j=0; j<brRed; j++)
            matrica[k][j]/=u;
        for(int i=0; i<brRed; i++){
            if(i!=k){
                u=matrica[i][k];
                matrica[i][k]=0.0;
                for(int j=0; j<brRed; j++)
                    matrica[i][j]-=u*matrica[k][j];
            }
        }
    }
    for(int j=brRed-1; j>=0; j--){
        int p=pivot[j];
        if(p!=j){
            for(int i=0; i<brRed; i++){
                std::swap(matrica[i][j], matrica[i][p]);
            }
        }
    }
    
}

Matrix Inverse(Matrix m){
    if(m.NRows()!=m.NCols())
        throw std::domain_error("Matrix is not square");
    try{
        m.Invert();
        return m;
    }
    catch(std::domain_error &e){
        throw e;
    }
}

void Matrix::ReduceToRREF(){
    int k=-1, l=-1, brRed=this->NRows(), brKol=this->NCols();
    std::vector<bool> w(brKol, false);
    double eps=this->GetEpsilon();

    while(k<brRed && l<brKol){
        l++; k++;
        double v=0;
        int p=0;
        while(v<eps && l<brKol){
            p=k;
            for(int i=k; i<brRed; i++){
                if(std::fabs(matrica[i][l])>v){
                    v=std::fabs(matrica[i][l]);
                    p=i;
                }
            }
            if(v<eps) l++;
        }
        if(l<brKol){
            w[l]=true;
            if(p!=k) std::swap(matrica[p], matrica[k]);
            long double u=matrica[k][l];
            for(int j=l; j<brKol; j++)
                matrica[k][j]/=u;
            for(int i=0; i<brRed; i++){
                if(i!=k){
                    u=matrica[i][l];
                    for(int j=l; j<brKol; j++)
                        matrica[i][j]-=u*matrica[k][j];
                }
            }
        }
    }
}

Matrix RREF(Matrix m){
    m.ReduceToRREF();
    return m;
}

int Matrix::Rank() const{
    Matrix pomocna=*this;
    pomocna.ReduceToRREF();
    int rank = 0, brRed=this->NRows(), brKol=this->NCols();
    double eps=this->GetEpsilon();

    for (int i=0; i<brRed; i++) {
        int j=i;
        while(j<brKol && pomocna[i][j]<std::numeric_limits<double>::epsilon()) 
            j++;
        if(j==brKol) break;
        rank++;
    }
    return rank;
}

int Rank(Matrix m){
    return m.Rank();
}

//LUDecoposer implementacije
LUDecomposer::LUDecomposer(Matrix m) : kompaktnaLU(m.NRows(),m.NCols()), permutacije(m.NRows()){
    if(m.NRows()!=m.NCols())
        throw std::domain_error("Matrix is not square");
    int brRedM=m.NRows(), brKol=m.NCols();
    double eps=m.GetEpsilon();

    for(int j=0; j<brRedM; j++){
        for(int i=0; i<=j; i++){
            long double s=m[i][j];
            for(int k=0; k<i; k++)
                s-=m[i][k]*m[k][j];
            m[i][j]=s;
        }
        int p=j;
        for(int i=j+1; i<brRedM; i++){
            long double s=m[i][j];
            for(int k=0; k<j; k++)
                s-=m[i][k]*m[k][j];
            m[i][j]=s;
            if(std::fabs(s)>std::fabs(m[p][j])) p=i;
        }
        if(std::fabs(m[p][j])<eps)
            throw std::domain_error("Matrix is singular");
        if(p!=j)
            for(int l=0; l<brKol; l++){
                double temp=m[j][l];
                m[j][l]=m[p][l];
                m[p][l]=temp;
            }
        permutacije[j]=p;
        long double u=m[j][j];
        for(int i=j+1; i<brRedM; i++)
            m[i][j]/=u;
    }
    kompaktnaLU=m;
}

void LUDecomposer::Solve(const Vector &b, Vector &x) const{
    if (b.NElems() != kompaktnaLU.NRows() || x.NElems() != kompaktnaLU.NRows())
        throw std::domain_error("Incompatible formats");

    x=Solve(b);
}

Vector LUDecomposer::Solve(Vector b) const{
    if (b.NElems() != kompaktnaLU.NRows())
        throw std::domain_error("Incompatible formats");

    int vel= b.NElems();
    Vector y(vel), x(vel);
    Matrix l(GetL()), u(GetU());

    // Rješavanje Ly = Pb
    for (int i = 0; i < vel; i++) {
        int p=permutacije[i];
        long double s=b[p];
        b[p]=b[i];
        for (int j = 0; j < i; j++) {
            s-= l[i][j] * y[j];
        }
        y[i]=s;
    }

    // Rješavanje Ux = y
    for (int i = vel - 1; i >= 0; i--) {
        long double s=y[i];
        for (int j = i + 1; j < vel; j++) {
            s -= u[i][j] * x[j];
        }
        x[i]= s/u[i][i];
    }
    return x;
}

void LUDecomposer::Solve(const Matrix &b, Matrix &x) const{
    if (b.NRows() != kompaktnaLU.NRows() || b.NCols() != x.NCols() 
        || b.NRows() != x.NRows())
        throw std::domain_error("Incompatible formats");
    x=Solve(b);
}

Matrix LUDecomposer::Solve(Matrix b) const{
    if (b.NRows() != kompaktnaLU.NRows())
        throw std::domain_error("Incompatible formats");

    int brRed = b.NRows();
    int brKol = b.NCols();
    Matrix rez(brRed, brKol);
    // Rješavanje za svaku kolonu matrice B
    for (int kol = 0; kol < brKol; kol++) {
        Vector bKol(brRed);
        for (int red = 0; red < brRed; red++)
            bKol[red] = b[red][kol];
        Solve(bKol, bKol);

        for (int red = 0; red <brRed; red++)
            rez[red][kol] = bKol[red];
    }
    return rez;
}
    
Matrix LUDecomposer::GetCompactLU() const{
    return kompaktnaLU;
}

Matrix LUDecomposer::GetL() const{
    int brRed = kompaktnaLU.NRows();
    Matrix L(brRed, brRed);
    for (int i = 0; i < brRed; i++) {
        L[i][i] = 1.0;
        for (int j = 0; j < i; j++) {
            L[i][j] = kompaktnaLU[i][j];
        }
    }
    return L;
}

Matrix LUDecomposer::GetU() const{
    int brRed = kompaktnaLU.NRows();
    Matrix U(brRed, brRed);
    for (int i = 0; i < brRed; i++)
        for (int j = i; j < brRed; j++)
            U[i][j] = kompaktnaLU[i][j];
    return U;
}

Vector LUDecomposer::GetPermuation() const{
    int vel=permutacije.NElems();
    Vector p(permutacije);
    for(int i=0; i<vel; i++)
        p[i]++;
    return p;
}


//QRDecomposer implementacije
QRDecomposer::QRDecomposer(Matrix m): R(m.NRows(), m.NCols()), diagR(m.NCols()){
    if (m.NRows() < m.NCols())
        throw std::domain_error("Invalid matrix format");

    int brKol=m.NCols(), brRed=m.NRows();
    double eps=m.GetEpsilon();

    for (int k = 0; k < brKol; k++) {
        long double norma = 0.0;
        for (int i = k; i < brRed; i++) {
            norma += m[i][k] * m[i][k];
        }
        norma = std::sqrt(norma);
        long double u=std::sqrt(norma*(norma+std::fabs(m[k][k])));

        if (u < eps)
            throw std::domain_error("Matrix is singular");


        if (m[k][k] < -std::numeric_limits<double>::epsilon())
            norma=-norma;

        m[k][k]=(m[k][k]+norma)/u;

        for (int i = k+1; i < brRed; i++)
            m[i][k]= m[i][k]/u;
        
        diagR[k]=-norma;

        for (int j = k + 1; j < brKol; j++) {
            long double s = 0.0;
            for (int i = k; i < brRed; i++)
                s += m[i][k] * m[i][j];

            for (int i = k; i < brRed; i++)
                m[i][j] -= s * m[i][k];
        }
    }
    R=m;
}

void QRDecomposer::Solve(const Vector &b, Vector &x) const{
    if (R.NRows() != R.NCols()) 
        throw std::domain_error("Matrix is not square");

    if (b.NElems() != R.NRows() || x.NElems()!=R.NRows())
        throw std::domain_error("Incompatible formats");

    x=Solve(b);
}

Vector QRDecomposer::Solve(Vector b) const{
    if (R.NRows() != R.NCols()) 
        throw std::domain_error("Matrix is not square");
        
    if (b.NElems() != R.NRows())
        throw std::domain_error("Incompatible formats");

    int brRed=R.NRows();
    Matrix r(GetR());
    Vector x(brRed);
    b=MulQTWith(b);

    for (int i = brRed - 1; i >= 0; i--) {
        long double s=b[i];
        for (int j = i + 1; j <brRed; j++)
            s-=r[i][j]*x[j];
        x[i]=s/r[i][i];  
    }
    return x;
}

void QRDecomposer::Solve(Matrix &b, Matrix &x) const{
    if (R.NRows() != R.NCols())
        throw std::domain_error("Matrix is not square");

    if (b.NRows() != R.NRows() || x.NRows()!=R.NRows() || x.NCols()!=b.NCols())
        throw std::domain_error("Incompatible formats");

    x = Solve(b);
}

Matrix QRDecomposer::Solve(Matrix b) const{
    if (R.NRows() != R.NCols())
        throw std::domain_error("Matrix is not square");

    if (b.NRows() != R.NRows())
        throw std::domain_error("Incompatible formats");

    int brRed=b.NRows(), brKol=b.NCols();
    Matrix rez(brRed, brKol);

    for (int j = 0; j < brKol; j++) {
        Vector kolonaB(brRed);
        for (int i = 0; i < brRed; i++)
            kolonaB[i] = b[i][j];

        Solve(kolonaB, kolonaB);
        for (int i = 0; i < brRed; i++)
            rez[i][j] = kolonaB[i];
    }
    return rez;
}

Vector QRDecomposer::MulQWith(Vector v) const{
    if (v.NElems() != R.NRows())
        throw std::domain_error("Incompatible formats");
    
    int brRed=R.NRows();
    for (int k = brRed-1; k >=0; k--) {
        long double skalarniProizvod = 0.0;
        for (int i = k; i < brRed; i++) 
            skalarniProizvod += R[i][k] * v[i];

        for (int i = k; i < brRed; i++) 
            v[i] -= skalarniProizvod * R[i][k];
    }
    return v;
}

Matrix QRDecomposer::MulQWith(Matrix m) const{
    if (m.NRows() != R.NRows())
        throw std::domain_error("Incompatible formats");

    int brKolM=m.NCols(), brRedM=m.NRows();
    for (int j = 0; j < brKolM; j++) {
        Vector kolonaM(brRedM);

        for (int i = 0; i < brRedM; i++)
            kolonaM[i] = m[i][j];

        kolonaM = MulQWith(kolonaM);

        for (int i = 0; i < brRedM; i++) {
            m[i][j] = kolonaM[i];
        }
    }
    return m;
}

Vector QRDecomposer::MulQTWith(Vector v) const{
    if (v.NElems() != R.NRows())
        throw std::domain_error("Incompatible formats");

    int brRed=R.NRows();
    for (int k = 0; k <brRed ; k++) {
        long double skalarniProizvod = 0.0;

        for (int i = k; i < brRed; i++)
            skalarniProizvod += R[i][k] * v[i];

        for (int i = k; i < brRed; i++) 
            v[i] -= skalarniProizvod * R[i][k];
    }
    return v;
}

Matrix QRDecomposer::MulQTWith(Matrix m) const{
    if (m.NRows() != R.NRows())
        throw std::domain_error("Incompatible formats");

    int brKolM=m.NCols(), brRedM=m.NRows();
    for (int j = 0; j < brKolM; j++) {
        Vector kolonaM(brRedM);

        for (int i = 0; i < brRedM; i++)
            kolonaM[i] = m[i][j];

        kolonaM = MulQTWith(kolonaM);

        for (int i = 0; i <brRedM; i++)
            m[i][j] = kolonaM[i];
    }
    return m;
}

Matrix QRDecomposer::GetQ() const{
    int brRed=R.NRows(), brKol=R.NCols();
    Matrix q(brRed, brRed);
    for(int j=0; j<brRed; j++){
        for(int i=0; i<brRed; i++)
            q[i][j]=0;
        q[j][j]=1.0;
        for(int k=brKol-1; k>=0; k--){
            long double s=0.0;
            for(int i=k; i<brRed; i++)
                s+=R[i][k]*q[i][j];
            for(int i=k; i<brRed; i++)
                q[i][j]-=s*R[i][k];
        }
    }
    return q;
}

Matrix QRDecomposer::GetR() const{
    int brRed = R.NRows(), brKol=R.NCols();
    Matrix R_vracanje(brRed, brKol);
    for (int j = 0; j < brKol; j++) {
        for (int i = 0; i < j; i++) {
            R_vracanje[i][j] = R[i][j];
        }
        R_vracanje[j][j]=diagR[j];
    }
    return R_vracanje;
}

int main(){
    try{
        Vector v1({1.2, 3, 5, 1e-20, 2.23});
        v1.Chop(); v1.Print(); std::cout<<std::endl;
        Vector v2({1.2, 3, 5, 0, 2.23});
        std::cout<<"Jednakost v1 i v2: "<<v2.EqualTo(v1)<<std::endl;
        Vector v3({1});
        try{
            std::cout<<"Jednakost v1 i v3: "<<v3.EqualTo(v1)<<std::endl<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }

        Matrix m1({{1.2, 1e-20, 2.23}, {1.23, 1.04, 2.23}});
        //PrintMatrix(m1); std::cout<<std::endl;
        m1.Chop(); PrintMatrix(m1); std::cout<<std::endl;
        Matrix m2({{1.2, 0, 2.23}, {1.23, 1.04, 2.23}});
        std::cout<<"Jednakost m1 i m2: "<<m2.EqualTo(m1)<<std::endl;
        Matrix m3(2,3);
        try{
            std::cout<<"Jednakost m2 i m3: "<<m2.EqualTo(m3)<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        //Lijevo dijeljenje
        try{
            m1=LeftDiv(m1,m2);
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        m3=Matrix({{1.2, 0, 2.23}, {1.23, 1.04, 2.23}, {1.2, 0, 2.23}});
        try{
            m1=LeftDiv(m3,m1);
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            m2=Matrix({{8, 1, 5},{-11, 2, 6},{-3, 3, 4}});
            Matrix rez=LeftDiv(m1,m2);
            std::cout<<"\nLijevo matricno m1 i m2:\n";
            rez.Print();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            v1=Vector({8, 1, 5});
            Vector rez=LeftDiv(m1,v1);
            std::cout<<"\nLijevo matricno m1 i v1:\n";
            rez.Print();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            m1/=2.3;
            std::cout<<"\nDijeljenje matrice skalarom:\n";
            m1.Print();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            m1/=0;
            std::cout<<"\nDijeljenje matrice skalarom:\n";
            m1.Print();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            m2=Matrix({{8, 1, 5},{-11, 2, 6},{-3, 3, 4}});
            m1/=m2;
            std::cout<<"\nDesno matricno m1 i m2:\n";
            m1.Print();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2}});
            std::cout<<"\nDeterminanta m1: "<<m1.Det();
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, 1, 2}});
            std::cout<<"\nDeterminanta m1: "<<m1.Det()<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
    try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2}});
            m1.Invert();
            std::cout<<"\nInvertovana m1\n";
            m1.Print(); std::cout<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2}, {1, 2, 3}});
            m1.Invert();
            std::cout<<"\nInvertovana m1\n";
            m1.Print(); std::cout<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1, 1},{-3, -1, 2, 2},{-2, 1, 2, -1}});
            m1.ReduceToRREF();
            std::cout<<"\nRREF m1\n";
            m1.Print(); std::cout<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1, 1},{-3, -1, 2, 2},{-2, 1, 2, -1}});
            std::cout<<"\nRang m1: "<<m1.Rank()<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{1, 2, 3},{4, 5, 6},{7, 8, 9}});
            m1.ReduceToRREF();
            std::cout<<"\nRREF m1\n";
            m1.Print(); std::cout<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{1, 2, 3},{4, 5, 6},{7, 8, 9}});
            std::cout<<"\nRang m1: "<<m1.Rank()<<std::endl;
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }

        //LUDecomposer
        try{
            m1=Matrix({{2, 1, -1, 1},{-3, -1, 2, 2},{-2, 1, 2, -1}});
            LUDecomposer lu1(m1);
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, 1, -1},{-3, -1, 2},{-2, -1, 1}});
            LUDecomposer lu1(m1);
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
            LUDecomposer lu1(m1);
            try{
                Vector b={1,2,3};
                Vector x(3);
                lu1.Solve(b,x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Vector x={1,2,3};
                lu1.Solve(x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        try{
            m1=Matrix({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
            LUDecomposer lu1(m1);
            try{
                Matrix b({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
                Matrix x(3,3);
                lu1.Solve(b,x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix x({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
                lu1.Solve(x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"Kompaktna LU za m1\n";
                m2=lu1.GetCompactLU(); 
                m2.Print();
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"L za m1\n";
                m2=lu1.GetL();
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"U za m1\n";
                m2=lu1.GetU();
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"Vektor permutacija za m1\n";
                v1=lu1.GetPermuation();
                v1.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
        //QRDecomposer
        try{
            try{
                m1=Matrix({{2, -1, 0}, {-1, 2, -1}});
                QRDecomposer qr1(m1);
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            m1=Matrix({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
            QRDecomposer qr1(m1);
            try{
                Matrix b({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
                Matrix x(3,3);
                qr1.Solve(b,x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix x({{2, -1, 0}, {-1, 2, -1},{0, -1, 2}});
                qr1.Solve(x);
                std::cout<<"Rjesnje sistema je:\n";
                x.Print(); std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
                Vector b = {1, 2, 3};
                QRDecomposer qr2(A);
                std::cout<<"MulQWith b vektor\n";
                v1=qr2.MulQWith(b); 
                v1.Print();
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
                Matrix B = {{1, 2}, {3, 4}, {5, 6}};

                QRDecomposer qr2(A);
                std::cout<<"MulQWith b matrica\n";
                m2 = qr2.MulQWith(B);
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
                Vector b = {1, 2, 3};
                QRDecomposer qr2(A);
                std::cout<<"MulQTWith b vektor\n";
                v1=qr2.MulQTWith(b); 
                v1.Print();
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                Matrix A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
                Matrix B = {{1, 2}, {3, 4}, {5, 6}};

                QRDecomposer qr2(A);
                std::cout<<"MulQTWith b matrica\n";
                m2 = qr2.MulQTWith(B);
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"Matrica Q\n";
                m2 = qr1.GetQ();
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
            try{
                std::cout<<"Matrica R\n";
                m2 = qr1.GetR();
                m2.Print(); 
                std::cout<<std::endl;
            }
            catch(std::exception &e){
                std::cout<<e.what()<<std::endl;
            }
        }
        catch(std::exception &e){
            std::cout<<e.what()<<std::endl;
        }
    }
    catch(std::exception &izuzetak){
        std::cout<<izuzetak.what()<<std::endl;
    }
    //testiranje kad je A[1][1]=0;
    try{
        Matrix m1({{0,1.23,4.23},{2, 3.4, 3.8}, {1.09, 5.34, 2.22}});
        try{
            Matrix m2=Inverse(m1); //Gauss-Jordan
            LUDecomposer lu1(m1); //LU
            std::cout<<"Testovi sa a[1][1]=0 prosli\n";
        }
        catch(std::exception &e){
            throw e;
        }
    }
    catch(std::exception &e){
        std::cout<<e.what()<<std::endl;
    }
    return 0;
}