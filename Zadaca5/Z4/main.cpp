#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

bool stepenDvice(unsigned int n) {
    return n!=0 && (n & (n - 1)) == 0;
}

bool prirodanBroj(double n){
    return n!=0 && n==static_cast<unsigned int>(n);
}

std::complex<double> dajW(int n){
    const double pi=4*std::atan(1);
    return std::exp(2*pi*std::complex<double>(0,1)/std::complex<double>(n));
}

void FFT(std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y, int n, int s=0, int d=0, int t=1) {
    if (n == 1) y[d]=x[s];
    else{
        FFT(x, y, n/2, s, d, 2*t);
        FFT(x, y, n/2, s+t, d+n/2, 2*t);

        std::complex<double> m(1), w(std::pow(dajW(n),-1)),u,v;

        for (int k = d; k < d+n/2; k++) {
            u=y[k];
            v=m*y[k+n/2];
            y[k]=u+v;
            y[k+n/2]=u-v;
            m*=w;
        }
    }
}

void inverseFFT(std::vector<std::complex<double>>& y, std::vector<std::complex<double>>& x, int n, int s=0, int d=0, int t=1){
    if (n == 1) x[d]=y[s];
    else{
        inverseFFT(y, x, n/2, s, d, 2*t);
        inverseFFT(y, x, n/2, s+t, d+n/2, 2*t);

        std::complex<double> m(1), w(dajW(n)),u,v;

        for (int k = d; k < d+n/2; k++) {
            u=x[k];
            v=m*x[k+n/2];
            x[k]=(u+v)/2.;
            x[k+n/2]=(u-v)/2.;
            m*=w;
        }
    }
}

std::vector<double> LossyCompress(std::vector<double> data, int new_size) {
    int vel=data.size();
    if (new_size <2 || new_size > data.size()) 
        throw std::range_error("Bad new size");
    if (!stepenDvice(data.size())) 
        throw std::range_error("Data size must be a power of two");

    std::vector<std::complex<double>> complexData(vel), complexCompressed(vel);
    
    for(int i=0; i<vel/2; i++) 
        complexData[i]=std::complex<double>(data[2*i],0.0);
    for(int i=vel/2; i<vel; i++)
        complexData[i]=std::complex<double>(data[2*(vel-i)-1],0.0);

    FFT(complexData, complexCompressed, vel);

    std::vector<double> kompresovano(new_size);
    std::complex<double> w(dajW(2*vel));

    for(int i=0; i<new_size-1; i++)
        kompresovano[i]=(std::pow(w, -i/2.)*complexCompressed[i]).real();

    kompresovano.back()=vel;

    return kompresovano;
}

std::vector<double> LossyDecompress(std::vector<double> compressed) {
    int vel=static_cast<int>(compressed.back());

    if (!stepenDvice(vel) || !prirodanBroj(vel) || vel < compressed.size()-1)
        throw std::logic_error("Bad compressed sequence");
    
    std::vector<double> podaci(compressed);
    podaci.back()=0;
    podaci.resize(vel);

    std::vector<std::complex<double>> complexCompressed(vel), complexData(vel);
    std::complex<double> w(dajW(2*vel));

    complexCompressed.front()=podaci.front();
    for(int i=1; i<vel; i++)
        complexCompressed[i]=2. * std::pow(w, i/2.)*podaci[i];
    
    inverseFFT(complexCompressed, complexData, vel);

    for(int i=0; i<vel; i+=2)
        podaci[i]=complexData[i/2].real();
    for(int i=1; i<vel; i+=2)
        podaci[i]=complexData[vel-(i+1)/2].real();
    
    return podaci;
}

int main() {
    const double pi=4*std::atan(1);
    try{
        //Ne treba baciti izuzetak
        std::vector<double> originalniPodaci;
        int N = 256;
        for (int i = 0; i < N; i++) 
            originalniPodaci.push_back(std::sin(2. * pi * i / N));

        int novaVel = 64;
        std::vector<double> kompresovano = LossyCompress(originalniPodaci, novaVel);
        std::vector<double> dekompresovano = LossyDecompress(kompresovano);

        std::cout<<("Poredjenje originalnih i dekompresovanih podataka:\n");
        for (int i = 0; i < N; i++)
            std::cout << "Original: " << originalniPodaci[i] << ", Dekompresovano: " << dekompresovano[i] << std::endl;
    }
    catch(std::exception &e){
        std::cout<<"Greska: "<<e.what()<< std::endl;
    }
    std::cout<<std::endl;
    try{
        //pogresna velicina kompresije
        std::vector<double> originalniPodaci(16, 0.0);
        int novaVel = 20;
        std::vector<double> kompresovano = LossyCompress(originalniPodaci, novaVel);
    }
    catch(std::range_error &e){
        std::cout<<"Greska: "<<e.what()<< std::endl;
    }
    std::cout<<std::endl;
    try{
        //nije stepen broja 2
        std::vector<double> originalniPodaci(15, 0.0);
        int novaVel = 8;
        std::vector<double> kompresovano = LossyCompress(originalniPodaci, novaVel);
    }
    catch(std::range_error &e){
        std::cout<<"Greska: "<<e.what()<< std::endl;
    }
    std::cout<<std::endl;
    try {
        //pogresna velicina dekompresije
        std::vector<double> kompresovano(10, 0.0);
        std::vector<double> dekompresovano = LossyDecompress(kompresovano);
    } 
    catch (std::logic_error &e) {
        std::cout<<"Greska: "<<e.what()<< std::endl;
    }

    return 0;
}