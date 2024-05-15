#include <matplot/matplot.h>
#include <cmath>
#include <math.h>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <complex>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include "AudioFile/AudioFile.h"
namespace py = pybind11;

int pochodna_sygnalu(std::string name)
{
    using namespace matplot;
    using namespace std;
    AudioFile<double> audioFile;
    vector<double> x;
    vector<double> y;
    audioFile.load (name);
    double samplingRate = audioFile.getSampleRate();
    vector<double> signal = audioFile.samples[0];
    int rozmiar=signal.size();
    if(rozmiar>500) rozmiar=500;
    vector<double> derivative(rozmiar);
    double dt = 1.0 / samplingRate;
    for (int i = 0; i < rozmiar; i++) {
        derivative[i] = (signal[i+1] - signal[i]) / dt;
        x.push_back(i);
    }
    //plot(x, signal)->line_width(2).color("red");
    plot(x, derivative)->line_width(2).color("red");
    xlabel("X");
    ylabel("Y");
    show();
    return 0;
}

std::vector<double> idft(std::vector<std::complex<double>> input)
{
    using namespace matplot;
    using namespace std;
    int N = input.size();
    int K=N;
    vector<double> output(N, 0);
    vector<double> x;
    //output.reserve(K);
    for(int k=0;k<K;k++)
    {
        for(int n = 0; n < N; n++)
        {
            output[k] += real(input[n])*cos((2 * M_PI * k * n) / N) + imag(input[n])*sin((2 * M_PI * k * n) / N);
        }
        output[k] /= N;
        x.push_back(k);
    }
    plot(x, output)->line_width(1).color("red");
    show();
    return output;
}

std::vector<double> cosinus(int freq)
{
    using namespace matplot;
    using namespace std;
    vector<double> x, y;
    double t=1000;
    for(double i = 0; i < 1; i+=0.001)
    {
        x.push_back(i);
        y.push_back(cos(2 * M_PI * freq * i));
    }

    // vector<double> x = linspace(0, a * pi);
    // vector<double> y=transform(x, [](auto x) { return cos(x); });
    plot(x, y)->line_width(1).color("red");
    show();
    return y;
}

std::vector<double> sinus(double freq)
{
    using namespace matplot;
    using namespace std;
    vector<double> x, y;
    double t=1000;
    for(double i = 0; i < 1; i+=0.001)
    {
        x.push_back(i);
        y.push_back(sin(2 * M_PI * freq * i));
    }

    // vector<double> x = linspace(0, a * pi);
    // vector<double> y=transform(x, [](auto x) { return cos(x); });
    plot(x, y)->line_width(1).color("red");
    show();
    return y;
}

int piloksztaltny(int a)
{
    using namespace matplot;
    std::vector<double> x;
    std::vector<double> y;
    for(int i = 0; i < 10; i++)
    {
       x.push_back(a*i);
       y.push_back(1);
       x.push_back(a*i+a);
       y.push_back(-1); 
    }
    plot(x, y)->line_width(1).color("red");
    show();
    return 0;
}

int prostokatny(int a)
{
    using namespace matplot;
    std::vector<double> x;
    std::vector<double> y;
    for(int i = 0; i < 10; i+=2)
    {
       x.push_back(a*i);
       y.push_back(1);
       x.push_back(a*i);
       y.push_back(-1); 
       x.push_back(a*i+a);
       y.push_back(-1);
       x.push_back(a*i+a);
       y.push_back(1);
       x.push_back(a*i+2*a);
       y.push_back(1);
    }
    plot(x, y)->line_width(1).color("red");
    show();
    return 0;
}

std::vector<std::complex<double>> dft(std::vector<double> input)
{
    using namespace std;
    using namespace matplot;
    int N = input.size();
    int K=N;
    vector<complex<double>> output(N, 0);
    //output.reserve(K);
    for(int k=0;k<K;k++)
    {
        for(int n = 0; n < N; n++)
        {
            output[k] += input[n] * (cos((2 * M_PI * k * n) / N) - 1i * sin((2 * M_PI * k * n) / N));
        }
    }
    return output;
}

std::vector<std::complex<double>> signals(std::vector<double> input)
{
    using namespace std;
    using namespace matplot;
    std::vector<double> x;
    std::vector<double> y;

    std::vector<std::complex<double>> Fx=dft(input);
    std::vector<double> z;
    //cout<< "\n"<<"k\t" << setw(12) <<"Real\t"<< setw(12) <<"Imag\t"<< setw(12) <<"signal"<<endl;
    int N = Fx.size();
    for(int i=0;i<=N;i++)
    {
        z.push_back(sqrt(Fx[i].real()*Fx[i].real() + Fx[i].imag()*Fx[i].imag()));
        //cout<<i<<"\t" <<setw(12)<<Fx[i].real() <<"\t" <<setw(12)<<Fx[i].imag() <<"\t" <<setw(12)<<signal[i] <<endl;
        x.push_back(i);
        y.push_back(z[i]);
    }
    
    plot(x, y)->line_width(1).color("red");
    show();
    return Fx;
}


PYBIND11_MODULE(tp_3, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("pochodna_sygnalu", &pochodna_sygnalu,"Funkcja generujaca pochodna sygnalu");
    m.def("idft", &idft,"Funkcja generujaca pochodna sygnalu");
    m.def("cosinus", &cosinus, "Funkcja generujaca wykres cosinusa");
    m.def("sinus", &sinus, "Funkcja generujaca wykres sinusa");
    m.def("prostokatny", &prostokatny, "Funkcja generujaca wykres prostokatny");
    m.def("piloksztaltny", &piloksztaltny, "Funkcja generujaca wykres piloksztaltny");
    m.def("signals", &signals, "A function that adds two numbers");

    py::class_<std::vector<double>>(m, "vector_float")
    .def(py::init<>())
    .def("clear", &std::vector<double>::clear)
    .def("reserve", &std::vector<double>::reserve)
    .def("resize", py::overload_cast<size_t>(&std::vector<double>::resize))
    .def("push_back", py::overload_cast<const double&>(&std::vector<double>::push_back))
    .def("pop_back", &std::vector<double>::pop_back)
    .def("__len__", [](const std::vector<double> &v) { return v.size(); }) 
    .def("__iter__", [](std::vector<double> &v) {
        return py::make_iterator(v.begin(), v.end());} , py::keep_alive<0, 1>());

    py::class_<std::vector<std::complex<double>>>(m, "vector_complex")
    .def(py::init<>())
    .def("clear", &std::vector<std::complex<double>>::clear)
    .def("reserve", &std::vector<std::complex<double>>::reserve)
    .def("resize", py::overload_cast<size_t>(&std::vector<std::complex<double>>::resize))
    //.def("push_back", py::overload_cast<const double&>(&std::vector<std::complex<double>>::push_back))
    .def("pop_back", &std::vector<std::complex<double>>::pop_back)
    .def("__len__", [](const std::vector<std::complex<double>> &v) { return v.size(); })
    .def("__iter__", [](std::vector<std::complex<double>> &v) {
        return py::make_iterator(v.begin(), v.end());} , py::keep_alive<0, 1>());
    
    // py::class_<std::complex<double>>(m, "vector_complex")
    // .def(py::init<>())
    // .def("real", &std::complex<double>::real)
    // .def("imag", &std::complex<double>::imag)
    
}