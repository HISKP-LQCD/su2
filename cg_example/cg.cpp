// cg.cpp

#include <complex>

#include "./LA.hpp"
#include "../cg.hpp"
// #include "./cg_invert.hpp"

typedef std::complex<double> Type;

typedef matrix<Type> LAmatrix;
typedef std::vector<Type> LAvector;

int main(int argc, char const *argv[]) {

        LAmatrix A = (std::vector<LAvector>) {
                {1.0, 2.0, 1.0},
                {2.0, 5.0, -6.0},
                {1.0, -6.0, -6.0}
                };
        
        std::cout << "A=" << '\n';
        cg::print_LAmatrix<Type, LAmatrix, LAvector>(A, ",");
        const std::vector<Type> x_star = {10.0, 11.0, -1.0};
        const std::vector<Type> x0 = {5e+5, -1e+7, 0.0};

        const std::vector<Type> b = A*x_star;


        std::cout << "The solution found by the cg is:" << '\n';

        cg::LinearCG<double, Type, LAmatrix, LAvector>  LCG(A, b);
        LCG.solve(x0, 1e-15, 2);

        std::cout << "The solution should be:" << '\n';
        cg::print_LAvector<Type, LAvector>(x_star, ",");

        // // std::cin.get();
        // std::cout << "Computing the inverse:" << '\n';
        // const int n_hits = 15;
        // cg::inverter<Type, LAmatrix, LAvector> CGI(A, n_hits, 1);
        // CGI.invert(1e-15, 1);
        // std::cout << "The solution found with the inverter is:" << '\n';
        // LAmatrix B = CGI.get_inverse();
        // cg::print_LAvector<Type, LAvector>(B*b, ",");






        return 0;
}
