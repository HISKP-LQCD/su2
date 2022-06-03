// *** EXPERIMENTAL AND DEPRECATED ***
// SOMETIMES THE SOLUYTION CONVERGES TO VECTOR OF NAN

// cg_invert.cpp
// invert a matrix using the cg algorithm

/**
 *  The same assumptions for cg.hpp apply here, with the following additions:
 *  1. LAmatrix::LAmatrix(const size_t& n, const size_t& m); // initialize knowing the number of rows and columns
 *  2. void LAmatrix::resize(const size_t& n, const size_t& m); // resize the matrix
 *  3. void LAmatrix::operator +=(const LAmatrix& B); // A = A+B
 *  4. template<class T> LAmatrix operator /(const LAmatrix& M, const T& lambda); // division by a scalar
*/

#pragma once
#include <random>
#include "../cg.hpp"


namespace cg
{


/** 
 * Inverter class for a matrix M using the cg algorithm for a given number HITS of stochastic sources. 
 * The steps are, for any source:
 * 1. Generate a vector eta_i with random components such that <eta_i>=0 and <eta_i^2>=1
 * 2. Solve numerically M_{ij}*v_{j} = eta_{i} using the cg
 * 3. Return the matrix Q_{ij} = v_{i}*eta_{j}
 * Q is an estimator of M_{ij}^{-1} , 
 * so the inverse returned by this class is the average of the Q obtained for each source.
*/ 
template<class T, class LAmatrix, class LAvector>
class inverter
{
private:
    LAmatrix M; // matrix to be inverted
    size_t HITS; // number of stochastic sources
    size_t SEED; // seed of  the random number generator
    
    std::vector<LAmatrix> Q; // inverse estimators for each stochastic source
    
    bool init_inv = false; // true after initialization is complete
    bool inverted = false; // true after the inversion has been done (for each source)
    
    LAmatrix Q_avr; // average of the Q[h]
    bool averaged = false; // true after the Q[h] have been averaged into Q_avr

    void check_bool(const bool& b) const
    {
        if(!b)
        {
            std::cerr << "Error. cg inverter hasn't been initialized. Aborting.\n";
           std::abort();
        }
    }

    void check_init() const { this->check_bool(init_inv); }
    void check_inverted() const { this->check_bool(inverted); }

public:
    inverter(){};
    ~inverter(){};

    // custom contructor
    inverter(const LAmatrix& m, const size_t& hits, const size_t& seed){ 
        M = m; // matrix to be inverted
        HITS = hits; // number of stochastic sources
        SEED = seed; // seed of the random number generator
        Q.resize(HITS); // storing the space for the HITS estimators of M^{-1}
        init_inv = true; // inverter initialized
    };

    // @tol = tolerance for the solution of the linear system 
    // @verb = verbosity of the solver's output
    void invert(const T& tol, const int& verb){
        this->check_init();

        const size_t N = M.rows();
        LAvector x0(N); // initial guess

        for (size_t h = 0; h < HITS; h++)
        {// loop over the stochasti sources

            if(verb > 0)
            {
                std::cout<< "Stochastic source number: "<< h <<"\n";
            }

            std::normal_distribution<T> dis{0.0, 1.0}; // <eta_i>=0 and <eta_i^2>=1
            std::mt19937 gen(SEED+h); // different seed for each source

            // generating the random vector eta
            Q[h].resize(N, N); // the inverse has the same size
            LAvector eta(N); // random vector
            for (size_t i = 0; i < N; i++){ eta[i] = dis(gen); }
            
            cg::LinearCG<T, LAmatrix, LAvector>  LCG(M, eta); // CG inverter

            // Though not necessary for particular reasons, 
            // we choose to start from the null vector
            print_LAvector<T>(eta);
            LCG.solve(x0, tol, verb); // solving
            LAvector v = LCG.get_solution(); // solution to the system
            
            // filling the elements of the matrices Q[h]
            for (size_t i = 0; i < N; i++){ 
                for (size_t j = 0; j < N; j++){
                    Q[h][i][j] = v[i]*eta[j]; 
                }
            }
        }
        inverted = true; // inversion has been done
        return;
    }

    void average()
    {
        this->check_inverted();
        const size_t N = Q[0].rows();
        Q_avr.resize(N,N);
        for (size_t h = 0; h < HITS; h++)
        {
            Q_avr += Q[h]/HITS; 
        }
        return;
    }


    LAmatrix get_inverse()
    {
        if(!averaged){ this->average(); }
        return Q_avr;
    }

};


}