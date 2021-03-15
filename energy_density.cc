#include"su2.hh"
#include"gaugeconfig.hh"
#include"energy_density.hh"
#include"tensors.hh"
#include<complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// symmetric definition of the energy density
// using the clover leaf
//
//         
//  <--   <--
// |  ^  |  ^
// |  |  |  | mu
// -->   -->
//
//  <--   <--
// |  ^  |  ^
// |  |  |  |
// -->   -->
//        nu
//
// checked for gauge invariance!
//
// energy density
// E = 1/4 G_{mu nu}^a G_{mu nu}^a = 1/2 tr(G_{mu nu} G_{mu nu})
//
// topological charge
// Q = 1./(32 pi^2) exp_\mu\nu\rho\sigma Trace[ G_\mu\nu G_\rho\sigma]
//
// from hep-lat/9603008 we take equation (6)
// G_\mu\nu = 1/4 sum_clover 1/2 (clover - h.c.)
//

void energy_density(gaugeconfig<su2> &U, double &res, double &Q) {
  res = 0.;
  Q = 0.;

  // Euclidean 4D totally anti-symemtric tensor 
  static epsilon4_t eps4 = new_epsilon4();
  
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> x1 = x;
          std::vector<size_t> x2 = x;
          std::vector<size_t> x3 = x;
          su2 G[4][4];
          for(size_t mu = 0; mu < 3; mu++) {
            for(size_t nu = mu+1; nu < 4; nu++) {
              x1[mu] += 1;
              x2[nu] += 1;
              su2 leaf = U(x, mu) * U(x1, nu) *
                U(x2, mu).dagger()*U(x, nu).dagger();
              x1[mu] -= 1;
              x2[nu] -= 1;

              x1[mu] -= 1;
              x1[nu] += 1;
              x2[mu] -= 1;
              leaf += U(x, nu) * U(x1, mu).dagger() *
                U(x2, nu).dagger()*U(x2, mu);
              x1[mu] += 1;
              x1[nu] -= 1;
              x2[mu] += 1;

              x1[mu] -= 1;
              x2[mu] -= 1;
              x2[nu] -= 1;
              x3[nu] -= 1;
              leaf += U(x1, mu).dagger() * U(x2, nu).dagger() *
                U(x2, mu)*U(x3, nu);
              x1[mu] += 1;
              x2[mu] += 1;
              x2[nu] += 1;
              x3[nu] += 1;
              
              x1[nu] -= 1;
              x2[nu] -= 1;
              x2[mu] += 1;
              leaf += U(x1, nu).dagger() * U(x1, mu) *
                U(x2, nu)*U(x, mu).dagger();
              x1[nu] += 1;
              x2[nu] += 1;
              x2[mu] -= 1;

              // traceless and anti-hermitian
              // here we include a factor 1/2 already
              G[mu][nu] =  su2(0.5*(leaf.geta()-std::conj(leaf.geta())), leaf.getb());
              // trace(G_{mu,nu}^a G_{mu,nu}^a)
              // averaged over four plaquette Wilson loops 1./4./4.
              res += trace(G[mu][nu]*G[mu][nu])/16.;
            }
          }

          // sum up the topological charge contribution now
          for( int i = 0; i < eps4.N; i++ ){
            int i1 = eps4.eps_idx[i][0];
            int i2 = eps4.eps_idx[i][1];
            int i3 = eps4.eps_idx[i][2];
            int i4 = eps4.eps_idx[i][3];
            
            // when Gmunu components from the lower triangle are to be used,
            // we can simply skip them and multiply our normalisation by a factor of two
            if( i2 < i1 ){
              continue;
            }
            if( i4 < i3 ){
              continue;
            }
            Q += eps4.eps_val[i]*trace(G[ i1 ][ i2 ]*G[ i3 ][ i4 ] );
          }
        }
      }
    }
  }
  // now we need to devide by 2, but we get a factor of two since we only
  // averaged mu < nu
  res = -res/U.getVolume();
  // factor 2 from summing only mu < nu
  // factor 1/4 from G_\mu\nu
  Q =  -2 * Q  / ( 4 * 32.0 * M_PI * M_PI );
}

void energy_density(gaugeconfig<_u1> &U, double &res, double &Q) {
  res = 0.;
  Q = 0.;

  // Euclidean 4D totally anti-symemtric tensor 
  static epsilon4_t eps4 = new_epsilon4();
  
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> x1 = x;
          std::vector<size_t> x2 = x;
          std::vector<size_t> x3 = x;
          Complex G[4][4];
          for(size_t mu = 0; mu < 3; mu++) {
            for(size_t nu = mu+1; nu < 4; nu++) {
              x1[mu] += 1;
              x2[nu] += 1;
              Complex leaf = 0.;
              leaf += U(x, mu) * U(x1, nu) *
                U(x2, mu).dagger()*U(x, nu).dagger();
              x1[mu] -= 1;
              x2[nu] -= 1;

              x1[mu] -= 1;
              x1[nu] += 1;
              x2[mu] -= 1;
              leaf += U(x, nu) * U(x1, mu).dagger() *
                U(x2, nu).dagger()*U(x2, mu);
              x1[mu] += 1;
              x1[nu] -= 1;
              x2[mu] += 1;

              x1[mu] -= 1;
              x2[mu] -= 1;
              x2[nu] -= 1;
              x3[nu] -= 1;
              leaf += U(x1, mu).dagger() * U(x2, nu).dagger() *
                U(x2, mu)*U(x3, nu);
              x1[mu] += 1;
              x2[mu] += 1;
              x2[nu] += 1;
              x3[nu] += 1;
              
              x1[nu] -= 1;
              x2[nu] -= 1;
              x2[mu] += 1;
              leaf += U(x1, nu).dagger() * U(x1, mu) *
                U(x2, nu)*U(x, mu).dagger();
              x1[nu] += 1;
              x2[nu] += 1;
              x2[mu] -= 1;

              // purely imaginary in U(1)
              // here we include a factor 1/2 already
              G[mu][nu] =  Complex(0, std::imag(leaf));
              // trace(G_{mu,nu}^a G_{mu,nu}^a)
              // averaged over four plaquette Wilson loops 1./4./4.
              res += trace(G[mu][nu]*G[mu][nu])/16.;
            }
          }

          // sum up the topological charge contribution now
          for( int i = 0; i < eps4.N; i++ ){
            int i1 = eps4.eps_idx[i][0];
            int i2 = eps4.eps_idx[i][1];
            int i3 = eps4.eps_idx[i][2];
            int i4 = eps4.eps_idx[i][3];
            
            // when Gmunu components from the lower triangle are to be used,
            // we can simply skip them and multiply our normalisation by a factor of two
            if( eps4.eps_idx[i][1] < eps4.eps_idx[i][0] ){
              continue;
            }
            if( eps4.eps_idx[i][3] < eps4.eps_idx[i][2] ){
              continue;
            }
            Q += eps4.eps_val[i]*trace(G[ i1 ][ i2 ]*G[ i3 ][ i4 ] );
          }
        }
      }
    }
  }
  // now we need to devide by 2, but we get a factor of two since we only
  // averaged mu < nu
  res = -res/U.getVolume();
  // factor 2 from summing only mu < nu
  // factor 1/4 from G_\mu\nu
  Q =  -2 * Q  / ( 4 * 32.0 * M_PI * M_PI );
}
