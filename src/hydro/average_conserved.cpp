//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file average_conserved.cpp
//  \brief Implements polar boundary averaging by Morgan MacLeod.

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../reconstruct/reconstruction.hpp"
#include <stdexcept>

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#include <math.h>
#endif

//MM
int POLAR_AVERAGE=2; // hard-coding this here for now

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::PhiAverageConserved
//  \brief averages conserved quantities in the phi direction

void Hydro::PhiAverageConserved(AthenaArray<Real> &u_in,AthenaArray<Real> &u_out)
{
  // check if we want polar average
  if(do_average_ == false) return;
  // check if we're on the x2-boundary
  MeshBlock *pmb=pmy_block;
  if((pmb->pbval->block_bcs[INNER_X2] == BLOCK_BNDRY) &&
     (pmb->pbval->block_bcs[OUTER_X2] == BLOCK_BNDRY)) {
    return;
  }

  Get_block_N_zone_avg(pmb);
      
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dphi = pmb->pcoord->dx3f(ks); // assume that cell size in phi is constant

  // Mode 1: average v_r, v_th, v_phi each as a scalar (v is momentum)
  // After average, each cell has the same IM1,2,3.
  if (POLAR_AVERAGE == 1) {
    for (int n=0; n<NHYDRO; ++n) {
      for (int j=js; j<=je; ++j) {
        int n_avg_loops = pmb->block_size.nx3/n_avg_(j);
        for (int l=1;l<=n_avg_loops;++l){
          int ks_avg = ks + (l-1)*n_avg_(j);
          int ke_avg = ks_avg + n_avg_(j) -1;
          for (int i=is; i<=ie; ++i) {
            Real u_k_avg = 0.0;
            for (int k=ks_avg; k<=ke_avg; ++k) {
              u_k_avg +=  u_in(n,k,j,i);
            }
            // set the new value
            for (int k=ks_avg; k<=ke_avg; ++k) {
              u_out(n,k,j,i) = u_k_avg/n_avg_(j);
            }
          }
        }
      }
    }
  }

  // Mode 2: average v as a vector
  // After average, each cell has the same v vector (but in general different IM1,2,3).
  else if (POLAR_AVERAGE == 2){
    for (int j=js; j<=je; ++j) {
      int n_avg_loops = pmb->block_size.nx3/n_avg_(j);
      Real sin_th = sin(pmb->pcoord->x2v(j));
      Real cos_th = cos(pmb->pcoord->x2v(j));
      for (int l=1;l<=n_avg_loops;++l){
        int ks_avg = ks + (l-1)*n_avg_(j);
        int ke_avg = ks_avg + n_avg_(j) -1;
        for (int i=is; i<=ie; ++i) {
          Real u_k_avg[NHYDRO]={};
          for (int k=ks_avg; k<=ke_avg; ++k) {
            u_k_avg[IDN] += u_in(IDN,k,j,i);
            u_k_avg[IEN] += u_in(IEN,k,j,i);
            Real phi = dphi * (k-ks_avg - (n_avg_(j)-1)/2.); // angle between the current cell and averaged cell center
            // below physically means u_k_avg += u_in
            u_k_avg[IM1] +=   u_in(IM1,k,j,i)*(1-SQR(sin_th)*(1-cos(phi)))
                            - u_in(IM2,k,j,i)*sin_th*cos_th*(1-cos(phi))
                            - u_in(IM3,k,j,i)*sin_th*sin(phi);
            u_k_avg[IM2] += - u_in(IM1,k,j,i)*sin_th*cos_th*(1-cos(phi))
                            + u_in(IM2,k,j,i)*(1-SQR(cos_th)*(1-cos(phi)))
                            - u_in(IM3,k,j,i)*cos_th*sin(phi);
            u_k_avg[IM3] +=   u_in(IM1,k,j,i)*sin_th*sin(phi)
                            + u_in(IM2,k,j,i)*cos_th*sin(phi)
                            + u_in(IM3,k,j,i)*cos(phi);
          } // k loop
          // set the new value
          for (int k=ks_avg; k<=ke_avg; ++k) {
            Real phi = - dphi * (k-ks_avg - (n_avg_(j)-1)/2.); // angle between averaged cell center and the current cell
            u_out(IDN,k,j,i) = u_k_avg[IDN]/n_avg_(j);
            u_out(IEN,k,j,i) = u_k_avg[IEN]/n_avg_(j);
            // below physically means u_out = u_k_avg/n_avg_(j)
            u_out(IM1,k,j,i) =(  u_k_avg[IM1]*(1-SQR(sin_th)*(1-cos(phi)))
                               - u_k_avg[IM2]*sin_th*cos_th*(1-cos(phi))
                               - u_k_avg[IM3]*sin_th*sin(phi))/n_avg_(j);
            u_out(IM2,k,j,i) =(- u_k_avg[IM1]*sin_th*cos_th*(1-cos(phi))
                               + u_k_avg[IM2]*(1-SQR(cos_th)*(1-cos(phi)))
                               - u_k_avg[IM3]*cos_th*sin(phi))/n_avg_(j);
            u_out(IM3,k,j,i) =(  u_k_avg[IM1]*sin_th*sin(phi)
                               + u_k_avg[IM2]*cos_th*sin(phi)
                               + u_k_avg[IM3]*cos(phi))/n_avg_(j);
          } // k loop
        } // i loop
      } // l loop
    } // j loop
  }

  // add new average algorithms here (as new modes)

  // if asking for a unsupported mode
  else {
    throw std::invalid_argument( "Invalid polar average mode!" );
  }
  return;
}


void Hydro::Get_block_N_zone_avg(MeshBlock *pmb){

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int j=js; j<=je; ++j) {
    n_avg_(j) = 1;
  }
  
  
  //check if we're on the inner x2-boundary
  if(pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY) {
    for (int j=js; j<=je; ++j) {
      //int n_avg_temp = pmb->block_size.nx3/pow(2,(j-js));
      int n_avg_temp = round(pmb->block_size.nx3/pow(2,round(log2(j-js+1))));
      if( n_avg_temp>1){     
	n_avg_(j) = n_avg_temp;
      }else{
	n_avg_(j) = 1;
      }
    }
  }
  if (pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY)  {
    // average the je index
    for (int j=js; j<=je; ++j) {
      //int n_avg_temp = pmb->block_size.nx3/pow(2,(je-j));
      int n_avg_temp = round(pmb->block_size.nx3/pow(2,round(log2(je-j+1))));
      if( n_avg_temp>1){     
	n_avg_(j) = n_avg_temp;
      }else{
	n_avg_(j) = 1;
      }
    }    
  }
  
}


