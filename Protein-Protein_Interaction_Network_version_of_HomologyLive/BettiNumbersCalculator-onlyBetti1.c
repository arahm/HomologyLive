///////////////////////////////////////////////////////////////////////////
// BettiNumbersCalculator code of HomologyLive: 
// C++ code by Alexander D. Rahm (based on rationalHomology.c Version 1.0, April 13, 2018),
// for importing the Vietoris-Rips complex truncated to three dimensions,
// and computing its rational homology using the LinBox library for the rank of sparse matrices.
//
// Compile with: g++ -O2 -Wall -g -DNDEBUG -U_LB_DEBUG -DDISABLE_COMMENTATOR -I/usr/local/include -fopenmp -fabi-version=6 -L/usr/local/lib -llinbox -fopenmp -lopenblas -lgivaro -lgmp -lgmpxx -o BettiNumbersCalculator.o BettiNumbersCalculator.c
///////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
// libraries needed specifically for the sparse matrix rank computation: 
#include <givaro/givrational.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/rank.h>
#define SP_STOR SparseMatrixFormat::SparseSeq
using namespace LinBox;

int main (int argc, char **addressfilename)
{ 
  LinBox::Timer tim ; tim.clear() ; tim.start();
  if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ./rationalHomology.o <address-file-linking-to-matrix-files> [<p>]" << std::endl;
        return -1;
  }
  std::ifstream addressFile (addressfilename[1]);
  if (!addressFile) {
        std::cerr << "Error opening address file: " << addressfilename[1] << std::endl;
        return -1;
  }
  string threshold; 
  string d1filename; 
  string d2filename; 
 /* string d3filename; */
  std::getline(addressFile, threshold);   
  std::getline(addressFile, d1filename);   
  std::getline(addressFile, d2filename);   
 /* std::getline(addressFile, d3filename);        */ 
  std::cout << "At threshold " << threshold << ", the Vietoris-Rips complex truncated to three dimensions has the following cell numbers and topological invariants:" << std::endl; 

  ////////////////////////////////////////////////////////////
  //Compute the rational homology of the Vietoris-Rips complex
  ////////////////////////////////////////////////////////////
    Givaro::QField<Givaro::Rational> ZZ;
    long unsigned int r1;
    std::ifstream input (d1filename);
    if (!input) {
        std::cerr << "Error opening d1 matrix file: " << std::endl;
    }
    MatrixStream<Givaro::QField<Givaro::Rational>> m1s( ZZ, input );
    SparseMatrix<Givaro::QField<Givaro::Rational>, SP_STOR> d1 ( m1s );
    // Compute the rank of the sparse matrix d1 over the rational numbers:
      if (d1.size() > 0){
        LinBox::rank (r1, d1);
      }
      else{ r1 = 0;}

    long unsigned int r2;
    std::ifstream secondinput (d2filename);
    if (!secondinput) {
        std::cerr << "Error opening d2 matrix file: " << std::endl;
    }
    MatrixStream<Givaro::QField<Givaro::Rational>> m2s( ZZ, secondinput );
    SparseMatrix<Givaro::QField<Givaro::Rational>, SP_STOR> d2( m2s );
    // Compute the rank of the sparse matrix d2 over the rational numbers:
      if (d2.size() > 0){
        LinBox::rank (r2, d2);
      }
      else{ r2 = 0;}
/*
    long unsigned int r3;
    std::ifstream thirdinput (d3filename);
    if (!thirdinput) {
        std::cerr << "Error opening d3 matrix file: " << std::endl;
    }
    MatrixStream<Givaro::QField<Givaro::Rational>> m3s( ZZ, thirdinput );
    SparseMatrix<Givaro::QField<Givaro::Rational>, SP_STOR> d3( m3s );
    // Compute the rank of the sparse matrix d3 over the rational numbers:
      if (d3.size() > 0){
        LinBox::rank (r3, d3);
      }
      else{ r3 = 0;}
 */
            std::cout << "There are " /*<< d3.rowdim() << " tetrahedra, " */<< d2.rowdim() << " triangles, " << d1.rowdim() << " edges and " << d1.coldim() << " vertices." << std::endl; 
std::cout << std::endl; 
/* std::cout << " Two-dimensional bubbles: " << d2.rowdim() -r2 -r3 << std::endl; // 2nd Betti number */
std::cout << " One-dimensional loops: " << d1.rowdim() -r1 -r2 << std::endl; // 1st Betti number
std::cout << " Connected components (\"clusters\"): " << d1.coldim() -r1 << std::endl; // 0th Betti number         
    tim.stop();
    std::cout << "Computing the ranks of the differential matrices took " << tim << " time for this Vietoris-Rips complex." << std::endl;

ofstream resultsfile("resultsOnFirstBettiNumber.txt", std::ios_base::app | std::ios_base::out);
  resultsfile   << "At threshold " << threshold << ", the Vietoris-Rips complex truncated to three dimensions has the following cell numbers and topological invariants:" << std::endl; 
  resultsfile<< "There are " /*<< d3.rowdim() << " tetrahedra, " */<< d2.rowdim() << " triangles, " << d1.rowdim() << " edges and " << d1.coldim() << " vertices." << std::endl; 
resultsfile<< std::endl; 
/* resultsfile<< " Two-dimensional bubbles: " << d2.rowdim() -r2 -r3 << std::endl; // 2nd Betti number */
resultsfile<< " One-dimensional loops: " << d1.rowdim() -r1 -r2 << std::endl; // 1st Betti number
resultsfile<< " Connected components (\"clusters\"): " << d1.coldim() -r1 << std::endl; // 0th Betti number         
    resultsfile<< "Computing the ranks of the differential matrices took " << tim << " time for this Vietoris-Rips complex." << std::endl;
  resultsfile.close();    
    return 0;
}

