/////////////////////////////////////////////////////////////////////////////////////////////////////////
// HomologyLive with csv input and LinBox output: C++ code by Alexander D. Rahm (January 15, 2019),
// subject to a GNU General Public License.
// Read in csv tables of protein interactions, and construct the Vietoris-Rips complex on their graph, 
// truncated to three dimensions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <list>
#include <vector>
#include <tuple>
#include <algorithm>
using namespace std;


std::vector<std::vector<float>> ConvertEntriesMatrixToDistanceMatrix(vector<vector<float>> the_entries){
////////////////////////////////////////////////////////////////////////////////
// Construct the distance matrix from the input matrix of genome multiplicities:
////////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<float>> distanceMatrix;
  std::vector<std::vector<float>> distanceTriangle; // We only need to compute the upper triangle
  // of the distance matrix, because the latter is symmetric.
  std::cout << "Converting the input DGE matrix into a distance matrix:" <<std::endl;
  int numberOfRows = the_entries.size();
  int numberOfColumns = the_entries.at(0).size();
  for (int i=0; i < numberOfColumns; i++)
  {
    vector<float> ithDistanceRow;    
    for (int j = i+1; j < numberOfColumns; j++)     
    {
        float thisDifference = 0;
        for (int k = 0; k < numberOfRows; k++)
        {
                thisDifference = thisDifference + std::abs(the_entries.at(k).at(i) -the_entries.at(k).at(j));
        }
        ithDistanceRow.push_back(thisDifference); 
    };
    distanceTriangle.push_back(ithDistanceRow);
  };
  // Enter the recorded values into the upper triangle of the distance matrix:
  for (int i=0; i < numberOfColumns; i++)
  {
    vector<float> full_ithDistanceRow(numberOfColumns);   
    for (int j = 0; j < distanceTriangle.at(i).size(); j++)     
    {  
        full_ithDistanceRow[i+1+j] = distanceTriangle.at(i).at(j);
    };
    distanceMatrix.push_back(full_ithDistanceRow);
   };
  // Mirror the upper triangle into the lower triangle to get the distance matrix:
  for (int i=0; i < numberOfColumns; i++)
  {
    for (int j = 0; j < i; j++)     
    {  
        distanceMatrix[i][j] = distanceMatrix[j][i];
    };
   };
//control output of the distance matrix:
//  for (int i=0; i < numberOfColumns; i++)
//  {
//    for (int j = 0; j < numberOfColumns; j++)     
//    {  
//       cout << distanceMatrix.at(i).at(j) << " , ";
//    };
//    cout << endl;
//   };
return distanceMatrix;
}


std::tuple<std::vector<std::vector<int>>, int>  ConstructEdgeVertexIncidences(std::vector<std::list<int>> adjacencies, int N_r){
//////////////////////////////////////////////////////////////////////////////
// Construct the edge-vertex incidence matrix (d_1) from the adjacency matrix:
//////////////////////////////////////////////////////////////////////////////
  int edgeNumber = 0; 
  std::vector<std::vector<int>> EdgeVertexIncidence(N_r);
  for (int i=0; i < N_r; i++){
    std::list<int> ithAdjacency = adjacencies[i];
    if(0 < ithAdjacency.size()){
      std::vector<int> EVIi(N_r);
      while (0 < ithAdjacency.size() ){
	  edgeNumber++;
	  int j = ithAdjacency.front();
	  EVIi[j] = edgeNumber;
	  // outFile << edgeNumber << "\t" << i+1 << "\t" << 1  <<  endl;
	  // outFile << edgeNumber << "\t" << j+1 << "\t" << -1 <<  endl;
	  ithAdjacency.pop_front();
      }
      EdgeVertexIncidence[i] = EVIi;
    }
  }
// std::cout << "EdgeVertexIncidence:\n";
/// control output, please deactivate this "for" loop:
//  for (int m=0; m<N_r; m++){
//    for (int n=0; n<EdgeVertexIncidence[m].size(); n++){
//        std::cout << EdgeVertexIncidence[m][n] << ",";
//     }	
//     std::cout << "\n";
//  }  
return std::make_tuple(EdgeVertexIncidence, edgeNumber);
}


std::tuple<std::vector<std::vector<std::vector<int>>>, std::list<std::vector<int>>, std::list<std::vector<int>>, int>  ConstructTriangleVertexIncidences(std::vector<std::list<int>> adjacencies, std::vector<std::vector<int>> EdgeVertexIncidence, int N_r){
////////////////////////////////////////////////////////////////////////////////////
// Construct the triangles-edges incidence matrix (d_2) from the adjacency matrix
// and the edge-vertex incidence matrix:
////////////////////////////////////////////////////////////////////////////////////
  int   triangleNumber = 0;
  std::vector<std::vector<std::vector<int>>> TriangleVertexIncidence(N_r); 
  std::list<std::vector<int>> Triangles;
  std::list<std::vector<int>> TriangleBoundaries;

  for (int i=0; i < N_r; i++){
      if(adjacencies[i].size() > 1){	

	// At vertex number i, abbreviate TriangleVertexIncidence by TVI:
	std::vector<std::vector<int>> TVIi(N_r); 

	// second vertex j:
	std::list<int> secondVertices = adjacencies[i];
	while(secondVertices.size() > 1){
       	 int j = secondVertices.front();
	 std::vector<int> TVIij(N_r); 

	 // Third vertex k:
	 std::list<int> thirdVertices = secondVertices; 
	 thirdVertices.pop_front();
	 while(thirdVertices.size() > 0){
       	  int k = thirdVertices.front();
       	  if(i < j && j < k){
	   // Check that vertex k is adjacent to vertex j:
	   if (EdgeVertexIncidence[j].size() > 0){
	    if (EdgeVertexIncidence[j][k] > 0){

	  	triangleNumber++;
		TVIij[k] = triangleNumber;

		std::vector<int> TriangleIndices(3);
		TriangleIndices[0] = i;
		TriangleIndices[1] = j;
		TriangleIndices[2] = k;
		Triangles.push_back(TriangleIndices);

		std::vector<int> TriangleBoundary(3);
		TriangleBoundary[0] = EdgeVertexIncidence[i][j];
		TriangleBoundary[1] = EdgeVertexIncidence[j][k];
		TriangleBoundary[2] = EdgeVertexIncidence[i][k];
		TriangleBoundaries.push_back(TriangleBoundary);
	     }
	    }
	   }
	   thirdVertices.pop_front();
	  }
	secondVertices.pop_front();
	TVIi[j] = TVIij;
	}
	TriangleVertexIncidence[i] = TVIi;
      }
  }
return std::make_tuple(TriangleVertexIncidence, Triangles, TriangleBoundaries, triangleNumber);
}



std::list<std::vector<int>>  ConstructTetrahedraTrianglesIncidences(std::vector<std::list<int>> adjacencies,  std::vector<std::vector<std::vector<int>>> TriangleVertexIncidence, std::list<std::vector<int>> Triangles){
////////////////////////////////////////////////////////////////////////////////////
// Construct the tetrahedra-triangles incidence matrix (d_3) 
// from the adjacency matrix, the Triangles list
// and the TriangleVertexIncidence sparse cubus matrix:
////////////////////////////////////////////////////////////////////////////////////
  int   tetrahedraNumber = 0;
  std::list<std::vector<int>> TetrahedraBoundaries;

  while (Triangles.size() > 0){

      std::vector<int> ijk = Triangles.front();
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];
      if(!(i < j && j < k)){cerr << "Vertex disorder \n";}
      else{
	// Potential fourth vertex l:
	std::list<int> potentialFourthVertices = adjacencies[k];
	while(potentialFourthVertices.size() > 0){
       	  int l = potentialFourthVertices.front();
	  if( l > k ){

	   // Check that vertex l forms a triangle with each edge of the triangle (i,j,k):
	   if (TriangleVertexIncidence[i].size() > 0 && TriangleVertexIncidence[j].size() > 0){
	    if (TriangleVertexIncidence[i][j].size() > 0 && TriangleVertexIncidence[j][k].size() > 0 && TriangleVertexIncidence[i][k].size() > 0){
	     if (TriangleVertexIncidence[i][j][l] > 0 && TriangleVertexIncidence[j][k][l] > 0 && TriangleVertexIncidence[i][k][l] > 0){

	  	tetrahedraNumber++;

		std::vector<int> TetrahedronBoundary(4);
		TetrahedronBoundary[0] = TriangleVertexIncidence[i][j][l];
		TetrahedronBoundary[1] = TriangleVertexIncidence[j][k][l];
		TetrahedronBoundary[2] = TriangleVertexIncidence[i][k][l];
		TetrahedronBoundary[3] = TriangleVertexIncidence[i][j][k];
		TetrahedraBoundaries.push_back(TetrahedronBoundary);
	      }
	     }
	    }
	  }
	potentialFourthVertices.pop_front();
	}
      }
      Triangles.pop_front();
  }
return TetrahedraBoundaries;
}




///////////////////////////////////////////
// Main Program: Extract edges from a csv file.
///////////////////////////////////////////

int main (int argc, char** threshold_and_distfilename)
{ 
  if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ./HomologyLiveDGE.o <threshold> <txt_file_with_distance_matrix>" << std::endl;
        return -1;
  }
  string distfilename = threshold_and_distfilename[2];
  std::string::size_type szend;
  float threshold = std::stof(threshold_and_distfilename[1], &szend);
  std::ifstream csvFile(distfilename);
  if (!csvFile) {
        std::cerr << "Error opening distance matrix file: " << distfilename << std::endl;
        return -1;
  }
  std::cout << "Opened input file: " << distfilename << std::endl;
  string row; 
  vector<vector<float>> the_entries;
  std::getline(csvFile, row); // discard the first row, which roughly describes the contents of the file.

  while( std::getline(csvFile, row))
  {
    	std::stringstream row_as_stringstream(row);
    	int i; i = 0;
	vector<float> row_as_vector;

	while( row_as_stringstream.good() )
	{
    		string substr;
    		getline( row_as_stringstream, substr, ';' );
                std::string::size_type sz;
                if(!row_as_stringstream.eof()){
                        float theEntry = std::stof(substr, &sz);
    		        row_as_vector.push_back( theEntry );
                };
	};
    	the_entries.push_back( row_as_vector );
  };
  // Now we have read in the input file as a matrix of strings. 
  // We convert it into the distance matrix:
  std::vector<std::vector<float>> distanceMatrix = the_entries;

/////////////////////////////////////////////////////////////////////////
// Construct the vertex-vertex adjacency matrix from the distance matrix.
/////////////////////////////////////////////////////////////////////////
  int N_r = distanceMatrix.size();
  std::vector<std::list<int>> adjacencies(N_r);
  int edgeNumber = 0; 
  for (int i=0; i<N_r; i++){
      std::list<int> ithadjacency; 
      for (int j=i+1; j<N_r; j++){
	if (distanceMatrix[i][j] < threshold ){
	  	edgeNumber++;
		ithadjacency.push_back(j);
        }
      }
      ithadjacency.sort();
      adjacencies[i] = ithadjacency;
  };
  // Construct the edge-vertex incidence matrix from the adjacency matrix:  
 std::tuple<std::vector<std::vector<int>>, int> returnedData  = ConstructEdgeVertexIncidences(adjacencies, N_r);
 std::vector<std::vector<int>> EdgeVertexIncidence = std::get<0>(returnedData);
 if (edgeNumber != std::get<1>(returnedData) ){cerr << "Edge counting error \n";}

  ///////////////////////////
  // Construct the triangles 
  ///////////////////////////
 std::tuple<std::vector<std::vector<std::vector<int>>>, std::list<std::vector<int>>, std::list<std::vector<int>>, int> d2data = ConstructTriangleVertexIncidences(adjacencies, EdgeVertexIncidence, N_r);
  std::vector<std::vector<std::vector<int>>> TriangleVertexIncidence = std::get<0>(d2data);  
  std::list< std::vector<int> > Triangles = std::get<1>(d2data);
  std::list< std::vector<int> > TriangleBoundaries = std::get<2>(d2data);
  int   triangleNumber = std::get<3>(d2data);
  ///////////////////////////
  // Construct the tetrahedra
  ///////////////////////////
  std::list<std::vector<int>> TetrahedraBoundaries = ConstructTetrahedraTrianglesIncidences(adjacencies, TriangleVertexIncidence, Triangles);
  int   tetrahedraNumber = TetrahedraBoundaries.size();



/////////////////////////////////////////////////////////////////////////////////////////
// Write out the edge-vertex incidence matrix in LinBox format. 
/////////////////////////////////////////////////////////////////////////////////////////
 char d1filename[60]; strcpy(d1filename, distfilename.c_str()); strcat(d1filename, "_d1LinBox"); string thrshld = std::to_string(threshold); strcat(d1filename, thrshld.c_str()); strcat(d1filename, ".txt");
  ofstream d1outFile(d1filename, ios::out);
  d1outFile << "#Edge-vertex incidence matrix of the protein-protein interaction network in Linbox format: " << endl;

  // Now we can write out the edge-vertex incidence matrix, which is the d1-boundary map.
  d1outFile << edgeNumber << "\t" << N_r << "\t" << "M" <<  endl;   
for (int i=0; i < N_r; i++){
    if(0 < EdgeVertexIncidence[i].size()){
	for (int j = 0; j < N_r; j++){
	 if(EdgeVertexIncidence[i][j] > 0){
	  d1outFile << EdgeVertexIncidence[i][j] << "\t" << i+1 << "\t" << "-1"  <<  endl;
	  d1outFile << EdgeVertexIncidence[i][j] << "\t" << j+1 << "\t" << "1"  <<  endl;
	 }
        }
    }
  }
  d1outFile << 0 << "\t" << 0 << "\t" << 0 <<  endl;
  d1outFile.close();

//////////////////////////////////////////////////////////////////////////////////////////////
// Write out the triangles-edges incidence matrix (d_2) in LinBox format:
//////////////////////////////////////////////////////////////////////////////////////////////
  char d2filename[60]; strcpy(d2filename, distfilename.c_str()); strcat(d2filename, "_d2LinBox");  strcat(d2filename, thrshld.c_str()); strcat(d2filename, ".txt");
  ofstream d2outFile(d2filename, ios::out);

  d2outFile << "#Triangles-Edges incidence matrix formed by the protein-protein interaction network: " << endl;

  d2outFile << triangleNumber << "\t" << edgeNumber << "\t" << "M" <<  endl;
  triangleNumber = 0;
  std::list< std::vector<int> > CheckedTriangleBoundaries = TriangleBoundaries;
  while (CheckedTriangleBoundaries.size() > 0){
	  triangleNumber++;
	  int ij = CheckedTriangleBoundaries.front()[0];
	  int jk = CheckedTriangleBoundaries.front()[1];
	  int ik = CheckedTriangleBoundaries.front()[2];
	  		
	  d2outFile << triangleNumber << "\t" << ij << "\t" << 1  <<  endl;
	  d2outFile << triangleNumber << "\t" << jk << "\t" << 1  <<  endl;
	  d2outFile << triangleNumber << "\t" << ik << "\t" <<-1  <<  endl;
	  
	  CheckedTriangleBoundaries.pop_front();
  }
  d2outFile << 0 << "\t" << 0 << "\t" << 0 <<  endl;
  d2outFile.close();

/////////////////////////////////////////////////////////////////////////////////////////////////
// Write out the tetrahedra-triangles incidence matrix (d_3) in LinBox format:
/////////////////////////////////////////////////////////////////////////////////////////////////

 char d3filename[60]; strcpy(d3filename, distfilename.c_str()); strcat(d3filename, "_d3LinBox");  strcat(d3filename, thrshld.c_str()); strcat(d3filename, ".txt");
  ofstream d3outFile(d3filename, ios::out);

  d3outFile << "#Tetrahedra-triangles incidence matrix:" << endl;

  d3outFile << tetrahedraNumber << "\t" << triangleNumber << "\t" << "M" <<  endl;
  tetrahedraNumber = 0;
  std::list< std::vector<int> > CheckedTetrahedraBoundaries = TetrahedraBoundaries;
  while (CheckedTetrahedraBoundaries.size() > 0){
	  tetrahedraNumber++;
	  d3outFile << tetrahedraNumber << "\t" << CheckedTetrahedraBoundaries.front()[0] << "\t" << 1  <<  endl;
	  d3outFile << tetrahedraNumber << "\t" << CheckedTetrahedraBoundaries.front()[1] << "\t" << 1  <<  endl;
	  d3outFile << tetrahedraNumber << "\t" << CheckedTetrahedraBoundaries.front()[2] << "\t" <<-1  <<  endl;
	  d3outFile << tetrahedraNumber << "\t" << CheckedTetrahedraBoundaries.front()[3] << "\t" <<-1  <<  endl;
	  CheckedTetrahedraBoundaries.pop_front();
  }
  d3outFile << 0 << "\t" << 0 << "\t" << 0 <<  endl;
  d3outFile.close();
 
///////////////////////////////////////////////////////
// Write an address file pointing to the matrix files,
/// and insert it into a shell script for being treated.   
////////////////////////////////////////////////////////
 
 char addressfilename[60]; strcpy(addressfilename, distfilename.c_str()); strcat(addressfilename, "_address");  strcat(addressfilename, thrshld.c_str()); strcat(addressfilename, ".txt");

  ofstream addressFile(addressfilename, ios::out);
  addressFile << thrshld.c_str() << endl;
  addressFile << d1filename <<  endl;
  addressFile << d2filename <<  endl;
  addressFile << d3filename <<  endl;
  addressFile.close();

  ofstream shellFile("getBettiNumbers.sh", std::ios_base::app | std::ios_base::out);
  shellFile << "./BettiNumbersCalculator.o " << addressfilename << endl;
  shellFile.close();

//////////////////////////////////////////////////////////////////////////////////////
// Write out the differential matrices in SAGE format.
///////////////////////////////////////////////////////////////////////////////////////
/* char d1sagefilename[60]; strcpy(d1sagefilename, distfilename.c_str()); strcat(d1sagefilename, "_d1_threshold");  strcat(d1sagefilename, thrshld.c_str()); strcat(d1sagefilename, ".sage");
  ofstream d1sageFile(d1sagefilename, ios::out);
  d1sageFile << "d1 = matrix(GF(997)," << edgeNumber << "," << N_r << ", sparse=True);\n" << endl;

  // Now we can write out the edge-vertex incidence matrix, which is the d1-boundary map.  
  for (int i=0; i < N_r; i++){
    if(0 < EdgeVertexIncidence[i].size()){
	for (int j = 0; j < N_r; j++){
	 if(EdgeVertexIncidence[i][j] > 0){
	  d1sageFile << "d1[" << EdgeVertexIncidence[i][j]-1 << "," << i << "] = -1;\n"  <<  endl;
	  d1sageFile << "d1[" << EdgeVertexIncidence[i][j]-1 << "," << j << "] = +1;\n"  <<  endl;
	 }
        }
    }
  }
  d1sageFile.close();
*/
//////////////////////////////////////////////////////////////////////////////////////////////
// Write out the triangles-edges incidence matrix (d_2) into a sage matrix:
//////////////////////////////////////////////////////////////////////////////////////////////
/*  std::list< std::vector<int> > sageCheckTriangleBoundaries = TriangleBoundaries;
 char d2sagefilename[60]; strcpy(d2sagefilename, distfilename.c_str()); strcat(d2sagefilename, "_d2_threshold");  strcat(d2sagefilename, thrshld.c_str()); strcat(d2sagefilename, ".sage");
  ofstream d2sageFile(d2sagefilename, ios::out);

  d2sageFile << "d2 = matrix(GF(997)," << triangleNumber << "," << edgeNumber << ", sparse=True)\n" << endl;
  
  triangleNumber = 0;
//  std::cout  << " triangleNumber = " << triangleNumber <<endl; 
  while (sageCheckTriangleBoundaries.size() > 0){
	  triangleNumber++;

	  int ij = sageCheckTriangleBoundaries.front()[0];
	  int jk = sageCheckTriangleBoundaries.front()[1];
	  int ik = sageCheckTriangleBoundaries.front()[2];
	  d2sageFile << "d2[" << triangleNumber-1 << "," << ij-1 << "] = 1;\n"  <<  endl;
	  d2sageFile << "d2[" << triangleNumber-1 << "," << jk-1 << "] = 1;\n"  <<  endl;
	  d2sageFile << "d2[" << triangleNumber-1 << "," << ik-1 << "] = -1;\n" <<  endl;	  
	  sageCheckTriangleBoundaries.pop_front();
  }
  d2sageFile.close();
*/
/////////////////////////////////////////////////////////////////////////////////////////////////
// Write out the tetrahedra-triangles incidence matrix (d_3) we obtain from the adjacency matrix:
/////////////////////////////////////////////////////////////////////////////////////////////////
/*  std::list< std::vector<int> > sageCheckTetrahedraBoundaries = TetrahedraBoundaries;
 char d3sagefilename[60]; strcpy(d3sagefilename, distfilename.c_str()); strcat(d3sagefilename, "_d3_threshold");  strcat(d3sagefilename, thrshld.c_str()); strcat(d3sagefilename, ".sage");
  ofstream d3sageFile(d3sagefilename, ios::out);

  d3sageFile << "d3 = matrix(GF(997)," << tetrahedraNumber << "," << triangleNumber << ", sparse=True);\n" << endl;

  tetrahedraNumber = 0;
  while (sageCheckTetrahedraBoundaries.size() > 0){
	  tetrahedraNumber++;
	  d3sageFile << "d3[" << tetrahedraNumber-1 << "," << sageCheckTetrahedraBoundaries.front()[0]-1 << "] = 1; "  <<  endl;
	  d3sageFile << "d3[" << tetrahedraNumber-1 << "," << sageCheckTetrahedraBoundaries.front()[1]-1 << "] = 1; "  <<  endl;
	  d3sageFile << "d3[" << tetrahedraNumber-1 << "," << sageCheckTetrahedraBoundaries.front()[2]-1 << "] =-1; "  <<  endl;
	  d3sageFile << "d3[" << tetrahedraNumber-1 << "," << sageCheckTetrahedraBoundaries.front()[3]-1 << "] =-1; "  <<  endl;
	  sageCheckTetrahedraBoundaries.pop_front();
  }
  d3sageFile.close();
  
   std::cout  << "Written the first three boundary matrices of the Vietoris-Rips complex into the files d1.sage, d2.sage and d3.sage." <<endl; 

 */ 
/*  //////////////////////////////////////////////////////////////////////////////////////
// Write out the differential matrices in PARI/GP format for a correctness check.
// This can be skipped by de-actvating the remaining source code in this procedure.
///////////////////////////////////////////////////////////////////////////////////////
  ofstream d1pariFile("d1pari.gp", ios::out);
  d1pariFile << "d1 = matrix(" << edgeNumber << "," << N_r << ", X, Y, 0);\n" << endl;

  // Now we can write out the edge-vertex incidence matrix, which is the d1-boundary map.  
  for (int i=0; i < N_r; i++){
    if(0 < EdgeVertexIncidence[i].size()){
	for (int j = 0; j < N_r; j++){
	 if(EdgeVertexIncidence[i][j] > 0){
	  d1pariFile << "d1[" << EdgeVertexIncidence[i][j] << "," << i+1 << "] = -1;\n"  <<  endl;
	  d1pariFile << "d1[" << EdgeVertexIncidence[i][j] << "," << j+1 << "] = +1;\n"  <<  endl;
	 }
        }
    }
  }
  d1pariFile.close();

//////////////////////////////////////////////////////////////////////////////////////////////
// Write out the triangles-edges incidence matrix (d_2) into a PARI/GP matrix:
//////////////////////////////////////////////////////////////////////////////////////////////
  std::list< std::vector<int> > PariCheckTriangleBoundaries = TriangleBoundaries;
  ofstream d2pariFile("d2pari.gp", ios::out);

  d2pariFile << "d2 = matrix(" << triangleNumber << "," << edgeNumber << ", X, Y, 0);\n" << endl;
  
  triangleNumber = 0;
//  std::cout  << " triangleNumber = " << triangleNumber <<endl; 
  while (PariCheckTriangleBoundaries.size() > 0){
	  triangleNumber++;

	  int ij = PariCheckTriangleBoundaries.front()[0];
	  int jk = PariCheckTriangleBoundaries.front()[1];
	  int ik = PariCheckTriangleBoundaries.front()[2];
	  d2pariFile << "d2[" << triangleNumber << "," << ij << "] = 1;\n"  <<  endl;
	  d2pariFile << "d2[" << triangleNumber << "," << jk << "] = 1;\n"  <<  endl;
	  d2pariFile << "d2[" << triangleNumber << "," << ik << "] = -1;\n" <<  endl; 
	  PariCheckTriangleBoundaries.pop_front();
  }
  d2pariFile.close();

/////////////////////////////////////////////////////////////////////////////////////////////////
// Write out the tetrahedra-triangles incidence matrix (d_3) we obtain from the adjacency matrix:
/////////////////////////////////////////////////////////////////////////////////////////////////
  std::list< std::vector<int> > PariCheckTetrahedraBoundaries = TetrahedraBoundaries;
  ofstream d3pariFile("d3pari.gp", ios::out);

  d3pariFile << "d3 = matrix(" << tetrahedraNumber << "," << triangleNumber << ", X, Y, 0);\n" << endl;

  tetrahedraNumber = 0;
  while (PariCheckTetrahedraBoundaries.size() > 0){
	  tetrahedraNumber++;
	  d3pariFile << "d3[" << tetrahedraNumber << "," << PariCheckTetrahedraBoundaries.front()[0] << "] = 1; "  <<  endl;
	  d3pariFile << "d3[" << tetrahedraNumber << "," << PariCheckTetrahedraBoundaries.front()[1] << "] = 1; "  <<  endl;
	  d3pariFile << "d3[" << tetrahedraNumber << "," << PariCheckTetrahedraBoundaries.front()[2] << "] =-1; "  <<  endl;
	  d3pariFile << "d3[" << tetrahedraNumber << "," << PariCheckTetrahedraBoundaries.front()[3] << "] =-1; "  <<  endl;
	  PariCheckTetrahedraBoundaries.pop_front();
  }
  d3pariFile.close();
  
   std::cout  << "Written the first three boundary matrices of the Vietoris-Rips complex into the files d1pari.gp, d2pari.gp and d3pari.gp." <<endl; 
 */ 
  ////////////////////////////////////////////////////////////
  // Write the edges matrix into a Sage file:
  ////////////////////////////////////////////////////////////
 // char outputFilename[60]; strcpy(outputFilename, distfilename.c_str()); strcat(outputFilename, "AtThreshold"); /* string thrshld = std::to_string(threshold); */ strcat(outputFilename, thrshld.c_str()); strcat(outputFilename, ".txt");

  return 0;
}

