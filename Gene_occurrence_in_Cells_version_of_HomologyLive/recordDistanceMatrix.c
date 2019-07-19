/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Preparing DGE input for HomologyLive: C++ code by Alexander D. Rahm (January 15, 2019),
// subject to a GNU General Public License.
// Read in a csv table with biological cells as columns and their genomes as rows.
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
  std::cout << "Converting the input DGE matrix into a distance matrix ..." <<std::endl;
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




///////////////////////////////////////////
// Main Program: Read DGE data from a csv file, and record distance matrix.
///////////////////////////////////////////

int main (int argc, char** input_csvfilename)
{ 
  if (argc < 2 || argc > 2) {
        std::cerr << "Usage: ./recordDistanceMatrix.o <csv_file_with_DGE_data>" << std::endl;
        return -1;
  }
  string csvfilename = input_csvfilename[1];
  std::ifstream csvFile(csvfilename);
  if (!csvFile) {
        std::cerr << "Error opening csv file: " << csvfilename << std::endl;
        return -1;
  }
  std::cout << "Opened input file: " << csvfilename << std::endl;
  string row; 
  vector<vector<float>> the_entries;


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
                float theEntry = std::stof(substr, &sz);
    		row_as_vector.push_back( theEntry );
	};
    	the_entries.push_back( row_as_vector );
  };
  // Now we have read in the input file as a matrix of strings. 
  // We convert it into the distance matrix:
  std::vector<std::vector<float>> distanceMatrix = ConvertEntriesMatrixToDistanceMatrix(the_entries);
 
  int N_r = distanceMatrix.size(); // Number of biological cells
  
/////////////////////////////////////////////////////////////////////////////////////////
// Write out the distance matrix. 
/////////////////////////////////////////////////////////////////////////////////////////
 char distfilename[60]; strcpy(distfilename, csvfilename.c_str()); strcat(distfilename, "_distances.txt");
  ofstream distancesFile(distfilename, ios::out);
  distancesFile << "#Genetic distances between cells recorded in the following matrix: " << endl;

  // Now we can write out the distance matrix and the thresholds.
  std::list<string> thresholds;   
  for (int i=0; i < N_r; i++){
     for (int j = 0; j < N_r; j++){
	  distancesFile << distanceMatrix[i][j] << ";";
          thresholds.push_back(to_string(distanceMatrix[i][j]));
     };
     distancesFile << endl;
  };
  distancesFile.close();
  thresholds.sort();
  thresholds.unique();
  cout << "Recorded distance matrix of size " << N_r << "x" << N_r << " with " << thresholds.size() << " critical values." << endl;
 
///////////////////////////////////////////////////////
// Write a shell script pointing to the distance matrix and threholds,  
////////////////////////////////////////////////////////
 
 char addressfilename[60]; strcpy(addressfilename, csvfilename.c_str()); strcat(addressfilename, "_threshold_variation.sh");
  ofstream shellFile(addressfilename, std::ios_base::app | std::ios_base::out);
  for (int j = 0; j < thresholds.size(); j++){ 
        shellFile << "./HomologyLiveDGE.o " << thresholds.front() << " " << distfilename << endl;
        thresholds.pop_front();
  };
  shellFile.close();

  return 0;
}

