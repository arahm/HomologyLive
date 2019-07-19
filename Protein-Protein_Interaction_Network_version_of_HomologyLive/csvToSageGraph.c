// csvToSageGraph: C++ code by Alexander D. Rahm (Version 1.0, November 2018),
// Convert csv tables of protein interactions to SAGE graph.
///////////////////////////////////////////////////////////////////////////

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


void writeGraphInSageFormat(string name, std::vector<std::vector<int>> TheEdges){
//////////////////////////////////////////////////////////////////////////////////////
// Write out the edges in SAGE format.
///////////////////////////////////////////////////////////////////////////////////////

  int edgeNumber = TheEdges.size();
  
  ofstream d1sageFile(name, ios::out);
  d1sageFile << "g = Graph([" << endl;

  for (int n=0; n < edgeNumber; n++){
    d1sageFile << "(" << TheEdges[n][0]+1  << "," << TheEdges[n][1]+1  << ")," <<  endl;
  }
  d1sageFile << "])" << endl;
  d1sageFile << "g.show()" << endl;
  d1sageFile.close();  
  std::cout << "SAGE graph written into the file "  << name << std::endl;    
}




std::vector<std::vector<int>> ConvertEntriesMatrixToEdges(vector<vector<string>> the_entries, int threshold){
//////////////////////////////////////////////////////////////////////////////
// Construct the edge-vertex incidence matrix (d_1) from the csv entries matrix:
//////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> proteinNames;
  std::vector<std::vector<int>> TheEdges;
  std::cout << "Registering only protein-protein interactions stronger than " << threshold << "." <<std::endl;
  for (int i=0; i< the_entries.size(); i++)
  {
    int j=3; 
    int thisInteractionValue = std::stoi(the_entries.at(i).at(j));

    if( thisInteractionValue > threshold){	
	j=1; 
	string thisProteinName = the_entries.at(i).at(j);
	std::vector<std::string>::iterator it;
	it = find (proteinNames.begin(), proteinNames.end(), thisProteinName);
  	if (it != proteinNames.end())  // if this protein in already among the recorded names,
	{ } 
	else {proteinNames.push_back(thisProteinName);};
	it = find (proteinNames.begin(), proteinNames.end(), thisProteinName);
	
	std::vector<int> CurrentEdge(2);
	CurrentEdge[0] = std::distance( proteinNames.begin(), it );

	j=2; 
	thisProteinName = the_entries.at(i).at(j);
	std::vector<std::string>::iterator edge_end_it;
	edge_end_it = find (proteinNames.begin(), proteinNames.end(), thisProteinName);
  	if (edge_end_it != proteinNames.end())  // if this protein in already among the recorded names,
	{ } 
	else {proteinNames.push_back(thisProteinName);};
	edge_end_it = find (proteinNames.begin(), proteinNames.end(), thisProteinName);
	
	CurrentEdge[1] = std::distance( proteinNames.begin(), edge_end_it );

	TheEdges.push_back(CurrentEdge);
    };
  };
  std::cout << "The numbering of the proteins in the graph will be as follows:" <<std::endl;
  for (int i=0; i< proteinNames.size(); i++)
  	std::cout << i+1 << ": " << proteinNames.at(i) << std::endl;
return TheEdges;
}





///////////////////////////////////////////
// Main Program: Extract edges from a csv file.
///////////////////////////////////////////

int main (int argc, char** threshold_and_csvfilename)
{ 
  if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ./csvToSageGraph.o <threshold> <csv_file_with_protein-protein_interaction_data>" << std::endl;
        return -1;
  }
  string csvfilename = threshold_and_csvfilename[2];
  int threshold = std::stoi(threshold_and_csvfilename[1]);
  std::ifstream csvFile(csvfilename);
  if (!csvFile) {
        std::cerr << "Error opening csv file: " << csvfilename << std::endl;
        return -1;
  }
  string row; 
  std::getline(csvFile, row); // discard the first row, which specifies the format of the file.
  vector<vector<string>> the_entries;

  while( std::getline(csvFile, row))
  {
    	std::stringstream row_as_stringstream(row);
    	int i; i = 0;
	vector<string> row_as_vector;

	while( row_as_stringstream.good() )
	{
    		string substr;
    		getline( row_as_stringstream, substr, ',' );
    		row_as_vector.push_back( substr );
	};
    	the_entries.push_back( row_as_vector );
  };
  ////////////////////////////////////////////////////////////
  // Now we assemble the entries to an edges matrix, and write it into a Sage file:
  ////////////////////////////////////////////////////////////
  std::vector<std::vector<int>> TheEdges = ConvertEntriesMatrixToEdges(the_entries, threshold);



  char outputFilename[60]; strcpy(outputFilename, csvfilename.c_str()); strcat(outputFilename, "AtThreshold"); string thrshld = std::to_string(threshold); strcat(outputFilename, thrshld.c_str()); strcat(outputFilename, ".txt");

  writeGraphInSageFormat(outputFilename, TheEdges);

  return 0;
}

