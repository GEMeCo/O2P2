// ================================================================================================
//
// This is the main file. A single object is created, FEAnalysis. It manages the entire process.
// 
// ================================================================================================

// Custom Header Files
#include "Common.h"
#include "FEAnalysis.h"

/**
  * @brief This main routine creates a single object, FEAnalysis. It manages the entire process.
  * @param argc The number of arguments in the call
  * @param args An array of arguments
  * @return Returns 0 if succeeded
  */
int main(int argc, char** args)
{
	std::string mi_stProj;		// Project
	std::string mi_stFile;		// Input File
	std::ifstream mi_File;		// Input File stream

	std::cout << INFO("Copyright(C) 2024 GEMeCO - All Rights Reserved.") << std::endl
		<< INFO("Structural Engineering Department.") << std::endl
		<< INFO("University of Sao Paulo at Sao Carlos School of Engineering.") << std::endl << std::endl
		<< WARN("This program is free: you can redistribute it under the terms of the License.")
		<< WARN("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;")
		<< WARN("Without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.")
		<< WARN("It is provided \"AS IS\".") << std::endl
		<< INFO("In no event shall the authors be liable for any claim, damages or other liability,") << std::endl
		<< INFO("whether in an action of contract, tort or otherwire, arising from, out of or in ") << std::endl
		<< INFO("connection with the software or the use of other dealing in the software.") << std::endl << std::endl
		<< INFO("Neither the name of the copyright holder nor the names of any other contributors") << std::endl
		<< INFO("may be used to endorse or promote products derived from this software without") << std::endl
		<< INFO("specific prior written permission.") << std::endl << std::endl
		<< INFO("If you are using, please cite our software in your reseach. We have a DOI: 10.5281/zenodo.8283439") << std::endl << std::endl;

	// --------------------------------------------------------------------------------------------
	// Check if program call has the project name
	if (argc > 1)
	{
		std::cout << WARN("The following input project was submitted: ") << args[1];
		mi_stProj = args[1];
	}
	else
	{
		std::cout << INPUT("Notice: Type exit if closing.\n");
		std::cout << INPUT("Project name: ");
		std::cin >> mi_stProj;

		if (mi_stProj.compare("exit") == 0) return 0;
	}
	mi_stFile = mi_stProj + ".txt";

	try
	{
#ifdef _DEBUG
		Instrumentor::beginSession("Profile");

		logFile.open(mi_stProj + ".log", std::ios::trunc);
		LOG("Project Name: " + mi_stProj);
#endif //_DEBUG

		// Any instruction inside a try structure may throw and exception
		mi_File.open(mi_stFile);
		if (!mi_File) {
			LOG("\n\n\nmain: An exception was thrown when opening file " + mi_stFile + "\n\n\n");
			throw std::invalid_argument("\n\n\nAn exception was thrown when opening file " + mi_stFile + "\n\n\n");
		}

		// Read line by line until find flag
		std::string mi_stLine;
		std::string mi_stFlag = "#DIM#";

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(mi_File, mi_stLine);
			if (mi_File.eof()) {
				LOG("\n\n\nmain: Reading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		int mi_nDim;
		mi_File >> mi_nDim;

		if (mi_nDim == 2)
		{
			{
				Timer timer("FEAnalysis.initComponents");
				O2P2::FEAnalysis<2>::initComponents(mi_File);
				mi_File.close();
			}
			{
				Timer timer("FEAnalysis.runAnalysis");
				O2P2::FEAnalysis<2>::runAnalysis();
			}
			{
				Timer timer("FEAnalysis.drawResults");
				O2P2::FEAnalysis<2>::drawResults(mi_stProj);
			}
		}
		else if (mi_nDim == 3)
		{
			{
				Timer timer("FEAnalysis.initComponents");
				O2P2::FEAnalysis<3>::initComponents(mi_File);
				mi_File.close();
			}
			{
				Timer timer("FEAnalysis.runAnalysis");
				O2P2::FEAnalysis<3>::runAnalysis();
			}
			{
				Timer timer("FEAnalysis.drawResults");
				O2P2::FEAnalysis<3>::drawResults(mi_stProj);
			}
		}
		else
		{
			LOG("\n\n\nmain: Wrong dimensionality\nCheck problem input file\n\n\n");
			throw std::invalid_argument("\n\nWrong dimensionality\nCheck problem input file\n\n\n");
		}

#ifdef _DEBUG
		Instrumentor::endSession();
		logFile.close();
#endif //_DEBUG
	}

	// Domain error (like negative in square root)
	catch (std::domain_error& e) {
		std::cerr << "\n\n\nDomain error in function: " << e.what() << "\n\n\n";
		system("pause");
		return 1;
	}

	// invalid_argument -> Something went wrong in the input files
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
		system("pause");
		return 2;
	}

	// length_error -> Required dimensions are not available
	catch (std::length_error& e) {
		std::cerr << "\n\n\nRequired size are invalid: " << e.what() << "\n\n\n";
		system("pause");
		return 3;
	}


	// out_of_range -> Some function tried to acess memory out of its limits
	catch (std::out_of_range& e) {
		std::cerr << "\n\n\nOut of range size: " << e.what() << "\n\n\n";
		system("pause");
		return 4;
	}

	std::cout << OK("\n\nClick a button to exit!");
	std::cin.get();

	return 0;
}