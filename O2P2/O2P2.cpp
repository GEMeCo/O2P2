// ================================================================================================
//
// This is the main file. A single object is created, FEAnalysis. It manages the entire process.
// 
// ================================================================================================
// 
// The following configuration is required:
// - Create an environment variable named $eigen_root, and address it to the path of eigen
// - In tab C++/General, additional include directories, type: $(SolutionDir)include;$(Eigen_root)
// 
// - In debug configuration, tab C++/Preprocessor, definitions: include _DEBUG
// 
// - In intel libraries for oneapi, use oneMKL
// 
// ================================================================================================
#define EIGEN_USE_MKL_ALL
#define VERSION 0.2

// Custom Header Files
#include "FEAnalysis.h"
#include "Common.h"
#include "Profiler.h"

/**
  * @brief This main routine creates a single object, FEAnalysis. It manages the entire process.
  * @param argc The number of arguments in the call
  * @param args An array of arguments
  * @return Returns 0 if succeeded
  */
int main(int argc, char** args)
{
	std::string stProj;			  // Project
	std::string stArquivo;		   // Input File
	std::ifstream file;			  // Input File Stream

	// Jabá
	std::cout << "Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved." << std::endl
		<< "Structural Engineering Department / University of Sao Paulo at Sao Carlos School of Engineering." << std::endl << std::endl
		<< "This program comes with ABSOLUTELY NO WARRANTY." << std::endl
		<< "This is a free software. You are welcome to redistribute it under conditions stablished in the license." << std::endl << std::endl;

	// --------------------------------------------------------------------------------------------
	// Check if program call has the project name
	if (argc > 1)
	{
		std::cout << "The following input project was submitted: " << args[1];
		stProj = args[1];
	}
	else
	{
		std::cout << "Project Name: ";
		std::cin >> stProj;

		if (stProj.compare("exit") == 0) return 0;
	}
	stArquivo = stProj + ".txt";

	// Check if there were any execution issue
	try {

#ifdef _DEBUG
		Instrumentor::beginSession("Profile");

		logFile.open(stProj + ".log", std::ios::trunc);
		LOG("Project Name: " + stProj);
#endif //_DEBUG

		// Any instruction inside a try structure may throw and exception
		file.open(stArquivo);
		if (!file) {
			LOG("\n\n\nmain: An exception was thrown when opening file " + stArquivo + "\n\n\n");
			throw std::invalid_argument("\n\n\nAn exception was thrown when opening file " + stArquivo + "\n\n\n");
		}

		// Line from file
		std::string stLine;
		std::string stFlag = "#VERSAO#";

		// Output log in debug mode
		LOG("\nmain: Reading flag " + stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(file, stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		float nVer;
		file >> nVer;

		if (nVer < VERSION) {
			throw std::invalid_argument("\n\n\nObsolete version of problem file.\nUpdate input file and check problem input.\nSomething may be missing\n\n");
		}

		// Look up for flags
		stFlag = "#DIM#";

		// Output log in debug mode
		LOG("\nmain: Reading flag " + stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(file, stLine);
			if (file.eof()) {
				LOG("\n\n\nmain: Reading Error!\nFlag " + stFlag + " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		int nDim;
		file >> nDim;

		if (nDim == 2)
		{
			{
				Timer timer("FEAnalysis.initComponents");
				O2P2::FEAnalysis<2>::initComponents(file);
				file.close();
			}
			{
				Timer timer("FEAnalysis.runAnalysis");
				O2P2::FEAnalysis<2>::runAnalysis();
			}
			{
				Timer timer("FEAnalysis.drawResults");
				O2P2::FEAnalysis<2>::drawResults(stProj);
			}
		}
		else if (nDim == 3) {
			{
				Timer timer("FEAnalysis.initComponents");
				O2P2::FEAnalysis<3>::initComponents(file);
				file.close();
			}
			{
				Timer timer("FEAnalysis.runAnalysis");
				O2P2::FEAnalysis<3>::runAnalysis();
			}
			{
				Timer timer("FEAnalysis.drawResults");
				O2P2::FEAnalysis<3>::drawResults(stProj);
			}
		}
		else {
			LOG("\n\n\nmain: Wrong dimensionality\nCheck problem input file\n\n\n");
			throw std::invalid_argument("\n\n\nWrong dimensionality\nCheck problem input file\n\n\n");
		}

		file.close();

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
		return 1;
	}

	// length_error -> Required dimensions are not available
	catch (std::length_error& e) {
		std::cerr << "\n\n\nRequired size are invalid: " << e.what() << "\n\n\n";
		system("pause");
		return 1;
	}

	// out_of_range -> Some function tried to acess memory out of its limits
	catch (std::out_of_range& e) {
		std::cerr << "\n\n\nOut of range size: " << e.what() << "\n\n\n";
		system("pause");
		return 1;
	}

	std::cout << "\n\nClick a button to exit!";
	std::cin.get();
	std::cin.get();

	return 0;
}
