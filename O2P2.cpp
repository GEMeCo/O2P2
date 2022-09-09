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
#define VERSION 1.2

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

/**
  * @mainpage O2P2, an object oriented environment for the positional finite element method.
  * @image html icon2.png width=100
  * @brief O2P2 is an C++ environment developed for non-linear coupled thermo-mechanical analyzes with the finite element method based on positions.
  *
  * @section first_sec Legal Information
  * @image html logo_inst.png width=300
  *
  * Copyright (C) 2022 Rogério Carrazedo - All Rights Reserved.
  * Structural Engineering Department / University of São Paulo at São Carlos School of Engineering
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  * SOFTWARE.
  *
  * You may use, distribute and even modify this code under terms of Creative Commons Attribution-NonCommerical 4.0 International license.
  * @image html CC-BY-NC.jpg width=150
  *
  * @section developers_sec Developers
  *  @author	[Rogério Carrazedo](http://lattes.cnpq.br/9715974030451423)
  *  @author	[Rafael Correa Salomão](http://lattes.cnpq.br/4408319800130401)
  *  @author	[Emerson Felipe Felix](http://lattes.cnpq.br/8352527462118419)
  *  @author	[Alexandre Ten Cate Matté](https://lattes.cnpq.br/8144116395864291)
  *  @author	[Chiara Pinheiro Teodoro](http://lattes.cnpq.br/6999948388655115)
  *  @author	[Thiago da Silva Costa Santos](http://lattes.cnpq.br/6048758348229035)
  *  @version   1.3.0.1
  *  @date	  2022.09.01
  * 
  * @warning None of the developers are software engineers of any sort. We are civil engineers, seeking to solve engineering problems. We know that this software has lots of bugs. We just don't have time to solve them all.
  *
  * @section citation_sec How to cite:
  *
  * Whether it was used in whole or parts, citation is a must! Our software is under development, and a proper presentation paper is underway.
  * For now, if you may, cite the following papers (since it was used for validations):
  * - Version 1: CARRAZEDO, R.; CODA, H. B. Triangular based prismatic finite element for the analysis of orthotropic laminated beams, plates and shells. Composite Structures, v. 168, p. 234-246, 2017.
  * DOI: <10.1016/j.compstruct.2017.02.027>
  * - Version 2: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B. Active face prismatic positional finite element for linear and geometrically nonlinear analysis of honeycomb sandwich plates and shells. Composite Structures, v. 200, p. 849-863, 2018.
  * DOI: <10.1016/j.compstruct.2018.06.009>
  * - Version 3: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B.; SALOMÃO, R. C. Vibration and stress analysis of orthotropic laminated panels by active face prismatic finite element. Composite Structures, v. 244, n. 112254, 2020.
  * DOI: <10.1016/j.compstruct.2020.112254>
  *
  * O2P2: an object oriented environment for the positional finite element method.
  * Current Version: 1.3.0.1 [S.l.]: SET - EESC - USP, 2022. Available at <https://github.com/GEMeCo/O2P2>
  * DOI:
  *
  * @copyright Licensed under Creative Commons Attribution-NonCommerical 4.0 International license
  *
  * @section akn_sec Funding and Acknowledgement
  *
  * @image html cnpq.png width=150
  * This software received research support from the Brazilian National Council for Scientific and Technological Development
  * (CNPq 428762/2018-2 and CNPq 310564/2018-2) which is gratefully acknowledged.
  *
  * @section install_sec Install Information
  *
  * This program was developed with the following libraries and technologies:
  *
  *	  - Eigen C++ Libraries, Version 3.4.0
  *	  - Microsoft Visual Studio Community 2022, Version 17.2.6
  *	  - Intel Libraries for oneAPI, Package ID: w_oneAPI_2021.2.0.256
  *	  - Intel C++ Compiler, Package ID: w_oneAPI_2021.2.0.256
  *	  - [AcadView](https://set.eesc.usp.br/?page_id=237)
  *
  * @warning Software under development. Improper use will crash your application.
  * @warning This is a research software. There are neither pre-processor nor post-processor. For now, output was made only for displacements, using AcadView.
  * @note A profiler will create a file named "results.json", but only in debug mode. Just open the google chrome in address chrome://tracing and drag/paste the file.
  *
  * @todo 1 - Incluir fluxogramas no detalhamento dos módulos (arquivo O2P2.cpp).
  * @todo 2 - Complementar documentação dos módulos.
  *
  * @defgroup Main_Module O2P2.
  * @brief An object oriented environment for the positional finite element method.
  * @details
  * - Starts the analysis;
  * - Requests to read files and populares containers;
  * - Request to begin the solution process; and
  * - Request to output solution.
  *
  * @defgroup PreProcessor_Module Pre-processor classes.
  * @ingroup Main_Module
  * @brief Read files and populate containers.
  * @details 
  *
  * @defgroup Material Material library
  * @ingroup PreProcessor_Module
  * @brief Materials Library.
  * @details
  *
  * @defgroup Elements Elements library
  * @ingroup PreProcessor_Module
  * @brief Elements Library.
  * @details The following elements are available:
  *
  * @section Bidimensional
  * @subsection tri Triangular
  * @image{inline} html Elem_Tri3.png "Triangular linear" height=150
  * @image{inline} html Elem_Tri6.png "Triangular quadratic" height=150
  * @image{inline} html Elem_Tri10.png "Triangular cubic" height=150
  * @subsection quad Quadrangular
  * @image{inline} html Elem_Quad4.png "Quadrangular linear" height=150
  * @image{inline} html Elem_Quad9.png "Quadrangular quadratic" height=150
  * @image{inline} html Elem_Quad16.png "Quadrangular cubic" height=150
  * @subsection ir_quad Rectangular
  * @image{inline} html Elem_Quad8.png "Quadrangular cubic/linear" height=150
  * @image{inline} html Elem_Quad12.png "Quadrangular cubic/quadratic" height=150
  * @image{inline} html Elem_Quad20png "Quadrangular cubic/quartic" height=150
  *
  * 
  * @section Tridimensional
  * @subsection tet Tetrahedral
  * @image{inline} html Elem_Tet4.png "Tetrahedral linear" height=150
  * @image{inline} html Elem_Tet10.png "Tetrahedral quadratic" height=150
  * @image{inline} html Elem_Tet20.png "Tetrahedral cubic" height=150
  * @subsection hexa Hexahedral
  * @image{inline} html Elem_Hex8.png "Hexahedral linear" height=150
  * @image{inline} html Elem_Hex27.png "Hexahedral quadratic" height=150
  * @image{inline} html Elem_Hex64.png "Hexahedral cubic" height=150
  * @subsection prism Prismatic
  * @image{inline} html Elem_Pri6.png "Prismatic linear" height=150
  * @image{inline} html Elem_Pri18.png "Prismatic quadratic" height=150
  * @image{inline} html Elem_Pri40.png "Prismatic cubic" height=150
  * @subsection ir_prism Irregular prismatic
  * @image{inline} html Elem_Pri20.png "Prismatic cubic/linear" height=150
  * @image{inline} html Elem_Pri30.png "Prismatic cubic/quadratic" height=150
  *
  * @defgroup Processor_Module Processor classes.
  * @ingroup Main_Module
  * @brief Begin the solution processes.
  * @details
  * 
  * @defgroup TimeStep Avaliable time step integration schemes
  * @ingroup Processor_Module
  * @brief Time step integration schemes
  * @details
  *
  * @defgroup NLSolver Avaliable nonlinear solver schemes
  * @ingroup Processor_Module
  * @brief Nonlinear schemes
  * @details

  * @defgroup PostProcessor_Module Post-processor classes.
  * @ingroup Main_Module
  * @brief Output solution for visualization files.
  * @details
  * 
  */
