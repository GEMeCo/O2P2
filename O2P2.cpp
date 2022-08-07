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
    std::string stProj;              // Project
    std::string stArquivo;           // Input File
    std::ifstream file;              // Input File Stream

    // Jabá
    std::cout << "Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved." << std::endl
        << "Structural Engineering Department / University of Sao Paulo at Sao Carlos School of Engineering." << std::endl << std::endl
        << "This program comes with ABSOLUTELY NO WARRANTY." << std::endl
        << "This is a free software. You are welcome to redistribute it under conditions stablished in the license." << std::endl << std::endl;

    // --------------------------------------------------------------------------------------------
    // Check if program call has the project name
    if (argc >= 1)
    {
        std::cout << "Project Name: ";
        std::cin >> stProj;
    }
    else
    {
        std::cout << "The following input project was submitted: " << args[1];
        stProj = args[1];
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
        std::string stFlag = "#DIM#";

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
                FEAnalysis<2>::initComponents(file);
                file.close();
            }
            {
                Timer timer("FEAnalysis.runAnalysis");
                FEAnalysis<2>::runAnalysis();
            }
            {
                Timer timer("FEAnalysis.drawResults");
                FEAnalysis<2>::drawResults(stProj);
            }
        }
        /*else if (nDim == 3) {
            {
                Timer timer("FEAnalysis.initComponents");
                FEAnalysis<3>::initComponents(file);
                file.close();
            }
            {
                Timer timer("FEAnalysis.runAnalysis");
                FEAnalysis<3>::runAnalysis();
            }
            {
                Timer timer("FEAnalysis.drawResults");
                FEAnalysis<3>::drawResults(stProj);
            }
        }*/
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
  * @brief Developed for non-linear coupled thermo-mechanical analyzes with the finite element method based on positions.
  *
  * @section first_sec Legal Information
  * @image html logo_inst.png width=300
  *
  * Copyright (C) 2022 Rogério Carrazedo - All Rights Reserved.
  * Structural Engineering Department / University of São Paulo at São Carlos School of Engineering
  *
  * This program comes with ABSOLUTELY NO WARRANTY.
  * This is free software, and you are welcome to redistribute it under certain conditions.
  *
  * You may use, distribute and even modify this code under terms of Creative Commons Attribution-NonCommerical 4.0 International license.
  * @image html CC-BY-NC.jpg width=150
  *
  * @section developers_sec Developers
  *  @author    Rogério Carrazedo
  *  @version   1.0.0.1
  *  @date      2022.08.01
  *
  * @subsection citation_sec How to cite:
  *
  * Whether it was used in whole or parts, citation is a must!
  *
  * O2P2: an object oriented environment for the positional finite element method.
  * Version: 1.0.0.1. [S.l.]: SET - EESC - USP, 2022. Available at <https://github.com/GEMeCo/O2P2>
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
  *      - Eigen C++ Libraries, Version 3.4.0
  *      - Microsoft Visual Studio Community 2022, Version 17.2.1
  *      - Intel Libraries for oneAPI, Package ID: w_oneAPI_2021.2.0.243
  *      - Intel C++ Compiler, Package ID: w_oneAPI_2021.2.0.243
  *
  * @warning Software under development. Improper use will crash your application.
  *
  * @bug     Os dynamic_cast não tem verificação posterior se estão retornando pointeiros nulos (poderiam ser static_cast - mais rápidos).
  * 
  * @todo No arquivo O2P2.cpp, incluir fluxogramas no detalhamento dos módulos.
  *
  * @defgroup Main_Module O2P2.
  * @brief An object oriented environment for the positional finite element method.
  * @details Starts the analysis;
  *
  * Requests to read files and populares containers;
  *
  * Request to begin the solution process; and
  *
  * Request to output solution.
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
  * @details
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
