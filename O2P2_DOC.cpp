
/**
  * @mainpage O2P2, an object oriented environment for the positional finite element method.
  * @image html icon2.png width=100
  * @brief O2P2 is an C++ environment developed for non-linear coupled thermo-mechanical analyzes with the finite element method based on positions.
  * @tableofcontents
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
  * @section citation_sec How to cite
  *
  * Whether it was used in whole or parts, citation is a must! Our software is under development, and a proper presentation paper is underway.
  * For now, if you may, please cite the following papers (since they were used for validations):
  * - CARRAZEDO, R.; CODA, H. B. Triangular based prismatic finite element for the analysis of orthotropic laminated beams, plates and shells. Composite Structures, v. 168, p. 234-246, 2017.
  * DOI: <10.1016/j.compstruct.2017.02.027>
  * - CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B. Active face prismatic positional finite element for linear and geometrically nonlinear analysis of honeycomb sandwich plates and shells. Composite Structures, v. 200, p. 849-863, 2018.
  * DOI: <10.1016/j.compstruct.2018.06.009>
  * - CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B.; SALOMÃO, R. C. Vibration and stress analysis of orthotropic laminated panels by active face prismatic finite element. Composite Structures, v. 244, n. 112254, 2020.
  * DOI: <10.1016/j.compstruct.2020.112254>
  *
  * Current Version: 1.3.0.1 [S.l.]: SET - EESC - USP, 2022. Available at <https://github.com/GEMeCo/O2P2>\n
  * DOI: Soon
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
  * @page inputfiles Input Files
  * @tableofcontents
  * You must provide five files containing all information required to create the domain and the analysis.
  * Information on files are divided by #FLAGS#, after which data must follow a required pattern.
  *
  * @subsection pf Project file
  * The following information must be provided:
  * - A control version of the input file, associated to current O2P2 version. Obsolete versions must be checked, for something might be missing.
  * @verbatim
    #VERSAO#
    1.2
    @endverbatim
  *
  * - Problem Dimensionality (either 2D or 3D)
  * @verbatim
    #DIM#
    2
    @endverbatim
  *
  * - Four files associated to the project (nodal, element, material and dof).\n
  * @verbatim
    #ARQUIVOS#
    nodal.txt
    element.txt
    material.txt
    dof.txt
    @endverbatim
  *
  * In a single line, input the following:
  * - Time integration scheme (1 for quasi-static, the only available for now)\n
  * - Nonlinear Solver: (1 - Newton-Raphson)\n
  * - Type of analysis: (1 - Mechanic)\n
  * - Tolerance (double)\n
  * - Minimum number of interactions\n
  * - Maximum number of interactions\n
  * - Number of load steps
  * @verbatim
    #ANALISE#
    1	1	1	0.0000001	3	10	1
    @endverbatim
  *
  * @subsection nf Node file
  * - Number of nodes
  * @verbatim
    #DADOS_ND#
    3
    @endverbatim
  *
  * - Number (index)
  * - Nodes coordinates (x and y for 2D; x, y and z for 3D)
  * @verbatim
    #NOS#
    1   0.0 0.0
    2   1.0 0.0
    3   0.5 1.0
    @endverbatim
  * 
  * @subsection ef Element file
  * - Number of elements and sections (if any)
  * @verbatim
    #DADOS_EL#
    1   1
    @endverbatim
  * 
  * - Type of cross section (1 - area; 2 - thickness)
  * - Value
  * - Plane state (1 - Plane stress; 2 - Plane strain)
  * @verbatim
    #SECOES#
    2   1.  1
    @endverbatim
  * 
  * - Number (index)
  * - Type of element (1 - bar; 2 - triangular; 3 - rectangular; 4 - tetrahedral; 5 - hexahedral; and 6 - prism)
  * - Order (1 - linear; 2 - quadratic; and 3 - cubic)
  * - Number of integration points (see Elem)
  * - Material number
  * - Section (if any)
  * - Conectivity (1 .. number of nodes)
  * @verbatim
    #ELEMENTOS#
    1   2   1   3   1   1   1   2   3
    @endverbatim
  * 
  * @subsection mf Material file
  * - Number of different materials
  * @verbatim
    #MATERIAIS#
    1
    @endverbatim
  * 
  * Material information
  * - material index
  * - Type (1 - elastic SVK)
  * Parameters
  * - Young modulus; Poisson ratio; Density; Damping coefficient
  * @verbatim
    #PARAMETROS#
    1   1
    1.  0.  0.  0.
    #
    @endverbatim
  * Notice that `#` is mandatory after each group of material information
  * 
  * @subsection df DOF file
  * Boundary conditions for each load step
  * - Number of time steps
  * - Time step
  * - Number of Dirichlet boundary contitions
  * - Number of Neumann boundary conditions
  * @verbatim
    #FC1#
    1   .5  3   1
    @endverbatim
  * 
  * Dirichlet boundary conditions
  * - Node
  * - Direction
  * - Value
  * - Time behaviour (var[0] + var[1].t + var[2].t^2)
  * @verbatim
    #DIR1#
    1   1   0.  1.  0.  0.
    1   2   0.  1.  0.  0.
    2   1   0.  1.  0.  0.
    @endverbatim
  * 
  * Neumann boundary condition
  * - Node
  * - Direction
  * - Value
  * - Time behaviour (var[0] + var[1].t + var[2].t^2)
  * @verbatim
    #NEU1#
    3   1   0.1  2.  0.  0.
    @endverbatim
  * 
  * These flags must be redefined for each and every load step (#FC2#; #DIR2#; #NEU2#; and so on).
  * 
  * @defgroup Main_Module O2P2.
  * @brief An object oriented environment for the positional finite element method.
  * @details
  * - Starts the analysis creating an FEAnalysis;
  * - Requests to read files and populares containers with FEAnalysis::initComponents;
  * - Request to begin the solution process with FEAnalysis::runAnalysis; and
  * - Request to output solution with FEAnalysis::drawResults.
  *
  * @image html O2P2_Class_Chart.png width=750
  *
  * @sa O2P2::FEAnalysis
  *
  * @defgroup PreProcessor_Module Pre-processor classes.
  * @ingroup Main_Module
  * @brief Read files and populate containers.
  * @details Once input files are defined, reads their contents and creates the geometry domain.
  *
  * @image html O2P2_Preprocessor_Module.png width=600
  *
  * @sa O2P2::ModelBuilder
  * @sa O2P2::Prep::Domain
  * @sa Elements
  * @sa Material
  *
  * @defgroup Material Material library
  * @ingroup PreProcessor_Module
  * @brief Materials Library.
  * @details The following constitutive models are available:
  * - Saint-Venant-Kirchhoff (isotropic).
  * - That´s all, for now :)
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
  * @image{inline} html Elem_Quad20.png "Quadrangular cubic/quartic" height=150
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
  * @details Aggregation of classes to evaluate the analysis.
  *
  * @image html O2P2_Processor_Module.png width=1100
  *
  * @sa O2P2::Proc::SolutionAlgorithm
  * @sa O2P2::Proc::TimeStepping
  * @sa O2P2::Proc::NonLinearSolver
  * @sa O2P2::Proc::Mesh
  *
  * @defgroup TimeStep Time step integration schemes
  * @ingroup Processor_Module
  * @brief Time step integration schemes
  * @details The following time integrantion schemes are available:
  * - Quasi-static time integration scheme.
  * - That's all, so far.
  *
  * @defgroup NLSolver Nonlinear solver schemes
  * @ingroup Processor_Module
  * @brief Nonlinear schemes
  * @details The following nonlinear solver schemes are available:
  * - Newton-Raphson.
  * - That's all, so far.
  *
  * @defgroup PostProcessor_Module Post-processor classes.
  * @ingroup Main_Module
  * @brief Output solution for visualization files.
  * @details Under development.
  *
  */
