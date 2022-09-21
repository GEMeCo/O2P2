# O2P2 [![GitHub license](https://img.shields.io/github/license/GEMeCo/O2P2?style=for-the-badge)](https://github.com/GEMeCo/O2P2/blob/main/LICENSE) [![GitHub issues](https://img.shields.io/github/issues/GEMeCo/O2P2?style=for-the-badge)](https://github.com/GEMeCo/O2P2/issues) [![GitHub stars](https://img.shields.io/github/stars/GEMeCo/O2P2?style=for-the-badge)](https://github.com/GEMeCo/O2P2/stargazers)
O2P2 is an Object Oriented Platform for the Positional Finite Element Method applied to Thermomechanical Problems. It is under development at Sao Carlos School of Engineering, University of Sao Paulo.

<h1 align="center">
  <img alt="Banner" title="#Banner" height="150" src="./images/Icon.png" />
</h1>

## Copyright Information:
:ballot_box_with_check: This program is free software: you can redistribute it and/or modify it under the terms of the Apache License 2.0 license.

:ballot_box_with_check: This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
It is provided "AS IS". In no event shall the authors be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwire, arising from, out of or in connection with the software or the use of other dealing in the software.

:ballot_box_with_check: This is a research software. There are neither pre-processor nor post-processor (so far). Output was made only for displacements, using AcadView. You may find it [here](https://set.eesc.usp.br/?page_id=237).

:triangular_ruler: None of the developers are **Software Engineers** of any sort. We are **Civil Engineers**, seeking to solve engineering problems. We know that this software has bugs and stuff to improve. We just don't have time to solve them all.

Documentation of this work is under the terms of the Creative Commons Attribution-NonCommerical 4.0 International license.

<h1 align="center">
  <img alt="Banner" title="#Banner" height="100" src="./images/CC-BY-NC.jpg" />
</h1>

## :warning: Under development - Use at your own risk.

Table of Contents
=================
<!--ts-->
   * [About](#about)
   * [Citation](#citation)
   * [Features](#features)
   * [Resources](#resources)
   * [Building and Running](#how-to-run)
   * [Documentation](#documentation)
   * [Project Manager](#project)
   * [Contributors](#contributors)
   * [How to contribute](#how-to-contribute)
   * [Acknowledgement](#acknowledgement)
<!--te-->

About
-----
O2P2 is an object oriented framework for the Positional Finite Element Method applied to Thermomechanical Problems. It was developed for non-linear coupled thermo-mechanical analyzes (:construction: We'll get there, eventually :construction:).

Citation
--------
  Whether it was used in whole or parts, citation is a must! Our software is under development, and a proper presentation paper is underway.
  For now, if you may, please cite the one following papers (as they were used for initial validations):
  - :green_book: CARRAZEDO, R.; CODA, H. B. Triangular based prismatic finite element for the analysis of orthotropic laminated beams, plates and shells. Composite Structures, v. 168, p. 234-246, 2017.
  DOI: [10.1016/j.compstruct.2017.02.027](https://doi.org/10.1016/j.compstruct.2017.02.027)
  - :green_book: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B. Active face prismatic positional finite element for linear and geometrically nonlinear analysis of honeycomb sandwich plates and shells. Composite Structures, v. 200, p. 849-863, 2018.
  DOI: [10.1016/j.compstruct.2018.06.009](https://doi.org/10.1016/j.compstruct.2018.06.009)
  - :green_book: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B.; SALOMÃO, R. C. Vibration and stress analysis of orthotropic laminated panels by active face prismatic finite element. Composite Structures, v. 244, n. 112254, 2020.
  DOI: [10.1016/j.compstruct.2020.112254](https://doi.org/10.1016/j.compstruct.2020.112254)

Features
--------
- v1: Geometric nonlinear mechanical analysis with SVK constitutive model.
	- :boom: only 2D elements;
	- :boom: 2D and 3D elements;
	- :boom: bar elements in 2D and 3D environments;
	- :boom: Parallel processing enabled;

Resources
---------
We are using the following libraries and resources:

- [Eigen C++ Libraries, version 3.4.0](https://eigen.tuxfamily.org/)
- [Microsoft Visual Studio Community 2022, Version 17.3.4](https://visualstudio.microsoft.com/)
- [Intel C++ Compiler, Package ID: w_oneAPI_2022.1.0.256](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html)
- [Doxygen, version 1.9.5](https://doxygen.nl/)

Building and Running
--------------------
1. Clone the source code
`git clone https://github.com/GEMeCo/O2P2.git`

2. Install development dependencies and resources. Just check above.

3. We provided the .vsxproj file. Some configuration may be needed. Good luck :innocent:.

4. Run the executable file and input the project file. Check the documentation below to create the project file.

Documentation
-------------
Software documentation is made directly from annotated sources by doxygen. Use doxywizard to create the documentation. See page `Input files`, and follow the provided instructions to create the project file.

Project Manager
---------------
Rogério Carrazedo

PhD obtained in 2009 from The University of Sao Paulo (USP) at Sao Carlos School of Engineering (EESC).

Joined the Federal University of Technology - Paraná (Brazil) as Assistant Professor in 2010.

Joined the University of Sao Paulo at Sao Carlos School of Engineering, in the Structural Engineering Department in 2015.

Associate Professor since 2020. Has worked on several research projects dealing with static and dynamic behaviour of composite structures using an alternative version of the Finite Element Method, based on Positions.

Contributors
------------
 - [Rafael Correa Salomão](http://lattes.cnpq.br/4408319800130401)
 - [Emerson Felipe Felix](http://lattes.cnpq.br/8352527462118419)
 - [Alexandre Ten Cate Matté](https://lattes.cnpq.br/8144116395864291)
 - [Chiara Pinheiro Teodoro](http://lattes.cnpq.br/6999948388655115)
 - [Thiago da Silva Costa Santos](http://lattes.cnpq.br/6048758348229035)


How to contribute
-----------------
```bash
Sorry, but we are not accepting external contributions, at least yet.
Only students supervised by the head of this project may contribute.
Nevertheless, you may use it as seen fit.
```

Acknowledgement
---------------
<h1 align="center">
  <img alt="Banner" title="#Banner" height="200" src="./images/logo_inst.png" />
</h1>

