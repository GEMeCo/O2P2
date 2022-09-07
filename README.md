# O2P2
## Object Oriented Platform for the Positional Finite Element Method applied to Thermomechanical Problems
### An object oriented environment for the positional finite element method.

<h1 align="center">
  <img alt="Banner" title="#Banner" height="200" src="./images/Icon.png" />
</h1>

## Copyright Information:
Content on this work is licensed under a Creative Commons Attribution-NonCommerical 4.0 International license.

<h1 align="center">
  <img alt="Banner" title="#Banner" height="100" src="./images/CC-BY-NC.jpg" />
</h1>

This program is free software: you can redistribute it and/or modify it under the terms of the Creative Commons Attribution-NonCommerical 4.0 International license..

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
It is provided "AS IS". In no event shall the authors be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwire, arising from, out of or in connection with the software or the use of other dealing in the software.

Also, this is a research software. There are neither pre-processor nor post-processor. For now, output was made only for displacements, using AcadView.

## Under development - Use at your own risk.

Table of Contents
=================
<!--ts-->
   * [About](#about)
   * [Citation](#citation)
   * [Current Features](#features)
   * [Resources](#resources)
   * [Author](#author)
   * [Contributors](#contributors)
   * [How to contribute](#how-to-contribute)
<!--te-->

About
-----
Object oriented platform for the Positional Finite Element Method applied to Thermomechanical Problems, an object oriented framework for the finite element method based on positions. Developed for non-linear coupled thermo-mechanical analyzes.

Citation
--------
  Whether it was used in whole or parts, citation is a must! Our software is under development, and a proper presentation paper is underway.
  For now, if you may, cite the following papers (since it was used for validations):
  - Version 1: CARRAZEDO, R.; CODA, H. B. Triangular based prismatic finite element for the analysis of orthotropic laminated beams, plates and shells. Composite Structures, v. 168, p. 234-246, 2017.
  DOI: <10.1016/j.compstruct.2017.02.027>
  - Version 2: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B. Active face prismatic positional finite element for linear and geometrically nonlinear analysis of honeycomb sandwich plates and shells. Composite Structures, v. 200, p. 849-863, 2018.
  DOI: <10.1016/j.compstruct.2018.06.009>
  - Version 3: CARRAZEDO, R.; PACCOLA, R. R.; CODA, H. B.; SALOMÃO, R. C. Vibration and stress analysis of orthotropic laminated panels by active face prismatic finite element. Composite Structures, v. 244, n. 112254, 2020.
  DOI: <10.1016/j.compstruct.2020.112254>


Features
--------
- v1: Geometric nonlinear analysis with SVK constitutive model.
	- v1.1: only 2D elements;
	- v1.2: 2D and 3D elements;
	- v1.3: bar elements in 2D and 3D environments;
	- v1.4: Parallel processing enabled;

Resources
---------
The following libraries have to be pre-instaled:

- [Eigen C++ Libraries, version 3.4.0](https://eigen.tuxfamily.org/)
- [Microsoft Visual Studio Community 2019, Version 16.11.14](https://visualstudio.microsoft.com/)
- [Intel C++ Compiler, Package ID: w_oneAPI_2021.2.0.243](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html)

Author
------
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
We are not yet accepting external contributions. For now, only students supervised by the head of this project may contribute. Nevertheless, you may use it as seen fit.
```

Acknowledgement
---------------
<h1 align="center">
  <img alt="Banner" title="#Banner" height="200" src="./images/logo_inst.png" />
</h1>

This software received research support from the Brazilian National Council for Scientific and Technological Development (CNPq 428762/2018-2 and CNPq 310564/2018-2) which is gratefully acknowledged.
