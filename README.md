Here we maintain an updateable catalog of open source molecular modeling software, initially taken from our [paper](https://www.authorea.com/users/26858/articles/66197).

Eventually we will deploy a less monolithic document with additional features (such as sorting and filtering), correct citations, and a better layout.  

***Please contribute edits by forking the repository and submitting a pull request***.


  * [Methods](#methods)
      * [Development Activity](#development-activity)
      * [Usage Activity](#usage-activity)
  * [Cheminformatics](#cheminformatics)
    * [Toolkits](#toolkits)
    * [Standalone Programs](#standalone-programs)
    * [Graphical Development Environments](#graphical-development-environments)
  * [Visualization](#visualization)
    * [2D Desktop Applications (Table [2ddesktopviz])](#2d-desktop-applications-table2ddesktopviz)
    * [3D Desktop Applications](#3d-desktop-applications)
    * [Web-Based Visualization](#web-based-visualization)
  * [QSAR/ADMET Modeling](#qsaradmet-modeling)
    * [Descriptor Calculators](#descriptor-calculators)
    * [Model Building](#model-building)
    * [Model Application](#model-application)
    * [Visualization](#visualization-1)
  * [Quantum Chemistry](#quantum-chemistry)
    * [<em>Ab initio</em> Calcuation](#ab-initio-calcuation)
    * [Helper Applications](#helper-applications)
    * [Visualization](#visualization-2)
  * [Ligand Dynamics and Free Energy Calculations](#ligand-dynamics-and-free-energy-calculations)
    * [Simulation Software](#simulation-software)
    * [Simulation Setup and Analysis](#simulation-setup-and-analysis)
  * [Virtual Screening and Ligand Design](#virtual-screening-and-ligand-design)
    * [Ligand-Based](#ligand-based)
    * [Docking and Scoring](#docking-and-scoring)
    * [Pocket Detection](#pocket-detection)
    * [Ligand Design](#ligand-design)


Methods
=======

For every identified software package, we report its primary URL and software license and assign it an activity code. For simplicity, BSD-like licenses (e.g. NCSA) are reported as BSD. Activity codes consist of a development activity level (alphabetical) and usage activity level (numerical). Activity codes were assigned as follows:

### Development Activity

A  
Substantial development (e.g. a new major release, the addition of new features, or substantial refinements of existing features) within the last 18 months. Note this includes all projects that were created in the last 18 months.

B  
Evidence of some development within the last 18 months such as a minor release or bug fixes to a development branch.

C  
No evidence of development (changes to the source code or documentation) within the last 18 months. Note that in cases where a package does not follow an open development model (i.e., source is only released with official releases) the estimate of development activity will be overly conservative.

### Usage Activity

1.  Substantial user usage within the last 18 months (more than 20 downloads a month on average from SourceForge, more than 20 stars or forks on GitHub, more than 10 citations a year, and/or a clearly active user community as indicated by traffic on mailing lists or discussion boards).

2.  Moderate user usage within the last 18 months.

3.  Minimal or no identifiable user usage within the last 18 months (fewer than 50 downloads total on SourceForge, three or fewer stars and/or forks on GitHub, or fewer than one citation a year).

We omit some packages with extended periods of inactivity (e.g. more than 10 years) where there is little evidence of any usage or packages that are referenced in the literature but for which we could not find a extant source code repository. We also omit packages that provide common and/or trivial functionality (e.g. molecular weight calculators) and those that require non-open source packages in order to function.

Cheminformatics
===============


Toolkits 
-----------------------------

| Name                 | URL                                                                      |  License | Activity | Citation |
|:---------------------|:-------------------------------------------------------------------------|:--------:|:--------:|:--------:|
| BALL                 | <http://www.ball-project.org/>                                           |   LGPL   |    A2    |          |
| CDK                  | <https://sourceforge.net/projects/cdk>                                   |   LGPL   |    A1    |          |
| Chem**<sup>*f*</sup> | <https://github.com/stefan-hoeck/chemf>                                  |    GPL   |    C2    |          |
| chemfp               | <http://chemfp.com>                                                      |    MIT   |    C3    |          |
| chemkit              | <https://github.com/kylelutz/chemkit>                                    |    BSD   |    B1    |          |
| ChemmineR            | <https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html> | Artistic |    A1    |          |
| Cinfony              | <https://github.com/cinfony/cinfony>                                     |  BSD/GPL |    B1    |          |
| CurlySMILES          | <http://www.axeleratio.com/csm/proj/main.htm>                            |    GPL   |    C2    |          |
| DisCuS               | <https://github.com/mwojcikowski/discus>                                 |    GPL   |    B3    |          |
| Fafoom               | <https://github.com/adrianasupady/fafoom>                                |   LGPL   |    A2    |          |
| fmcsR                | <http://www.bioconductor.org/packages/fmcsR>                             | Artistic |    A1    |          |
| frowns               | <http://frowns.sourceforge.net>                                          |  Python  |    C2    |          |
| Helium               | <http://www.moldb.net/helium.html>                                       |    BSD   |    B2    |          |
| Indigo               | <http://lifescience.opensource.epam.com/indigo>                          |    GPL   |    A1    |          |
| JoeLib               | <http://sourceforge.net/projects/joelib>                                 |    GPL   |    C1    |          |
| LICSS                | <https://github.com/KevinLawson/excel-cdk>                               |    GPL   |    A2    |          |
| MayaChemTools        | <http://www.mayachemtools.org>                                           |   LGPL   |    A2    |          |
| Mychem               | <http://mychem.sourceforge.net>                                          |    GPL   |    B2    |          |
| ODDT                 | <https://github.com/oddt/oddt>                                           |    BSD   |    A2    |          |
| Open Babel           | <http://openbabel.org>                                                   |    GPL   |    A1    |          |
| OPSIN                | <http://opsin.ch.cam.ac.uk>                                              | Artistic |    A1    |          |
| OrChem               | <http://orchem.sourceforge.net>                                          |   LGPL   |    C2    |          |
| osra                 | <http://sourceforge.net/projects/osra>                                   |    GPL   |    A1    |          |
| OUCH                 | <https://github.com/odj/Ouch>                                            |    GPL   |    C2    |          |
| pybel                | <http://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html>          |    GPL   |    A1    |          |
| rcdk                 | <https://cran.r-project.org/web/packages/rcdk>                           |   LGPL   |    B2    |          |
| RDKit                | <http://www.rdkit.org>                                                   |    BSD   |    A1    |          |
| RInChI               | <http://www-rinchi.ch.cam.ac.uk>                                         |  Apache  |    A3    |          |
| rpubchem             | <https://r-forge.r-project.org/projects/rpubchem>                        |    GPL   |    C3    |          |
| rubabel              | <https://github.com/princelab/rubabel>                                   |    MIT   |    C2    |          |
| SMSD                 | <http://www.ebi.ac.uk/thornton-srv/software/SMSD>                        |   CCAL   |    B2    |          |
| Som-it<sup>TM</sup>  | <http://silicos-it.be>                                                   |   LGPL   |    C3    |          |
| webchem              | <https://github.com/ropensci/webchem>                                    |    MIT   |    A2    |          |


The Biochemical Algorithms Library (BALL) provides an object-oriented C++ library for structural bioinformatics, and its capabilities include molecular mechanics, support for reading and writing a variety of file formats, protein-ligand scoring, docking, and QSAR modeling.

The Chemistry Development Kit (CDK) is a cheminformatics toolkit written in Java. Its capabilities include support for reading and writing a variety of chemical formats, descriptor and fingerprint calculation, force field calculations, substructure search, and structure generation.

Chem**<sup>*f*</sup> is a minimal cheminformatics toolkit written in the functional language Scala.

chemfp is a high-performance library with a Python interface for generating and searching for molecular fingerprints.

chemkit is a C++ cheminformatics toolkit that includes support for visualization with the Qt framework and molecular modeling.

ChemmineR is a cheminformatics package for the R statistical programming language that is built using Open Babel. Its capabilities include property calculations, similarity search, and classification and clustering of compounds.

Cinfony provides a single, simple standardized interface to other cheminformatics toolkits, including Open Babel, RDKit, the CDK, Indigo, JChem, OPSIN, and several web services.

CurlySMILES provides parsing functionality for an extension of the SMILES format that supports the description of complex molecular systems.

DisCuS (Database System for Compound Selection) provides support for analyzing the results of a high-throughput screen.

Fafoom (flexible algorithm for optimization of molecules) is a Python library for identifying low energy conformers using a genetic algorithm.

fmcsR is an R package that efficiently performs flexible maximum common substructure matching that allows minor mismatches between atoms and bonds in the common substructure.

Frowns is a cheminformatics toolkit mostly written in Python that provides basic support for SMILES and SD files, SMARTS search, fingerprint generation, and property perception.

Helium is a cheminformatics toolkit written using modern C++ idioms that provides support for SMILES files, fingerprints generation, and SMARTS and SMIRKS.

Indigo is a cheminformatics toolkit written in C++ with C, Python, Java (including a KNIME node), and C\# bindings. Its capabilities include general support for manipulating molecules, property calculation, combinatorial chemistry, scaffold detection and R-group decomposition, reaction processing, substructure matching and similarity search.

JOELib is a cheminformatics toolkit written in Java. Its capabilities include SMARTS substructure search, descriptor calculation, and processing/filtering pipes.

LICSS integrates with the CDK to provide representations and analysis of chemical data embedded within Microsoft Excel.

MayChemTools is a collection of Perl scripts for manipulating chemical data, interfacing with databases, generating fingerprints, performing similarity search, and computing molecular properties.

Mychem is built using OpenBabel and provides an extension to the MySQL database package that adds the ability to search, analyze, and convert chemical data within a MySQL database.

The Open Drug Discovery Toolkit (ODDT) is entirely written in Python, is built on top of RDKit and Open Babel, and is focused on providing enhanced functionality for managing and implementing drug discovery workflows, such as making it easy to implement a docking pipeline.

Open Babel is substantial cheminformatics toolkit written in C++ with Python, Perl, Java, Ruby, R, PHP, and Scala bindings. Its capabilities include support for more than 100 chemical file formats, fingerprint generation, property determination, similarity and substructure search, structure generation, and molecular force fields. It has absorbed the Confab conformer generator which produces 3D structures through the systematic enumeration of torsions and energy minimization.

OPSIN , the Open Parser for Systematic IUPAC nomenclature, converts plain-text chemical nomenclature to machine readable CML or InChi formats.

OrChem is built using the CDK and provides an extension to Oracle databases that adds the ability to incorporate and search chemical data.

OSRA provides optical structure recognition. It takes as input an image and generates a SMILES string.

Ouch (Ouch Uses Chemical Haskell) is a minimal cheminformatics toolkit written in the functional language Haskell.

Pybel provides the full functionality of Open Babel, but common routines are provided in a simplified, more ‘pythonic’ interface.

rcdk provides an R interface to the CDK and working with fingerprints.

RDKit is a substantial cheminformatics toolkit written in C++ with Python, Java and C\# bindings. Its capabilities include file handling, manipulation of molecular data, chemical reactions, substantial support for fingerprinting, substructure and similarity search, 3D conformer generation, property determination, force field support, shape-based alignment and screening, and integration with PyMOL, KNIME, and PostgreSQL.

RInChI provides tools for creating and manipulating reaction InChIs, a unique string for describing a reaction.

rpubchem is an R package for interfacing with the PubChem database.

rubabel is similar to Pybel in that it provides a native Ruby interface to Open Babel.

The Small Molecule Subgraph Detector (SMSD) is a Java library for calculating the maximum common subgraph between small molecules.

Som-it<sup>TM</sup> is an R package for creating and visualizing self-organizing maps from large datasets.

webchem is an R package for interfacing with a dozen different on-line resources for chemical data.

Standalone Programs 
-------------------

| Name                   | URL                                                         | License | Activity | Citation |
|:-----------------------|:------------------------------------------------------------|:-------:|:--------:|:--------:|
| cApp                   | <http://www.structuralchemistry.org/pcsb>                   |   GPL   |    A3    |          |
| checkmol/matchmol      | <http://merian.pch.univie.ac.at/~nhaider/cheminf/cmmm.html> |   GPL   |    C3    |          |
| ConvertMAS             | <http://sourceforge.net/projects/convertmas>                |   GPL   |    A3    |          |
| Filter-it<sup>TM</sup> | <http://silicos-it.be>                                      |   LGPL  |    C3    |          |
| Frog2                  | <https://github.com/tuffery/Frog2>                          |   GPL   |    B2    |          |
| LMR                    | <https://github.com/IanAWatson/Lilly-Medchem-Rules>         |   GPL   |    B2    |          |
| Molpher                | <https://www.assembla.com/spaces/molpher/wiki>              |   GPL   |    C2    |          |
| MoSS                   | <http://www.borgelt.net/moss.html>                          |   MIT   |    A2    |          |
| OMG                    | <http://sourceforge.net/projects/openmg>                    |   GPL   |    C1    |          |
| sdf2xyz2sdf            | <http://sdf2xyz2sdf.sourceforge.net>                        |   GPL   |    C2    |          |
| sdsorter               | <https://sourceforge.net/projects/sdsorter>                 |   GPL   |    B3    |          |
| Shape                  | <http://sourceforge.net/projects/shapega>                   |   GPL   |    C3    |          |
| Strip-it<sup>TM</sup>  | <http://silicos-it.be>                                      |   LGPL  |    C3    |          |


cApp is a Java application that provides tools for evaluating physico-chemical properties, performing similarity searches, and querying the PubChem database.

The utilities checkmol and matchmol decompose and compare functional groups of input molecules.

ConvertMAS is a utility for converting between formats and merging and splitting multi-molecule files.

Filter-it<sup>TM</sup> filters a set of molecules based on their properties such as physicochemical parameters and graph-based properties.

Frog2 uses a two stage Monte Carlo approach coupled with energy minimization to rapidly generate 3D conformers.

The Lilly MedChem Rules (LMR) apply filters to avoid reactive and promiscuous compounds.

Molpher generates a virtual chemical library that represents the chemical space between two input molecules as it consists of the path found by morphing one molecule to another.

MoSS (Molecular Subsstructure miner) finds common molecular substructures and discriminative fragments within a compound library.

The Open Molecule Generator (OMG) enumerates all possible chemical structures given constraints on their composition.

sdf2xyz2sdf converts between SDF and TINKER XYZ files.

sdsorter provides convenient routines for manipulating, sorting, and filtering the contents of sdf molecular data files based on the embedded sd data tags.

Shape employs a genetic algorithm to generate conformations of carbohydrates.

Strip-it<sup>TM</sup> is built using Open Babel and extracts molecular scaffolds.


Graphical Development Environments 
----------------------------------

| Name        | URL                                        |  License | Activity | Citation |
|:------------|:-------------------------------------------|:--------:|:--------:|:--------:|
| AMBIT       | <http://ambit.sourceforge.net>             |    GPL   |    A1    |          |
| Bioclipse   | <http://www.bioclipse.net>                 |  Eclipse |    B1    |          |
| Galaxy Tool | <https://github.com/bgruening/galaxytools> | Academic |    A1    |          |
| KNIME       | <https://www.knime.org>                    |    GPL   |    A1    |          |
| Orange      | [orange.biolab.si](orange.biolab.si)       |    BSD   |    A1    |          |
| SA2         | <http://sa2.sourceforge.net>               |    GPL   |    A1    |          |
| Taverna     | <http://www.taverna.org.uk>                |   LGPL   |    A1    |          |
| Weka        | <https://sourceforge.net/projects/weka>    |    GPL   |    A1    |          |



Ambit integrates with the CDK to provide web-based applications for chemical search and analysis and includes a tautomer generation algorithm .

Bioclipse is a workbench, based on the Eclipse framework, for manipulating and analyzing biochemical data and databases. It integrates with the CDK and Jmol to provide cheminformatic functionality and also has modules for bioinformatics (primarly sequence analysis) and QSAR modeling.

Galaxy is a web platform for exploring biomedical data and includes as a component a Chemical Toolbox that integrates a number of other cheminformatics tools to offer an array of molecular search, property calculation, clustering, and manipulation capabilities.

The Konstanz Information Miner (KNIME) is a general workflow environment that includes a number of plugins for cheminformatics, such as CDK and RDKit modules, as well as bioinformatics and machine learning modules.

Orange is a graphical interface for construction interactive workflows and performing data analysis and visualization.

Screening Assistant 2 (SA2) is a GUI written in Java that integrates with other toolkits to help manage, analyze, and visualize libraries of compounds.

Taverna is a graphical workflow editor that includes support for integrating with web services and the CDK .

Weka is a platform for data mining and machine learning that can be adapted for cheminformatics.


Visualization
=============


2D Desktop Applications (Table \[2ddesktopviz\])
------------------------------------------------

| Name        | URL                                                                  | License | Activity | Citation |
|:------------|:---------------------------------------------------------------------|:-------:|:--------:|:--------:|
| BKchem      | <http://bkchem.zirael.org>                                           |   GPL   |    C3    |          |
| chemfig     | <https://www.ctan.org/pkg/chemfig>                                   |  LaTeX  |    A2    |          |
| Chemtool    | <http://ruby.chemie.uni-freiburg.de/~martin/chemtool>                |   GPL   |    B3    |          |
| JChemPaint  | <http://jchempaint.github.io>                                        |   LGPL  |    B1    |          |
| LeView      | <http://www.pegase-biosciences.com/leview-ligand-environment-viewer> |   GPL   |    B3    |          |
| mol2chemfig | <http://chimpsky.uwaterloo.ca/mol2chemfig>                           |  LaTeX  |    C3    |          |
| Molsketch   | <http://sourceforge.net/projects/molsketch>                          |   GPL   |    A1    |          |
| SketchEl    | <http://sketchel.sourceforge.net>                                    |   GPL   |    A1    |          |


BKChem is a 2D molecular editor written in python that uses the Tk GUI toolkit.

chemfig is a tool for embedding chemical drawings in <span>LaTeX</span> documents.

Chemtool is a 2D molecular editor for Linux systems that uses the GTK toolkit.

JChemPaint is a Java-based 2D molecular editor built using the CDK toolkit.

LeView generates 2D representations of ligand-protein interactions that highlight features such as hydrogen bonds.

mol2chemfig converts SMILES files into <span>LaTeX</span> source code.

Molsketch is a 2D molecular editor written in C++ with the Qt toolkit that includes support for the Windows and Android operating systems.

SketchEl is a Java-based 2D molecular editor that includes support for a datasheet view for handling multi-molecule files.

3D Desktop Applications 
-----------------------

| Name           | URL                                                  |  License | Activity | Citation |
|:---------------|:-----------------------------------------------------|:--------:|:--------:|:--------:|
| Avogadro       | <http://avogadro.cc>                                 |    GPL   |    A1    |          |
| BALLView       | <http://www.ball-project.org/ballview>               |   LPGL   |    A2    |          |
| gMol           | <https://github.com/tjod/gMol/wiki>                  |    GPL   |    A3    |          |
| Jamberoo       | <https://sourceforge.net/projects/jbonzer>           |   LGPL   |    A3    |          |
| LPMV           | <https://sourceforge.net/projects/lpmolecularviewer> |   LGPL   |    B3    |          |
| Luscus         | <https://sourceforge.net/projects/luscus>            | Academic |    A1    |          |
| Molecular Rift | <https://github.com/Magnusnorrby/MolecularRift>      |    GPL   |    A3    |          |
| OpenStructure  | <http://www.openstructure.org>                       |   LGPL   |    A2    |          |
| PLIP           | <https://github.com/ssalentin/plip>                  |  Apache  |    A2    |          |
| PyMOL          | <https://sourceforge.net/projects/pymol>             |  Python  |    A1    |          |
| RasTop         | <https://sourceforge.net/projects/rastop>            |    GPL   |    C1    |          |
| OpenRasMol     | <https://sourceforge.net/projects/openrasmol>        |    GPL   |    C1    |          |
| SPADE          | <http://www.spadeweb.org>                            |    BSD   |    C3    |          |
| QuteMol        | <http://qutemol.sourceforge.net>                     |    GPL   |    C1    |          |

Avogadro is a 3D molecular viewer and editor with a modular plugin architecture with both Python and C++ bindings that includes interactive structure optimization for real-time editing.

BALLView provides interactive 3D visualizations as part of the BALL cheminformatics toolkit.

gMol provides basic interactive 3D visualizations of molecular data readable by Open Babel.

Jamberoo provides a basic Java-based 3D molecular viewer and editor.

LP Molecular Viewer is an ActiveX/ATL control for embedding interactive 3D representations of molecular data in Microsoft products.

Luscus is a 3D viewer and editor that is designed with a focus on electronic structure information.

Molecular Rift integrates with the Oculus Rift virtual reality headset to provide immersive visualization of 3D molecular data.

OpenStructure is a computational structural biology framework that provides a 3D viewer for manipulating structural information and includes an interactive Python shell.

PLIP (Protein-Ligand Interaction Profiler) runs as a web application and analyzes and visualizes protein-ligand interactions in 3D.

PyMOL is a substantial 3D molecular viewer that includes a full Python interface to support scripting and plugin development.

RasTop and OpenRasMol are based off the venerable RasMol software and provide basic 3D visualization.

SPADE (Structural Proteomics Application Development Environment) is a graphical Python interface for structural informatics.

QuteMol provides high-quality, visually engaging renderings of 3D molecular data.

Web-Based Visualization
-----------------------

| Name      | URL                                       | License | Activity | Citation |
|:----------|:------------------------------------------|:-------:|:--------:|:--------:|
| 3Dmol.js  | [3dmol.csb.pitt.edu](3dmol.csb.pitt.edu)  |   BSD   |    A1    |          |
| CH5M3D    | <https://sourceforge.net/projects/ch5m3d> |   GPL   |    C1    |          |
| Chemozart | <https://chemozart.com>                   |  Apache |    A1    |          |
| CWC       | <https://web.chemdoodle.com>              |   GPL   |    A1    |          |
| JSME      | <http://peter-ertl.com/jsme>              |   BSD   |    A1    |          |
| Jmol      | <http://jmol.sourceforge.net>             |   LGPL  |    A1    |          |
| JSmol     | <https://sourceforge.net/projects/jsmol>  |   LGPL  |    A1    |          |
| NGL       | <http://proteinformatics.charite.de/ngl>  |   MIT   |    A1    |          |
| PV        | <https://biasmv.github.io/pv>             |   MIT   |    A1    |          |

3Dmol.js is a JavaScript library that provides WebGL-accelerated interactive 3D visualizations of molecular structures and surfaces.

CH5M3D uses JavaScript and HTML5 to provide visualization and editing of 3D structures of small molecules.

Chemozart is a WebGL-based web application for 3D editing of small molecules.

CWC (ChemDoodle Web Components) provides a suite of web-based visualizers and editors for 2D and 3D molecular data.

JSME is a pure JavaScript 2D molecular editor that can export and import SMILES data.

Jmol is a Java applet for interactive 3D visualization that provides significant cheminformatics support and a custom scripting language.

JSmol is the JavaScript port of Jmol and does not require the Java plugin to run.

NGL is a WebGL-accelerated viewer and JavaScript library for interactive 3D visualization of macromolecules.

PV (Protein Viewer) is a WebGL-accelerated viewer for interactive 3D visualization of macromolecules with a functional-style API.


QSAR/ADMET Modeling
===================


Descriptor Calculators
----------------------

| Name             | URL                                                   |    License    | Activity | Citation |
|:-----------------|:------------------------------------------------------|:-------------:|:--------:|:--------:|
| 4D-FAP           | <http://www.ra.cs.uni-tuebingen.de/software/4DFAP>    |      LGPL     |    C2    |          |
| BlueDesc         | <http://www.ra.cs.uni-tuebingen.de/software/bluedesc> |      GPL      |    C3    |          |
| MolSig           | <http://molsig.sourceforge.net>                       |      GPL      |    C2    |          |
| PaDEL-descriptor | <http://www.yapcwsoft.com/dd/padeldescriptor>         | Public Domain |    C1    |          |
| TMACC            | <http://comp.chem.nottingham.ac.uk/download/tmacc>    |      GPL      |    C2    |          |


4D Flexible Atom-Pair Kernel (4D FAP) computes a ‘4D’ similarity measure from the molecular graphs of an ensemble of conformations which can be incorporated into QSAR models.

The BlueDesc descriptor calculator is a command-line tool that converts an MDL SD file into ARFF and LIBSVM format using CDK and JOELib2 for machine learning and data mining purposes. It computes 174 descriptors taken from both libraries.

MolSig computes molecular graph descriptors that include stereochemistry information.

PaDEL-Descriptor calculates molecular descriptors and fingerprints. It computes 1875 descriptors (1444 1D, 2D descriptors and 431 3D descriptors) and 12 types of fingerprints.

Topological maximum cross correlation descriptors (TMACC) generates 2D autocorrelation descriptors that are low dimensional and interpretable and appropriate for QSAR modeling.

Model Building 
--------------

| Name       | URL                                       | License | Activity | Citation |
|:-----------|:------------------------------------------|:-------:|:--------:|:--------:|
| AZOrange   | <https://github.com/AZCompTox/AZOrange>   |   LGPL  |    C2    |          |
| Bioalerts  | <https://github.com/isidroc/bioalerts>    |   GPL   |    A3    |          |
| camb       | <https://github.com/cambDI>               |   GPL   |    B2    |          |
| eTOXlab    | <https://github.com/manuelpastor/eTOXlab> |   GPL   |    B3    |          |
| Open3DGRID | <http://open3dgrid.sourceforge.net>       |   GPL   |    B1    |          |
| Open3DQSAR | <http://open3dqsar.sourceforge.net>       |   GPL   |    B1    |          |
| QSAR-tools | <https://github.com/dkoes/qsar-tools>     |   BSD   |    A3    |          |

AZOrange is a machine learning package that supports QSAR model building in a full work flow from descriptor computation to automated model building, validation and selection. It promotes model accuracy by using several high performance machine learning algorithms for efficient data set specific selection of the statistical approach.

Bioalerts uses RDKit fingerprints to create models from discrete (e.g., toxic/non-toxic) and continuous data. It includes the capability to visualize problematic functional groups.

Chemistry aware model builder (camb) is an R package for the generation of quantitative models. Its capabilities include descriptor calculation (including 905 two-dimensional and 14 fingerprint type descriptors for small molecules, 13 whole protein sequence descriptors, and 8 types of amino acid descriptors), model generation, ensemble modeling, and visualization.

eTOXLab provides a portable modeling framework embedded in a self-contained virtual machine for easy deployment.

Open3DGrid and Open3DQSAR are a suite of related tools that build 3D QSAR models. Open3DGrid generates molecular interaction fields (MIFs) in a variety of formats, and Open3DQSAR builds predictive models from the MIFs of aligned molecules. Calculations can be visualized in real time in PyMOL.

QSAR-tools is a set of Python scripts that use RDKit to build linear QSAR models from 2D chemical data.

Model Application 
-----------------

| Name     | URL                                                 | License | Activity | Citation |
|:---------|:----------------------------------------------------|:-------:|:--------:|:--------:|
| SMARTCyp | <http://www.farma.ku.dk/smartcyp>                   |   LGPL  |    C1    |          |
| Toxtree  | <http://toxtree.sourceforge.net>                    |   GPL   |    A1    |          |
| UG-RNN   | <http://cdb.ics.uci.edu/cgibin/tools/AquaSolWeb.py> |  Apache |    C1    |          |

SMARTCyp is a QSAR model that predicts the sites of cytochrome P450-mediated metabolism of drug-like molecules directly from the 2D structure of a molecule using fragment-based energy rules.

Toxtree is a Java GUI application for estimating the “toxic hazard” of molecules using a variety of toxicity prediction modules, such as oral toxicity, skin and eye irritation prediction, covalent protein binding and DNA binding, Cytochrome P450-mediated drug metabolism (using SMARTCyp) and more.

UG-RNN/AquaSol is an undirected graph recursive neural network that has been trained to predict aqueous solubility from molecular graphs.

Visualization
-------------


| Name            | URL                                                  | License | Activity | Citation |
|:----------------|:-----------------------------------------------------|:-------:|:--------:|:--------:|
| CheS-Mapper     | <http://ches-mapper.org>                             |   GPL   |    A2    |          |
| DataWarrior     | <http://www.openmolecules.org/datawarrior>           |   GPL   |    A1    |          |
| DecoyFinder     | <http://urvnutrigenomica-ctns.github.io/DecoyFinder> |   GPL   |    A1    |          |
| Scaffold Hunter | <http://scaffoldhunter.sf.net>                       |   GPL   |    A1    |          |
| Synergy Maps    | <https://github.com/richlewis42/synergy-maps>        |   MIT   |    A2    |          |
| VIDEAN          | <https://github.com/jimenamartinez/VIDEAN>           |   BSD   |    A3    |          |
| WCSE            | <http://www.cheminfo.org/wikipedia>                  |   BSD   |    A2    |          |
| WebChemViewer   | <http://sourceforge.net/projects/webchemviewer>      |   BSD   |    C3    |          |

CheS-Mapper (chemical space mapper) . is a 3D-viewer for small compounds in chemical datasets. It embeds a dataset into 3D space by performing dimensionality reduction on the properties of the compounds.

DataWarrior is a data visualization and analysis tool for chemical data with a rich set of available property calculations, similarity metrics, modeling capabilities, and data set representations.

DecoyFinder provides a GUI for selecting a set of decoy compounds from a large library that are appropriate matches to a given set of actives.

Scaffold Hunter provides a Java-based GUI for visualizing the relationship between compounds in a dataset.

Synergy Maps visualizes synergistic activity extracted from screens of drug combinations.

VIDEAN (visual and interactive descriptor analysis) is a visual tool for iteratively choosing a subset of descriptors appropriate for predicting a target property with the aid of statistical methods.

WCSE (Wikipedia chemical structure explorer) runs as a web application and provides a 2D interface for visualizing and searching for 2D molecules.

WebChemViewer is an online viewer for viewing and interacting with lists of compounds and their associated data.

Quantum Chemistry
=================


*Ab initio* Calcuation
----------------------

| Name             | URL                                         | License | Activity | Citation |
|:-----------------|:--------------------------------------------|:-------:|:--------:|:--------:|
| ABINIT           | <http://www.abinit.org>                     |   GPL   |    A1    |          |
| ACES             | <http://www.qtp.ufl.edu/aces>               |   GPL   |    A1    |          |
| BigDFT           | <http://bigdft.org>                         |   GPL   |    A1    |          |
| CP2K             | <http://www.cp2k.org>                       |   GPL   |    A1    |          |
| DACAPO           | <https://wiki.fysik.dtu.dk/dacapo>          |   GPL   |    C1    |          |
| ErgoSCF          | <http://www.ergoscf.org>                    |   GPL   |    C2    |          |
| ERKALE           | <https://github.com/susilehtola/erkale>     |   GPL   |    B2    |          |
| GPAW             | <https://wiki.fysik.dtu.dk/gpaw>            |   GPL   |    A1    |          |
| HORTON           | <http://theochem.github.io/horton>          |   GPL   |    A1    |          |
| JANPA            | <http://janpa.sourceforge.net>              |   BSD   |    A1    |          |
| MPQC             | <http://www.mpqc.org>                       |   LGPL  |    B1    |          |
| NWChem           | <http://www.nwchem-sw.org>                  |   ECL   |    A1    |          |
| Octopus          | <http://www.tddft.org/programs/octopus>     |   GPL   |    A1    |          |
| OpenMX           | <http://www.openmx-square.org>              |   GPL   |    A1    |          |
| Psi4             | <http://www.psicode.org>                    |   GPL   |    A1    |          |
| pyquante         | <http://sourceforge.net/projects/pyquante>  |   BSD   |    A1    |          |
| PySCF            | <https://github.com/sunqm/pyscf>            |   BSD   |    A1    |          |
| QMCPACK          | <http://qmcpack.org>                        |   BSD   |    A2    |          |
| Quantum espresso | <http://www.quantum-espresso.org>           |   GPL   |    A1    |          |
| RMG              | <http://rmgdft.sourceforge.net>             | BSD/GPL |    A1    |          |
| SQ               | <https://sites.google.com/site/siamquantum> |   GPL   |    A2    |          |

ABINIT can calculate the total energy, charge density and electronic structure of molecules and periodic solids with density functional theory (DFT) and Many-Body Perturbation Theory (MBPT), using pseudopotentials and a planewave or wavelet basis. ABINIT also can optimize the geometry, perform molecular dynamics simulations, or generate dynamical matrices, Born effective charges, and dielectric tensors and many more properties.

ACES performs calculations such as single point energy calculations, analytical gradients, and analytical Hessians, and is highly parallelized, including support for GPU computing. A focus of ACES is the use of MBPT and the coupled-cluster approximation to reliable treat electron correlation.

BigDFT performs *ab initio* calculations using Daubechies wavelets and has the capability to use a linear scaling method. Periodic systems, surfaces and isolated systems can be simulated with the proper boundary conditions. It is included as part of ABINIT.

CP2K performs simulations of solid state, liquid, molecular and biological systems. Its particular focus is massively parallel and linear scaling electronic structure methods and state-of-the-art *ab-initio* molecular dynamics (AIMD) simulations. It is optimized for the mixed Gaussian and Plane-Waves method using pseudopotentials and can run on parallel and on GPUs.

Dacapo is a total energy program that uses density functional theory. It can do molecular dynamics/structural relaxation while solving the Schrödinger equations. It has support for parallel execution and is used through the Atomic Simulation Environment (ASE)

ErgoSCM is a quantum chemistry program for large-scale self-consistent field calculations. It performs electronic structure calculations using Hartree-Fock and Kohn-Sham density functional theory and achieves linear scaling for both CPU usage and memory utilization.

ERKALE is designed to compute X-ray properties, such as ground-state electron momentum densities and Compton profiles, and core (x-ray absorption and x-ray Raman scattering) and valence electron excitation spectra of atoms and molecules.

GPAW is a DFT code that uses the projector-augmented wave (PAW) technique and integrates with the atomic simulation environment (ASE) .

HORTON (Helpful Open-source Research TOol for N-fermion systems) has as a primary design goal ease of extensibility for researching new methods in *ab initio* electronic structure theory.

JANPA computes natural atomic orbitals from a reduced one-particle density matrix.

MPQC (massively parallel quantum chemistry program) offers many features including closed shell, unrestricted and general restricted open shell Hartree-Fock energies and gradients, closed shell, unrestricted and general restricted open shell density functional theory energies and gradients, second order open shell perturbation theory and Z-averaged perturbation theory energies.

NWChem provides a full suite of methods for modeling both classical and QM systems. Its capabilities include molecular electronic structure, QM/MM, pseudopotential plane-wave electronic structure, and molecular dynamics and is designed to scale across hundreds of processors.

Octopus pervorms *ab initio* calculations using time-dependent DFT (TDDFT) and pseudopotentials. Included in the project is libxc which is a standalone library of exchange-correlation functionals for DFT (released under the LGPL).

OpenMX (Open source package for material eXplorer) is designed for nano-scale material simulations based on DFT, norm-conserving pseudopotentials, and pseudo-atomic localized basis functions. OpenMX is capable of performing calculations of physical properties such as magnetic, dielectric, and electric transport properties and is optimized for large-scale parallelism.

Psi4 is a suite of *ab initio* quantum chemistry programs that supports a wide range of computations (e.g., Hartree–Fock, MP2, coupled-cluster) and general procedures such as geometry optimization and vibrational frequency analysis with more than 2500 basis functions.

PyQuante is a collection of modules, mostly written in Python, for performing Hartree-Fock and DFT calculations with a focus on providing a well-engineered set of tools. A new version is under development (<https://github.com/rpmuller/pyquante2>).

PySCF is also written primarily in Python and supports several popular methods such as Hartree-Fock, DFT, and MP2. It also has easy of use and extension as primary design goals.

QMCPACK is a many-body *ab initio* quantum Monte Carlo implementation for computing electronic structure properties of molecular, quasi-2D and solid-state systems. The standard file formats utilized for input and output are in XML and HDF5.

QUANTUM ESPRESSO is designed for modeling at the nanoscale using DFT, plane waves, and pseudopotentials and its capabilities include ground-state calculations, structural optimization, transition states and minimum energy paths, *ab initio* molecular dynamics, DFT perturbation theory, spectroscopic properties, and quantum transport.

RMG is a DFT code that uses real space grids to provide high scalability across thousands of processors and GPU acceleration for both structural relaxation and molecular dynamics.

Siam Quantum (SQ) is optimized for parallel computation and its capabilities include the calculation of Hartree-Fock and MP2 energies, minimum energy crossing point calculations, geometry optimization, population analysis, and quantum molecular dynamics.

Helper Applications
-------------------

| Name      | URL                                        | License | Activity | Citation |
|:----------|:-------------------------------------------|:-------:|:--------:|:--------:|
| FragIt    | <https://github.com/FragIt>                |   GPL   |    A2    |          |
| cclib     | <https://github.com/cclib/cclib>           |   LGPL  |    A1    |          |
| GaussSum  | <http://sourceforge.net/projects/gausssum> |   GPL   |    A1    |          |
| Geac      | <https://github.com/LaTruelle/Geac>        |   GPL   |    A3    |          |
| Nancy\_EX | <http://sourceforge.net/projects/nancyex>  |   GPL   |    A2    |          |
| orbkit    | <http://orbkit.sourceforge.net>            |   GPL   |    A2    |          |


FragIt generates fragments of large molecules to use as input files in quantum chemistry programs that support fragment based methods.

cclib provides a consistent interface for parsing and interpreting the results of a number of quantum chemistry packages.

GaussSum uses cclib to extract useful information from the results of quantum chemistry programs (ADF, GAMESS, Gaussian, Jaguar) including monitoring the progress of geometry optimization, the UV/IR/Raman spectra, molecular orbital (MO) levels and MO contributions.

Geac (Gaussian ESI Automated Creator) extracts data from Gaussian log files.

Nancy\_EX post-processes Gaussian output and analyzes excited states including natural transition orbitals, detachment and attachment density matrices, and charge-transfer descriptors.

orbkit is a post-processing tool for the results of quantum chemistry programs. It has native support for a number of programs (MOLPRO, TURBOMOLE, GAMESS-US, PROAIMS/AIMPAC, Gaussian) and additionally interfaces with cclib for additional file format support. It can extract grid-based quantities such as molecular orbitals and electron density, as well as Muliken population charges and other properties.

Visualization
-------------

| Name        | URL                                                 | License | Activity | Citation |
|:------------|:----------------------------------------------------|:-------:|:--------:|:--------:|
| CCP1GUI     | <http://www.scd.stfc.ac.uk/research/app/40501.aspx> |   GPL   |    C3    |          |
| ccwatcher   | <http://sourceforge.net/projects/ccwatcher>         |   GPL   |    B2    |          |
| Gabedit     | <http://gabedit.sourceforge.net>                    |   BSD   |    C1    |          |
| J-ICE       | <http://j-ice.sourceforge.net>                      |   GPL   |    A1    |          |
| QMForge     | <http://qmforge.sourceforge.net>                    |   GPL   |    A1    |          |
| wxMacMolPlt | <http://brettbode.github.io/wxmacmolplt>            |   GPL   |    A1    |          |

CCP1GUI provides a graphical user interface to various computational chemistry codes with an emphasis on integration with the GAMESS-UK quantum chemistry program.

ccwatcher provides a graphical interface for the monitoring of computational chemistry programs.

Gabedit is a graphical user interface to a large number of quantum chemistry packages. It can create input files and graphically visualize calculation results.

J-ICE is a Jmol-based viewer for crystallographic and electronic properties that can be deployed as a Java applet embedded in a web browser.

QMForge provides a graphical user interface for analyzing and visualizing results of quantum chemistry DFT calculations (Gaussian, ADF, GAMESS, Jaguar, NWChem, ORCA, QChem). Analyses include a number of population analyses, Mayer’s bond order, charge decomposition, and fragment analysis.

wxMacMolPlt is a multi-platform GUI for setting up and visualizing input and output files for the GAMESS quantum chemistry software.

Ligand Dynamics and Free Energy Calculations
============================================


Simulation Software
-------------------

| Name             | URL                                          | License | Activity | Citation |
|:-----------------|:---------------------------------------------|:-------:|:--------:|:--------:|
| Campari          | <http://campari.sourceforge.net>             |   GPL   |    B1    |          |
| DL\_POLY Classic | <http://www.ccp5.ac.uk/DL_POLY_CLASSIC/>     |   BSD   |    C3    |          |
| GALAMOST         | <http://galamost.ciac.jl.cn>                 |   GPL   |    A2    |          |
| Gromacs          | <http://www.gromacs.org>                     |   LGPL  |    A1    |          |
| Iphigenie        | <https://sourceforge.net/projects/iphigenie> |   GPL   |    A2    |          |
| LAMMPS           | <http://lammps.sandia.gov>                   |   GPL   |    A1    |          |
| MDynaMix         | <http://www.fos.su.se/~sasha/mdynamix>       |   GPL   |    A1    |          |
| MMTK             | <http://dirac.cnrs-orleans.fr/MMTK>          |  CeCILL |    C1    |          |
| OpenMM           | <https://simtk.org/home/openmm>              | GPL/MIT |    A1    |          |
| ProtoMol         | <http://protomol.sourceforge.net>            |   GPL   |    C1    |          |
| ProtoMS          | <http://www.essexgroup.soton.ac.uk/ProtoMS>  |   GPL   |    A2    |          |
| Sire             | <http://siremol.org>                         |   GPL   |    C2    |          |
| WESTPA           | <https://westpa.github.io/westpa>            |   GPL   |    A2    |          |
| yank             | <http://getyank.org>                         |   LGPL  |    A2    |          |

Campari conducts flexible Monte Carlo sampling of biopolymers in internal coordinate space, with built-in analysis routines to estimate structural properties and support for replica exchange and Wang-Landau sampling.

DL\_POLY Classic is a general purpose molecular dynamics simulation package that can run in parallel and includes a Java graphical user interface.

GALAMOST (GPU accelerated large-scale molecular simulation toolkit) uses GPU computing to perform traditional molecular dynamics with a special focus on polymeric systems at mesoscopic scales.

Gromacs is a complete and well-established package for molecular dynamics simulations that provides high performance on both CPUs and GPUs. It can be used for free energy and QM/MM calculations and includes a comprehensive set of analysis tools.

Iphigenie is a molecular mechanics program that features polarizable force fields, the HADES reaction field, and QM/(P)MM hybrid simulations.

LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is a highly modular classical molecular dynamics simulator that includes a diverse array of energy potentials and integrators.

MDynaMix is a basic general purpose molecular dynamics package.

MMTK (Molecular Modelling Toolkit) is a library written in Python (with some time critical parts written in C) for constructing and simulating molecular systems. Its capabilities include molecular dynamics, energy minimization, and normal mode analysis and it is well-suited for methods development.

OpenMM is a substantial toolkit for high performance molecular dynamics simulations that includes support for GPU acceleration.

ProtoMol , and the associated MDLab Python bindings , provides an object-oriented framework for prototyping algorithms for molecular dynamics simulations and includes an interface to OpenMM.

ProtoMS is a Monte Carlo biomolecular simulation program which can be used to calculate relative and absolute free energies and water placement with the GCMC and JAWS methodologies.

Sire is a collection of modular libraries intended to facilitate fast prototyping and the development of new algorithms for molecular simulation and molecular design. It has apps for system setup, simulation, and analysis.

WESTPA (The Weighted Ensemble Simulation Toolkit with Parallelization and Analysis) is a library for performing weighted ensemble simulations to sample rare events and compute rigorous kinetics.

yank is built off of OpenMM and provides a Python interface for performing alchemical free energy calculations.

Simulation Setup and Analysis 
-----------------------------


| Name       | URL                                               |  License | Activity | Citation |
|:-----------|:--------------------------------------------------|:--------:|:--------:|:--------:|
| AmberTools | <http://ambermd.org>                              |    GPL   |    A1    |          |
| LOOS       | <http://loos.sourceforge.net>                     |    GPL   |    A1    |          |
| lsfitpar   | <http://mackerell.umaryland.edu/~kenno/lsfitpar>  |    GPL   |    A2    |          |
| MDAnalysis | <http://mdanalysis.org>                           |    GPL   |    A1    |          |
| MDTraj     | [mdtraj.org](mdtraj.org)                          |   LGPL   |    A1    |          |
| MEMBPLUGIN | <https://sourceforge.net/projects/membplugin>     |    GPL   |    C1    |          |
| MEPSA      | <http://bioweb.cbm.uam.es/software/MEPSA>         |    GPL   |    A3    |          |
| MSMBuilder | <http://msmbuilder.org>                           |   LGPL   |    A1    |          |
| packmol    | <http://www.ime.unicamp.br/~martinez/packmol>     |    GPL   |    A1    |          |
| PDB2PQR    | <http://www.poissonboltzmann.org>                 |    BSD   |    A1    |          |
| PLUMED     | <http://www.plumed.org>                           |   LGPL   |    A1    |          |
| ProDy      | <http://prody.csb.pitt.edu>                       |    MIT   |    A1    |          |
| Pteros     | <http://pteros.sourceforge.net>                   | Artistic |    B2    |          |
| PyEMMA     | <http://www.emma-project.org>                     |   LGPL   |    A1    |          |
| PyRED      | <http://upjv.q4md-forcefieldtools.org>            |    GPL   |    C1    |          |
| PYTRAJ     | <https://github.com/Amber-MD/pytraj>              |    GPL   |    A1    |          |
| simpletraj | <https://github.com/arose/simpletraj>             |    GPL   |    A2    |          |
| WHAM       | <http://membrane.urmc.rochester.edu/content/wham> |    BSD   |    C1    |          |

AmberTools is an open source component of the non-open source Amber package and provides a large suite of analysis programs. As of Amber15, AmberTools includes the lower performance, but readily extendable, sander molecular dynamics code.

LOOS (Lightweight Object-Oriented Structure library) is a C++ library (with Python bindings) for reading and analyzing molecular dynamics trajectories that also includes a number of standalone programs.

lsfitpar derives bonded parameters for Class I force fields by performing a robust fit to potential energy scans provided by the user.

MDAnalysis is a Python library for reading and analyzing molecular dynamics simulations with some time critical sections written in C.

MDTraj provides high-performance reading, writing, and analysis of molecular dynamics trajectories in a diversity of formats from a Python interface.

MEMBPLUGIN analyzes molecular dynamics simulations of lipid bilayers and is most commonly used as a VMD plugin.

MEPSA (Minimum Energy Pathway Analysis) provides tools for analyzing energy landscapes and pathways.

MSMBuilder is an application and Python library for building Markov models of high-dimensional trajectory data.

packmol is a utility for setting molecular systems by realistically packing molecules to obey a variety of constraints and can create solvent mixtures and lipid bilayers.

PDB2PQR prepares structures for electrostatics calculations by adding hydrogens, calculating sidechain pKa, adding missing heavy atoms, and assigning force field-dependent parameters; users can specify an ambient pH.

PLUMED interfaces with an assortment of molecular dynamics software packages to provide a unified interface for performing free energy calculations using methods such as metadynamics, umbrella sampling and steered MD (Jarzynski).

ProDy is a Python toolkit for analyzing proteins and includes facilities for trajectory analysis and druggability predictions using simulations of molecular probes .

Pteros is a C++ library (with Python bindings) for reading and analyzing molecular dynamics trajectories.

PyEMMA is a Python library for performing kinetic and thermodynamic analyses of molecular dynamics simulations using Markov models.

PyRED generates RESP and ESP charges for the AMBER, CHARMM, OPLS, and Glycam and force fields.

PYTRAJ is a Python interface to the `cpptraj` tool of AmberTools.

simpletraj is a lightweight Python library for parsing molecular dynamics trajectories.

WHAM (Weighted Histogram Analysis Method) calculates the potential of mean force (PMF) from umbrella sampling simulations.

Virtual Screening and Ligand Design
===================================


Ligand-Based
------------

| Name        | URL                                      | License | Activity | Citation |
|:------------|:-----------------------------------------|:-------:|:--------:|:--------:|
| ACPC        | <https://github.com/UnixJunkie/ACPC>     |   BSD   |    B2    |          |
| Align-it    | <http://silicos-it.be>                   |   LGPL  |    C3    |          |
| Open3DALIGN | <http://open3dalign.sourceforge.net>     |   GPL   |    B2    |          |
| PAPER       | <https://simtk.org/home/paper>           |   GPL   |    C2    |          |
| Pharmer     | <http://pharmer.sf.net>                  |   GPL   |    B1    |          |
| Pharmit     | <http://pharmit.sf.net>                  |   GPL   |    A3    |          |
| Shape-it    | <http://silicos-it.be>                   |   LGPL  |    C3    |          |
| USRCAT      | <https://bitbucket.org/aschreyer/usrcat> |   MIT   |    C2    |          |

ACPC (AutoCorrelation of Partial Charges) computes ligand similarity based on a rotation and translation invariant electrostatic descriptor.

Align-it<sup>TM</sup> is a successor of Pharao and aligns and scores 3D representations of molecules based on their pharmacophore features. It includes a plugin for integration with PyMOL.

Open3DALIGN performs unsupervised rigid-body molecular alignment.

PAPER performs GPU accelerated alignment of molecular shapes using Gaussian overlays.

Pharmer uses efficient data structures to rapidly screen large libraries for ligand conformations that match a pharmacophore.

Pharmit is a successor of Pharmer that also incorporates shape matching and energy minimization (if a receptor structure is available) as part of the screen. It is primarily intended to be used as a backend to a web service.

Shape-it<sup>TM</sup> uses Gaussian volumes to align and score molecular shapes.

USRCAT performs “ultra-fast shape recognition” with the addition of pharmacophoric information to rapidly screen compound libraries for similar molecules.


Docking and Scoring
-------------------

| Name                     | URL                                                    | License | Activity | Citation |
|:-------------------------|:-------------------------------------------------------|:-------:|:--------:|:--------:|
| ADplugin                 | <https://github.com/ADplugin>                          |   LGPL  |    A2    |          |
| APBS                     | <http://www.poissonboltzmann.org>                      |   BSD   |    A1    |          |
| AutoDock                 | <http://autodock.scripps.edu>                          |   GPL   |    C1    |          |
| AutoDock Vina            | <http://vina.scripps.edu>                              |  Apache |    C1    |          |
| Clusterizer-DockAccessor | <http://cheminf.com/software/clusterizer_dockaccessor> |   GPL   |    A3    |          |
| DockoMatic               | <https://sourceforge.net/projects/dockomatic>          |   LGPL  |    B1    |          |
| DOVIS                    | <http://bhsai.org/software>                            |   GPL   |    C2    |          |
| idock                    | <https://github.com/HongjianLi/idock>                  |  Apache |    A2    |          |
| MOLA                     | <http://www.esa.ipb.pt/~ruiabreu/mola>                 |   GPL   |    C2    |          |
| NNScore                  | <http://nbcr.ucsd.edu/data/sw/hosted/nnscore>          |   GPL   |    C1    |          |
| Paradocks                | <https://github.com/cbaldauf/paradocks>                |   GPL   |    A2    |          |
| PyRx                     | <http://pyrx.sourceforge.net>                          |   BSD   |    A1    |          |
| rDock                    | <http://rdock.sourceforge.net>                         |   LGPL  |    C1    |          |
| RF-Score                 | <https://github.com/HongjianLi/RF-Score>               |  Apache |    A2    |          |
| smina                    | <https://sourceforge.net/projects/smina>               |   GPL   |    A1    |          |
| VHELIBS                  | <http://urvnutrigenomica-ctns.github.io/VHELIBS>       |   GPL   |    A2    |          |
| VinaLC                   | <http://mvirdb1.llnl.gov/static_catsid/vina>           |  Apache |    C2    |          |
| VinaMPI                  | <http://cmb.ornl.gov/~sek>                             |  Apache |    C2    |          |
| Zodiac                   | <https://sourceforge.net/projects/zodiac-zeden>        |   GPL   |    C1    |          |

ADplugin is a plugin for PyMOL for interfacing with AutoDock and AutoDock Vina.

APBS performs solvation free energy calculations using the Poisson-Boltzmann implicit solvent method.

AutoDock is an automated docking program that uses a physics-based semiempirical scoring function mapped to atom type grids to evaluate poses and a genetic algorithm to explore the conformational space. It includes the ability to incorporate sidechain flexibility and covalent docking.

AutoDock Vina is an entirely separate code base and approach from Autodock that was developed with a focus on runtime performance and ease of system setup. It uses a fully empirical scoring function and an iterated local search global optimizer to produce docked poses. It includes support for multi-threading and flexible sidechains.

Clusterizer-DockAccessor are tools for accessing the quality of docking protocols. It interfaces with a number of open source and free tools.

DockoMatic provides a graphical user interface for setting up and analyzing AutoDock and AutoDock Vina docking jobs, including when run on a cluster. It also includes the ability to run inverse virtual screens (find proteins that bind a given ligand) and support for homology model construction.

DOVIS is an extension of AutoDock 4.0 that provides more efficient parallelization of large virtual screening jobs.

idock is a multi-threaded docking program that includes support for the AutoDock Vina scoring function and a random forest scoring function. I can output per-atom free energy information for hotspot detection.

MOLA is a pre-packaged distribution of AutoDock and AutoDock Vina for deployment on multi-platform computing clusters.

NNScore uses a neural network model to score protein-ligand poses.

Paradocks is a parallelized docking program that includes a number of population-based metaheuristics, such as particle swarm optimization, for exploring the space of potential poses.

PyRx is a visual interface for AutoDock and AutoDock Vina that simplifies setting up and analyzing docking workflows. Its future as an open-source solution is in doubt due to a recent shift to commercialization.

rDock is designed for docking against proteins or nucleic acids and can incorporate user-specified constraints. It uses an empirical scoring function that includes solvent accessible surface area terms. A combination of genetic algorithms, Monte Carlo, and simplex minimization is used to explore the conformational space. Distinct scoring functions are provided for docking to proteins and nucleic acids.

RF-Score uses a random forest classifier to score protein-ligand poses.

smina is a fork of AutoDock Vina designed to better support energy minimization and custom scoring function development (scoring function terms and atom type properties can be specified using a run-time configuration file). It also simplifies the process of setting up a docking run with flexible sidechains.

VHELIBS (Validation HElper for LIgands and Binding Sites) assists the non-crystallographer in validating ligand geometries with respect to electron density maps.

VinaLC is a fork of AutoDock Vina designed to run on a cluster of multiprocessor machines.

VinaMPI is a wrapper for AutoDock Vina that uses OpenMPI to run large-scale virtual screens on a computing cluster.

Zodiac is a visual interface for structure-based drug design that includes support for haptic feedback.

Pocket Detection
----------------

| Name         | URL                                                          | License | Activity | Citation |
|:-------------|:-------------------------------------------------------------|:-------:|:--------:|:--------:|
| eFindSite    | <http://brylinski.cct.lsu.edu/efindsite>                     |   GPL   |    C2    |          |
| fpocket      | <http://fpocket.sourceforge.net>                             |   GPL   |    C1    |          |
| KVFinder     | <http://lnbio.cnpem.br/facilities/bioinformatics/software-2> |   GPL   |    B1    |          |
| mcvol        | <http://www.bisb.uni-bayreuth.de/data/mcvol/mcvol.html>      |   GPL   |    C2    |          |
| PAPCA        | <https://sourceforge.net/projects/papca>                     |   BSD   |    C2    |          |
| PCS          | <https://sourceforge.net/projects/cavity-search>             |   GPL   |    C2    |          |
| PocketPicker | <http://gecco.org.chemie.uni-frankfurt.de/pocketpicker>      |   BSD   |    C1    |          |
| POVME        | <https://sourceforge.net/projects/povme>                     |   GPL   |    C1    |          |


eFindSite using homology modeling and machine learning predicts ligand binding sites in a protein structure.

fpocket detects and delineates protein cavities using Voronoi tessellation and is able to process molecular dynamics simulations.

KVFinder is a PyMOL plugin for identifying and characterizing pockets.

mcvol calculates protein volumes and identifying cavities using a Monte Carlo algorithm.

PAPCA (PocketAnalyzerPCA) is a pocket detection utility designed to analyze ensembles of protein conformations.

PCS (Pocket Cavity Search) measures the volume of internal cavities and evaluates the environment of ionizable residues.

PocketPicker is a PyMOL plugin that automatically identifies potential ligand binding sites using a grid-based shape descriptor.

POVME (POcket Volume MEasurer) is a tool for measuring and characterizing pocket volumes that includes a graphical user interface.



Ligand Design
-------------

| Name          | URL                                              | License | Activity | Citation |
|:--------------|:-------------------------------------------------|:-------:|:--------:|:--------:|
| AutoClickChem | <https://sourceforge.net/projects/autoclickchem> |   GPL   |    C2    |          |
| AutoGrow      | <http://autogrow.ucsd.edu>                       |   GPL   |    A1    |          |
| igrow         | <https://github.com/HongjianLi/igrow>            |  Apache |    A3    |          |
| OpenGrowth    | <http://opengrowth.sf.net>                       |   GPL   |    A1    |          |

AutoClickChem performs *in silico* click chemistry to generate large libraries of synthetically accessible compounds.

AutoGrow uses a genetic algorithm to explore the space of reactants and reactions accessible via AutoClickChem and identifies compounds that dock well using AutoDock Vina.

igrow, like AutoGrow, uses a genetic algorithm but transforms ligands using branch exchange and uses idock as the underlying docking evaluation protocol.

OpenGrowth assembles candidate ligands by connecting small organic fragments in the active site of proteins. It includes a graphical user interface.

