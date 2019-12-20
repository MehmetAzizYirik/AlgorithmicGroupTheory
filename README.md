# AlgorithmicGroupTheory

Copyright 2019 Mehmet Aziz Yirik

## Introduction

Algorithmic group theory functions shared in this repository are especially developed for the implementations in mathematical chemistry. However, the functions can be used also in different fields of applied group theory. The current algorithms are for the molecular structure generation problem.

## Method

This is an ongoing project for molecular structure generation. Already developed functions are regularly shared on this repository. The prototype algorithm is the demo function of the jar file. The theoretical basis of the algorithm can be found in [1]. The outlines of the functions can be found also in [1,2].

## Download jar File

The current jar file  provides only the connectivity matrix generator. The file can be downloaded from [here](https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory/blob/master/cangen.jar)

## Download Source Code

It is assumed that users have git on their system and have initialised their local directory. For more information [set-up-git](https://help.github.com/articles/set-up-git/ )

To download AlgorithmicGroupTheory source code:

```
$ git clone https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory.git
```
## Compiling

To compile AlgorithmicGroupTheory, Apache Maven and Java 1.8 (or later) are required.
```
AlgorithmicGroupTheory/$ mvn package
```
This command will create jar file named specifically as "jar-with-dependencies" under target folder.

## Usage

AlgorithmicGroupTheory.jar can be run from command line with the specified arguments. An example command is given below.

The definitions of the arguments are given below:

```
usage: java -jar AlgorithmicGroupTheory.jar -f <arg> [-v] -d <arg>

Generates adjacency matrices for a given molecular formula. The input is 
a molecular formula string. For example 'C2H2O'. Besides this formula, the
directory is needed to be specified for the output file.

 -f,--molecularFormula <arg>   Molecular formula as a string (required)
 
 -v,--verbose                  Print messages about matrix generation
 
 -d,--filedir <arg>            Creates and store the output txt file 
                               in the directory (required)

Please report issues at https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory
```

```
java -jar AlgorithmicGroupTheory.jar -f C2H2O -v -d C:\Users\UserName\Desktop\AlgorithmicGroupTheory/output.txt
```

## Running the Tests

For the Generator class, a test class called Test-Generator is built. This test class includes the tests of the main functions. The outputs of the the functions are tested based on the size ( or the length) of the expected output files. 

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/MehmetAzizYirik/AlgorithmicGroupTheory/blob/master/LICENSE) file for details

## Authors

 - Mehmet Aziz Yirik - [MehmetAzizYirik](https://github.com/MehmetAzizYirik)
 
## Acknowledgements
![YourKit](https://camo.githubusercontent.com/97fa03cac759a772255b93c64ab1c9f76a103681/68747470733a2f2f7777772e796f75726b69742e636f6d2f696d616765732f796b6c6f676f2e706e67)

The developer uses YourKit to profile and optimise code.

YourKit supports open source projects with its full-featured Java Profiler. YourKit, LLC is the creator of YourKit Java Profiler and YourKit .NET Profiler, innovative and intelligent tools for profiling Java and .NET applications.

![cdk](https://github.com/MehmetAzizYirik/HMD/blob/master/cdk.png)

This project relies on the Chemistry Development Project (CDK), hosted under [CDK GitHub](http://cdk.github.io/). Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is available though our [Github site](http://cdk.github.io/cdk/).

## References

1- Kerber, A., Laue, R., Meringer, M., Rücker, C. and Schymanski, E., 2013. Mathematical chemistry and chemoinformatics: structure generation, elucidation and quantitative structure-property relationships. Walter de Gruyter.

2- Grund, R. and Müller, R., 1995. Konstruktion molekularer Graphen mit gegebenen Hybridisierungen und überlappungsfreien Fragmenten. Lehrstuhl II für Mathematik.

