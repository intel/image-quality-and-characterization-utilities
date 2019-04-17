# Image quality and characterization utilities

## Introduction

Image quality and characterization utilities is a minimum overhead utility 
library for image quality and raw image characterization. 

Its purpose is to provide:
-Commonly used tools like color error calculations
-Chart detectors
-Image pre-processing
-Commonly known RAW image characterization methods

## License

BSD 3-Clause License

Copyright (c) 2019, Intel Corporation
All rights reserved.

## Dependencies

Image quality utilities uses Catch2 for tests and Eigen for certain math
operations. For new user convenience (and as the library is targeted for imaging 
experts) source code of both Eigen and Catch are located in the repository.

* Eigen - https://github.com/libigl/eigen
* Catch - https://github.com/catchorg/Catch2

## Getting started

Building Image quality and characteriation utilities requires cmake 3.10 or later

If you are using Visual Studio 2013, 2015, or 2017, you can create the dev env
by running corresponding vs20xx_64bit.cmd located in the main folder. This is a
convenience utility for people who are not familiar to cmake, not a requirement.
Alternative development tools can be used.

Best way to get acquinted with the code is to follow behavioral tests in 
.\image-quality-utilities\tests\ and it's subfolders, as those work also as an
specification to the libraries.

Tests can be identified as specs_*.cpp

## Intended use

As you can see the repo contain header only libraries, e.g. MacbethDetector.hpp
The intended use is to include the needed header only library, if it needs other
functionalities from the repo it will include them. The idea is inspired by Eigen
library (https://github.com/libigl/eigen).


