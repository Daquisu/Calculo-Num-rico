# Number classification with Non-Negative Matrix Factorization (NNMF).

This project uses the concept of Non-Negative Matrix Factorization to learn a base of images to represent each MNIST number.
Using NNMF, we classify each new digit according to the closest digit created using a combinations of base images.

Thus we are able to interpret the results as well as visualize the closest digit representation.

# How to use

We compiled using gcc -std=c++11

If desired, one can run the program using Code::Blocks by opening the file ep1.cbp and, in Setting options, checking the box "Have g++ follow the C++11 ISO C++ language standard [-std=c++11]" and then compile.

The entire source code used is in the file main.cpp
