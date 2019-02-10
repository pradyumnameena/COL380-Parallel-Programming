Sometimes the code returns a few errors related to malloc. Kindly instead of categorizing it as invalid output run it multiple times or see the text file which contains the mapping of the point with cluster. In case of large datasets pthread tends to give more errors than sequential or openmp.

Also if running on mac use the following command
clang++ -Xpreprocessor -fopenmp main_sequential.c lab1_io.c lab1_sequential.cpp -o a.out -lomp