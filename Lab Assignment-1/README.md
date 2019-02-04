# Lab Assignment 1

## Sequential Program
### Commands
Use this command to run the sequential code </br>
`g++ main_sequential.c lab1_io.c lab1_sequential.cpp`</br>
`./a.out 1 data.txt data2.txt data3.txt` </br>

### Files Description
`data.txt` contains the input points of the data. `data2.txt` will have the centroid indexes of all the input points. `data3.txt` will have all the computed centroids from each iteration of the algorithm.</br></br>

Use the following command to run the openmp version of the code </br>
`clang++ -Xpreprocessor -fopenmp main.cpp -o a.out -lomp`</br></br>
Click [here](https://stackoverflow.com/questions/39979836/using-openmp-with-c11-on-mac-os) for more help.</br>


Click [here](https://docs.google.com/document/d/1WxGl2QPuQrBRN0awb77xorafp1oeeL-lia-EWp_pA3k/edit) for the report.</br>