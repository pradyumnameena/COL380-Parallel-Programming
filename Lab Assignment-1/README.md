# Lab Assignment 1

## Sequential Program
### Commands

Use this command to run the sequential code </br>
`clang++ -Xpreprocessor -fopenmp main_X.c lab1_io.c lab1_X.cpp -o a.out -lomp`</br>

`./a.out number_of_clusters number_of_threads data.txt data2.txt data3.txt` </br>

#### Avoid number of threads if using sequential code

### Files Description
`data.txt` contains the input points of the data. `data2.txt` will have the centroid indexes of all the input points. `data3.txt` will have all the computed centroids from each iteration of the algorithm.</br></br>

Click [here](https://docs.google.com/document/d/1WxGl2QPuQrBRN0awb77xorafp1oeeL-lia-EWp_pA3k/edit) for the report.</br>