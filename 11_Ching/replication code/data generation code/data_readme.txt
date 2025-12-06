This folder contains the files necessary to create the estimation data from the raw data.  All the code to create the data is written in R and C++ and the code can be run directly from the R prompt.

Please follow the steps below to create the estimation data:

1. Obtain the IRI Academic Dataset, following the instructions in Bronnenberg, Kruger and Mela (2008).  Once you have obtained the data, extract it into a folder.

2. Create a folder where the output data will be stored.

3. Open the file, data_master.r, and set the variable rawdatapath on line 20 to point to the folder containing the IRI data.  Set the variable outputpath on line 21 to point to the folder you created to contain the output data.

4. Run data_master.r. It will take some time for the files to process.