# TODO for extracting gaze traces
Your task is to modify the current file in such a way that it will also extract 
gaze traces from the SpikeGLX analog files and save them in the neurodata without borders file. Lines x-y of the current file show how to open the analog files. In zero based indexing, the gaze traces are in channel 0 and 1. In the same way as the script separates 
the spiking data into trials and further factors the data based on weather the laser was 
on or not, the gaze data needs to be separated based on those criteria.
