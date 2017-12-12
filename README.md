# MSDTracking

To install the scripts on your computer, navigate to the MSDTracking folder and run
```
cp *.R *.py *.sh /usr/local/bin/
```

To use, simply run 
```
MSD.sh Tracking_file.txt Output_prefix Frame_duration
```
All the files/plots will be created in your current directory.

# Combining tracking files

To combine several tracking files into one, run the following command with all the files you need to merge
```
combineFiles.sh Tracking_file.txt name_of_output_file.txt
```

Running the command in sequence with different tracking files will append trajectories and update the output file.
