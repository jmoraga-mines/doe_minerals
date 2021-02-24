#### DOE Project: DE-EE0008760

### USE

python doe_ann3.py -a -d brady_samples_19x3d -l tmp/doe_19x3d200b.l -m tmp/doe19x3d200b.h5 -k 19 -b 200 -e 200 -c 3 -r

(-r : reset) Creates from scratch

(-a) Data augmentation

(-d) Samples file directory 
For use with the AI
“-d brady_samples_19x3d”

(-l) Labels file output directory
“-l tmp/doe_19x3d200b.l”

This will contain the network architecture and all weights
(-m) Model file output directory
“-m tmp/doe19x3d200b.h5” 

This will contain the network architecture and all weights

(-k) The kernel to use is 19x19

(-b: batch size) it will run batches of size 200

(-e: epochs) and will run for 200 full iterations

(-c) The Network will use 3 data bands

