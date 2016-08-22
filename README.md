# fib-representation
Simulation code for efficient FIB representations

In order to build the basic program `reduce_fib`, you have to have boost libraries installed; e.g., on an Ubuntu system run
```
sudo apt-get install libboost-dev libboost-random-dev libboost-system-dev libboost-program-options-dev libboost-filesystem-dev libboost-date-time-dev
```

Then, to produce the tables similar to ICNP 2016 submission, run
```
python3 run_icnp2016.py -i data -t 8
```
where `data` is the directory where input FIBs are stored (in *.txt format) and `8` is the (max) number of threads to use for experiments.


