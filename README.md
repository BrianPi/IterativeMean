# IterativeMean

A class from which to derive iterative means. There are three means included as examples, a normal arithmetic-geometric mean (accepting only positive values), a custom mean 'Clocks-at-Sea', described below, and an arithmetic-geometric mean extended to accept negative values. I used gprof to make sure each mean behaved efficiently when compiled with -O3.

The 'Clocks-at-Sea' mean works by modelling having a certain number of 'clocks' which have different relative displacements ('times') each 'day'. After each 'day', the clock that's the farthest from the collective mean is reset to the collective mean. It involves carrying a state for each 'clock' over a number of iterations, and is therefore an optimization challenge.
