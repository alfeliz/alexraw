# alexraw
Octave script to extract electrical data from ALEX scope. Initial works from 2015.

- alexraw04.m Octave script to extract electrical data from scope channels.
- baseline.m Function used by "alexraw04.m" to remove the baseline on a signal.
- chan.m Function used by "alexraw04.m" to extract data from scope channel files.
- display_rounded_matrix.m Function used by "alexraw04.m" to save data in a readable and nice way.
- supsmu.m Function used by "alexraw04.m" to smooth correctly data (It uses a filter) 
