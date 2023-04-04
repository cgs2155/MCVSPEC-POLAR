-------MCVSPEC User Manual---------

MCVSPEC has five input parameters: B,
MDOT0, M, ABUND, and SIGMA_S. This version of 
this model primarily serves to find the mass of
IPs and Polars. 

The first parameter, B, is the B-field in MG.

The second parameter, MDOT0, is the specific accretion
rate of the compact object.

The third parameter, M, is the mass of the white dwarf
in solar masses.

The fourth parameter, ABUND, is the abundance of
the white dwarf in terms of solar abundance 
and uses AtomDB for the modelling of spectral
lines.

The fifth parameter, SIGMA_S, is the ratio of
electron pressure to ion pressure at the shock height.
It is taken from Saxton, 2005, and it accounts for 
two-temperature effects in polars. The higher 
B of polars leads to significantly more efficient
bremmstraulung energy loss than electrons can regain
through collisions, leading to unequal electron and ion
temperatures. It should be used for sources with high B
fields and will decrease mass overestimation.
 