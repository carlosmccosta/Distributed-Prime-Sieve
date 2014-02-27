# Distributed prime sieve

## Overview
This project provides an efficient implementation of a ditributed prime sieve and has the following associated paper:

[Distributed prime sieve in heterogeneous computer clusters](https://github.com/carlosmccosta/Distributed-Prime-Sieve/blob/master/DistributedPrimeSieve-Report/DistributedPrimeSieve.pdf?raw=true)


**Abstract:**
Prime numbers play a pivotal role in current encryption algorithms  and given the rise of 
cloud computing, the need for larger primes never been so high. This increase in available computation 
power can be used to either try to break the encryption or to strength it by finding larger primes. With 
this in mind, this paper provides an analysis of different sieves  implementations that can be used  to 
generate primes up to 2^64. It starts by analyzing cache friendly sequential sieves with wheel factorization, then expands to multicore architectures and ends with a cache friendly segmented hybrid implementation of a distributed prime sieve, designed to  efficiently use all the available computation resources of  heterogeneous computer clusters with variable workload  and to scale very well to any cluster 
size.
