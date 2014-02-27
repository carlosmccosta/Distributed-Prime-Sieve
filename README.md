# [Distributed prime sieve](http://carlosmccosta.github.io/Distributed-Prime-Sieve/)

## Overview
This project provides several efficient implementations of ditributed prime sieves and has the following associated paper:

[Distributed prime sieve in heterogeneous computer clusters](https://www.researchgate.net/publication/257213711_Distributed_prime_sieve_in_heterogeneous_computer_clusters)


**Abstract:**
Prime numbers play a pivotal role in current encryption algorithms  and given the rise of 
cloud computing, the need for larger primes never been so high. This increase in available computation 
power can be used to either try to break the encryption or to strength it by finding larger primes. With 
this in mind, this paper provides an analysis of different sieves  implementations that can be used  to 
generate primes up to 2^64. It starts by analyzing cache friendly sequential sieves with wheel factorization, then expands to multicore architectures and ends with a cache friendly segmented hybrid implementation of a distributed prime sieve, designed to  efficiently use all the available computation resources of  heterogeneous computer clusters with variable workload  and to scale very well to any cluster 
size.



## Results

![Fig. 1 - Global performance comparison](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/GlobalPerformance.png)
Fig. 1 - Global performance comparison


![Fig. 2 - Real speedup](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/RealSpeedup.png)
Fig. 2 - Real speedup


![Fig. 3 - Efficiency](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/Efficiency.png)
Fig. 3 - Efficiency


![Fig. 4 - Scalability_shared_memory_algorithm_efficiency](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/Scalability_shared_memory_algorithm_efficiency.png)
Fig. 4 - Scalability of shared memory algorithm (efficiency vs number of cores)


![Fig. 5 - Scalability_shared_memory_algorithm_time](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/Scalability_shared_memory_algorithm_time.png)
Fig. 5 - Scalability shared memory algorithm (running time vs number of cores)


![Fig. 6 - Scalability_distributed_algorithm_efficiency](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/Scalability_distributed_algorithm_efficiency.png)
Fig. 6 - Scalability distributed algorithm (efficiency vs number of nodes - each node is a quadcore processor with different processing capabilities)


![Fig. 7 - Scalability_distributed_algorithm_time](https://raw.github.com/carlosmccosta/Distributed-Prime-Sieve/master/DistributedPrimeSieve/docs/Scalability_shared_memory_algorithm_time.png)
Fig. 7 - Scalability distributed algorithm (running time vs number of nodes - each node is a quadcore processor with different processing capabilities)


## Usage
A brief usage guide is available [here](https://github.com/carlosmccosta/Distributed-Prime-Sieve/blob/master/DistributedPrimeSieve/docs/usage.txt)
