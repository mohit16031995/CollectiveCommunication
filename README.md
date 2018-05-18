We describe our work on improving the performance of collective communication
operations. For each collective operation, we can use multiple algorithms depending
on the message size, with the goal of minimizing latency for short messages and
minimizing bandwidth use for long messages. We will implement algorithms for all
MPI (Message Passing Interface) collective operations such as Broadcat, Reduce ,
Scan, and also neighborhood collectives.

In this thesis, we propose, develop and implement a tool for choosing an appropriate
algorithm for a particular collective operation under particular parameters. Also
the analysis could lead to results which might help to develop other collective operations
or algorithms.
