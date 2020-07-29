The current folder contains scale free network constructed according to the algorithm outlined in Rozenfeld, 2002, Scale free networks on lattices.

Constraints used: minimum no:of nodes = 7, maximum no:of nodes = 117

The Barabasi-Albert algorithm, without modifications, resulted in very sparse networks, with a large number of nodes with no connections. Also, the hubs tended to be towards the top left corner of the network, the way I was constructing the network.

In order to keep the number of connections comparable to that of the other two cases, the variance of the Gaussian has been changed.

2.0 - 6 to 116 - Gaussian, d
2.5 - 7 to 117 - Gaussian, 1.5*d
3.0 - 8 to 118 - Gaussian, 1.5*d
3.5 - 9 to 119 - Gaussian, 1.5*d
4.0 - 10to 120 - Gaussian, 1.5*d



