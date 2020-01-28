class ABString:
    """ The class implements encoding scheme for alpha-beta strings
    representing Slater determinants """

    def __init__(self, nel, norb, verbose):
        self.nel = nel
        self.norb = norb
        self._calculate_vertex_weights(verbose)
        self._calculate_edge_weights(verbose)

    def _pair2idx(self, p):
        return p[0] * (self.nel + 1) + p[1]

    def _idx2pair(self, i):
        p1 = i % ( self.nel + 1) 
        p0 = (i - p1) // (self.nel + 1)
        return (p0, p1)

    def _get_parents(self, pvertex):
        # Input : pvertex - vertex as a tuple
        # Output : parents - array of encoded parent vertices
        k, m = pvertex
        parents = []
        if m > 0 and k > 0 and m != k:
            parents.extend([self._pair2idx((k - 1, m - 1)), self._pair2idx((k - 1, m))])
        elif m == 0 and  k > 0:
            parents.append(self._pair2idx((k - 1, m)))
        elif m == k and m != 0:
            parents.append(self._pair2idx((k - 1, m - 1)))

        return parents

    def _vertex_weight(self, vertex, vweights):
        # Recursive function that calculates the weights
        # of vertices with memoization;
        # The function does not return weight - it updates
        # vweights instead

        if vertex not in vweights:
           k, m = self._idx2pair(vertex)
           # Any m > 0 implies that there is a diagonal parent

           vweights[vertex] = 0
           for p in self._get_parents((k, m)):
               self._vertex_weight(p, vweights)
               vweights[vertex] += vweights[p]


    def _calculate_vertex_weights(self, verbose=False):
        vertex_labels = []
        for m in range(0, self.nel + 1):
            for k in range(m, self.norb - self.nel + m + 1):
                vertex_labels.append(self._pair2idx((k, m)))

        self.vertex_weights = {0 : 1} # important
        for l in vertex_labels:
            #print("Before {} --> {}".format(l, vertex_weights))
            self._vertex_weight(l, self.vertex_weights)

        if verbose:
            print("Vertex weights for the graph ({}, {})".format(self.nel, self.norb))
            for i, v in self.vertex_weights.items():
                print(" {} : {} ".format(i, v))

    def _calculate_edge_weights(self, verbose=False):
        self.edge_weights = {}
        for v, w in self.vertex_weights.items():
            k, m = self._idx2pair(v)
            k, m = k - 1, m - 1
            if k >= 0 and m >= 0:
                self.edge_weights[v] = w - self.vertex_weights[self._pair2idx((k, m))]
            else:
                self.edge_weights[v] = 0

        if verbose:
            print("Edge weights for the graph ({}, {})".format(self.nel, self.norb))
            for i, v in self.edge_weights.items():
                print(" -{} : {} ".format(i, v))

    def str2address(self, s):
        # Orbital numbers start with 0
        # Classical map reduce but will do in a more simple 
        # way
        address = 1
        for pos, o in enumerate(s):
            address += self.edge_weights[self._pair2idx((o + 1, pos + 1))]

        return address

    def _a2s(self, pvertex, remaining_weight):
        # Input : pvertex - vertex as a (k, m) pair
        #         remaining_weight - index of a str/substr encoded
        #                            by the graph defined by pvertex
        # Output : det - list of orbital indeces
        k, m = pvertex
        ivertex = self._pair2idx(pvertex) # vertex code
        parents = self._get_parents(pvertex)
        orbitals = []
        for p in parents:
            k1, m1 = self._idx2pair(p)
            diagonal = ((k - k1) == 1 and (m - m1) == 1)
            if diagonal: 
                new_remaining_weight = remaining_weight - self.edge_weights[ivertex]
                if  new_remaining_weight <= self.vertex_weights[p] and new_remaining_weight > 0: # it can never be 0 if we are doing things right
                    print("Diagonal call with {} and weight {}".format((k1, m1), new_remaining_weight))
                    orbitals_ = self._a2s((k1, m1), new_remaining_weight)
                    orbitals.extend(orbitals_)
                    orbitals.append(k - 1) # should it be -1?
            else:
                if remaining_weight <= self.vertex_weights[p]:
                    print("Vertical(up) call with {} and weight {}".format((k1, m1), remaining_weight))
                    orbitals_ = self._a2s((k1, m1), remaining_weight)
                    orbitals.extend(orbitals_)

        # Need to think about how to handle exceptions properly in this function
        #if len(parents) > 0 and len(orbitals) == 0:
        #    assert False, "Decoding exception (we should never get here!)"

        return orbitals


    def address2str(self, a):
        # the request will be handled by a private recursive method
        return self._a2s((self.norb, self.nel), a)








