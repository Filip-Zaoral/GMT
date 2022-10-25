
# Old algorithm to identify surface loops ("islands") of the STL geometry:    

E = gmsh.model.mesh.getElementsByType(2,-1,0)[1]
E = np.reshape(E,(-1,3),'C')
nI = 1
Ie = []; Ie.append([0])
nea = 0
ea = []
for i in range(nE)[1:]:                                                        # Loops over all triangles in the STL geometry in order to find the total number of surface loops ("islands") in the model.
    for j in range(nI):
        eInE = any(np.isin(E[Ie[j]],E[i],True,False).any(0))
        if eInE:
            nea += 1
            ea.append(j)
            Ie[j].extend([i])
    ea.sort(reverse = True)
    if nea == 0:
        nI += 1
        Ie.append([i])
    else:
        for j in range(nea)[:-1]:
            Ie[ea[-1]] += Ie[ea[j]][:-1]
            Ie.pop(ea[j])
        nI -= nea - 1
    nea = 0
    ea = []
del E
