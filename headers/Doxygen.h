//
//  documentation mainpage
//
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __DOXYGEN__
#define __DOXYGEN__

/**
 @mainpage A Simple object oriented finite element code
 
 Object oriented programming combines well with finite element computations :
 elements are objects, nodes are objects, matrices are objects etc.
 
 The proper definition of classes, however, is not straightforward. In this code we separate the geometric map (\ref geometrymap) from the definition of the approximation space (\ref shape). This makes the code inherently extensible: other approximation spaces can be added without changing the geometric map. The geometric map can be extended without affecting the definition of the approximation space.
 
 
 
 Another separate module is the definition of the variational statement. Finite elements applied to partial differential equations invariably apply a Galerkin approximation to a variational statement associated with a (system of) partial differential equation(s). Keeping the classes that compute the integration point contributions seperate allow to use the same class to the Poisson equation and/or elasticity equation.
 */
#endif
