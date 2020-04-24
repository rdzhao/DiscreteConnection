#ifndef __UTILITY__
#define __UTILITY__

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

const double PI = 3.1415926535897932384626433;
const double EPS = 0.0000000001;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef Kernel::Vector_3                                        Vector_3;

typedef CGAL::Polyhedron_3<Kernel>                              Polyhedron;
typedef Polyhedron::Vertex_iterator                             Vertex_iterator;
typedef Polyhedron::Edge_iterator                               Edge_iterator;
typedef Polyhedron::Halfedge_iterator                           Halfedge_iterator;
typedef Polyhedron::Facet_iterator                              Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator           Halfedge_around_vertex_circulator;

typedef Eigen::SparseMatrix<double>                             SparseMatrix;
typedef Eigen::VectorXd                                         VectorXd;
typedef Eigen::Triplet<double>                                  Triplet;

#endif