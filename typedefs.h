#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <iostream>
#include <fstream>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <stdio.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<4> LCC_4;
typedef LCC_4::Dart_handle Dart_handle;
typedef LCC_4::Dart_const_handle Dart_const_handle;
typedef LCC_4::Point Point;
typedef LCC_4::Vector Vector;
typedef LCC_4::FT FT;

#endif // TYPEDEFS_H
