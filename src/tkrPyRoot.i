%module tkrPyRoot
%{
#include <vector>
#include "../src/tkrPyRoot/tkrPyRoot.cxx"
%}
%include "carrays.i"
%array_class(double, doubleArray);

%include "std_string.i"
%include "std_vector.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;

%include ../src/tkrPyRoot/tkrPyRoot.cxx
