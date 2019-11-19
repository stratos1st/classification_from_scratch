#ifndef INIT_K
#define INIT_K

#include <iostream>
#include <list>

#include "my_vector.hpp"
#include "util.hpp"

using namespace std;


list<my_vector>* initialization2(list<my_vector>*,unsigned int,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
void updatedist(my_vector*,list<pair<double,my_vector*>> &,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
my_vector* newcentroid(list<pair<double,my_vector*>> &);
my_vector* rangebinarysearch(double ,pair<double,my_vector*>* ,int);


#endif
