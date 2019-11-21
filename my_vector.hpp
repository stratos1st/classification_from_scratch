#ifndef MY_VECTOR
#define MY_VECTOR

#include <iostream>

class my_vector{
public:
  double *coordinates;
  unsigned int id;
  unsigned int dim;

  my_vector(unsigned int dimentions);
  ~my_vector();
  my_vector& operator=(const my_vector &other);//copies id
  my_vector(const my_vector &p2);//copies id

  unsigned int get_dimentions() const;
  void print_vec(unsigned int until=0) const;//prints x1,x1,...,xuntil. if until==0 prints all
  friend bool operator==(const my_vector &other,const my_vector &other2);
};

// //structs i need to compaire and hash *my_vector. used in unordered_set and find_if
// struct PointedMy_vectorEq{
//   my_vector *var;
//
//   PointedMy_vectorEq(my_vector* a=NULL):var(a){}
//
//   bool operator () (my_vector const * cls) const{
//     return *cls==*var;
//   }
//
//   bool operator () (my_vector const * lhs, my_vector const * rhs ) const{
//     return *lhs==*rhs;
//   }
// };
//
// struct PointedMy_vectorHash{
//   size_t operator() (my_vector* const& c) const {
//     std::hash<my_vector*> my_vector_hasher;
//     return my_vector_hasher(c);
//   }
// };


#endif
