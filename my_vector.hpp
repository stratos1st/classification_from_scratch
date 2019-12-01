#ifndef MY_VECTOR
#define MY_VECTOR

#include <iostream>

class my_vector{
public:
  double *coordinates;
  unsigned int id;
  unsigned int dim;
  static double vector_tol;

  my_vector(unsigned int dimentions);
  ~my_vector();
  my_vector& operator=(const my_vector &other);//copies id
  my_vector(const my_vector &p2);//copies id

  unsigned int get_dimentions() const;
  void print_vec(unsigned int until=0) const;//prints x1,x1,...,xuntil. if until==0 prints all
  bool operator==(const my_vector &other);
};

// custom specialization of std::hash
namespace std{
  template<>
  struct hash<my_vector*>{
    size_t operator()(const my_vector* const s) const noexcept{
      size_t h1 = hash<unsigned int>{}(s->id);
      return h1;
    }
  };
}

#endif
