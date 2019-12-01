#ifndef MY_CURVE
#define MY_CURVE

#include "my_vector.hpp"

class my_curve{
 public:
  my_vector** vectors;
  unsigned int numofvectors;
  unsigned int vectordimentions;
  unsigned int id;
  static double curve_tol;

  my_curve(unsigned int points,unsigned int dimentions=2);//how many points the curve will have and their dimentions
  ~my_curve();

  my_vector& get_vector(unsigned int);
  unsigned int get_dimentions() const;//returns dimentions
  void print_vec(unsigned int until=0);//prints x1,x1,...,xuntil. if until==0 prints all
  bool operator==(const my_curve &other);
  my_curve& operator=(const my_curve &other);
  my_curve(const my_curve &p2);
};

// custom specialization of std::hash
namespace std{
  template<>
  struct hash<my_curve*>{
    size_t operator()(const my_curve* const s) const noexcept{
      size_t h1 = hash<unsigned int>{}(s->id);
      return h1;
    }
  };
}

#endif
