#include <gtest/gtest.h>
#include "my_vector.hpp"
#include "my_curve.hpp"
#include "my_vector.cpp"
#include "my_curve.cpp"
#include "clustering_funcs.hpp"
#include "clustering_funcs.cpp"
#include "util.hpp"
#include "util.cpp"

using namespace std;




//-------------Binsearch--------------------//
TEST(RangeBinarySearch, left){
  pair<double,my_curve*>* prob_array = new pair<double,my_curve*> [100];
  for(unsigned int i = 0 ;i<100;i++){
    prob_array[i].second = new my_curve(10);
    prob_array[i].first = 2*i;
  }
  ASSERT_EQ(prob_array[49].second,rangebinarysearch(97.5,prob_array,99));
}

TEST(RangeBinarySearch, right){
  pair<double,my_curve*>* prob_array = new pair<double,my_curve*> [100];
  for(unsigned int i = 0 ;i<100;i++){
    prob_array[i].second = new my_curve(10);
    prob_array[i].first = 2*i;
  }
  ASSERT_EQ(prob_array[51].second,rangebinarysearch(102.00000,prob_array,99));
}

TEST(RangeBinarySearch, first){
  pair<double,my_curve*>* prob_array = new pair<double,my_curve*> [100];
  for(unsigned int i = 0 ;i<100;i++){
    prob_array[i].second = new my_curve(10);
    prob_array[i].first = 2*i;
  }
  //std::cout << prob_array[1].second << '\n';//if error happens
  my_curve* probone = prob_array[1].second;
  ASSERT_EQ(probone,rangebinarysearch(0.1,prob_array,99));
}

TEST(RangeBinarySearch, last){
  pair<double,my_curve*>* prob_array = new pair<double,my_curve*> [100];
  for(unsigned int i = 0 ;i<100;i++){
    prob_array[i].second = new my_curve(10);
    prob_array[i].first = 2*i;
  }
  ASSERT_EQ(prob_array[99].second,rangebinarysearch(198,prob_array,99));
}
//---------------MeanVector----------------//
TEST(MeanVector,simple){
  my_vector A(10);//dim 10
  for (size_t i = 0; i < 10; i++)
  A.coordinates[i] = i;
  my_vector B(10);//dim 10
  for (size_t i = 0; i < 10; i++) {
    B.coordinates[i] = 1;
  }
  vector<my_vector*> L;
  L.push_back(&A);
  L.push_back(&B);
  my_vector expected(10);
  for (size_t i = 0; i < 10; i++) {
    expected.coordinates[i] = (double)i/2.0 + 0.5;
  }
  expected.print_vec();
  get_mean(10,L).print_vec(0);
  std::cout << '\n';
  ASSERT_TRUE(expected == get_mean(10,L));
}

TEST(MeanVector,multiple){
  vector<my_vector*> L ;
  my_vector *t1,*t2;//dim 10
  unsigned int count = 3;
  unsigned int zerocount = 7;

  //----------zero vector---------//
  t1 = new my_vector(10);
  for (size_t i = 0; i < 10; i++) {
    t1->coordinates[i] = 0;
  }
  for (unsigned int j = 0; j < zerocount; j++) {
    L.push_back(t1);
  }
  //--------non zero vector--------//
  t2 = new my_vector(10);
  for (size_t i = 0; i < 10; i++) {
    t2->coordinates[i] = i;
  }
  for (unsigned int j = 0; j < count; j++) {
    L.push_back(t2);
  }
  //-----------expected--------------//
  my_vector expected(10);
  for (size_t i = 0; i < 10; i++) {
    expected.coordinates[i] = (double)i*((double) count/(count+zerocount));
  }
  //-------------cout-----------//
  expected.print_vec();
  get_mean(10,L).print_vec(0);
  std::cout << '\n';
  ASSERT_TRUE(expected == get_mean(10,L));
}
//-------------K++-------------------//
TEST(kmeansplusplus,vectors){
  //from smallvectors test
  unsigned int k = 2;
  list<my_vector> *data_tmp = read_vector_file("./Input/testfar");
  vector<my_vector>* data=new vector<my_vector>;//to metatrepo se vector gt i palia sinartisi diabasmatos epestrefe list
  for(auto i: *data_tmp)
  data->push_back(i);
  data_tmp->clear();
  delete data_tmp;
  vector<my_vector>* result=initialization2(data, k, manhattan_distance);
  std::cout << "points were :";
  (*result)[0].print_vec();
  (*result)[1].print_vec();
  std::cout <<"d="<<manhattan_distance((*result)[0],(*result)[1])<< '\n';

  EXPECT_FALSE(manhattan_distance((*result)[0],(*result)[1]) <= 4.0);
}
//----------------------------------//
int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
