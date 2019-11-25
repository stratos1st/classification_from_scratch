#include "my_curve.hpp"
#include <cmath>
#include <iostream>

using namespace std;

#define DEBUG 0

double my_curve::tol = 0.01;

my_curve::my_curve(unsigned int points,unsigned int dimentions){
  #if DEBUG
  cout<<"Constructing "<<dimentions<<"d curve with "<<points<<" points"<<'\n';
  #endif
  // initializing local variables
  numofvectors=points;
  vectordimentions=dimentions;
  id=0;
  //making an empty array of (pointers to ) vectors
  vectors = new my_vector* [numofvectors];
  for(unsigned int i=0;i<=numofvectors-1;i++)
    vectors[i]=new my_vector(vectordimentions);
}

my_curve::~my_curve(){
  #if DEBUG
  cout<<"Destructing "<<vectordimentions<<"d curve with "<<numofvectors<<" points"<<'\n';
  #endif
  //deleting first each vector and then the whole array
  for(unsigned int i = 0; i<=numofvectors-1;i++)
    delete vectors[i];
  delete[] vectors;
}

my_vector& my_curve::get_vector(unsigned int index){
  return *vectors[index];
}

my_curve& my_curve::operator=(const my_curve &other){
  if(numofvectors!=other.numofvectors){
    cerr<<"\n\n!! ERROR my_curve coppy construstor !!\n\n";
    exit(1);
  }
  #if DEBUG
    cout<<"assinging values with = opperator id curve"<<other.id<<'\n';
  #endif

  for(unsigned int i=0;i<numofvectors;i++)
    *vectors[i]=*other.vectors[i];  //this calls vector '=' operator
  id=other.id;

  vectordimentions=other.vectordimentions;
  return *this;
}

bool my_curve::operator==(const my_curve &other){
  if(numofvectors!=other.numofvectors){
    // cerr<<"\n\n!! ERROR my_curve equal operator !!\n\n";
    // exit(1);
    //OR are not equal
    return false;
  }
  #if DEBUG
    cout<<"compairing values with == opperator id curve"<<other.id<<'\n';
  #endif

  for(unsigned int i=0;i<numofvectors;i++)
    for (unsigned int j = 0; j < vectordimentions; j++)
      if(abs(get_vector(i).coordinates[j]-other.vectors[i]->coordinates[j])>tol)
        return false;
  return true;
}

my_curve::my_curve(const my_curve &p2){
  #if DEBUG
  cout<<"Copy Constructing "<<p2.numofvectors<<" points curve in("<<p2.vectordimentions<<"d)"<<'\n';
  #endif
  // this creates an empty space for the newcurve that will be filled with the '='operator
  vectordimentions=p2.vectordimentions;
  numofvectors=p2.numofvectors;
  vectors = new my_vector* [numofvectors];
  for(unsigned int i=0;i<=numofvectors-1;i++)
    vectors[i]=new my_vector(vectordimentions);
  *this=p2;   //this calls curve '='operator
}

unsigned int my_curve::get_dimentions() const{
  return vectordimentions;//sizeof(coordinates)/sizeof(coordinates[0]);
}

void my_curve::print_vec(unsigned int until){
  cout<<numofvectors<<" curve "<<id<<endl;
  if(until>numofvectors){
    cerr<<"!! ERROR print_vec";
    return;
  }
  if(until==0)
    until=numofvectors;
  for(unsigned int i=0;i<until;i++){
    cout<<"\t";
    vectors[i]->print_vec();
    cout<<", ";
  }
  cout<<endl;
}
