#include <iostream>
#include <list>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <float.h>

#include "util.hpp"
#include "my_vector.hpp"
#include "my_curve.hpp"
#include "lsh.hpp"

using namespace std;

#define DEBUG 0

void init1_ass1_update2(list <my_vector>* data,unsigned int k);
template<class T>
list<T>* initialization1(list <T>* data, unsigned int k);
template list<my_vector>* initialization1<>(list <my_vector>* data, unsigned int k);
template list<my_curve>* initialization1<>(list <my_curve>* data, unsigned int k);
list<my_vector*>* assigment1(list<my_vector>* data, list<my_vector>* centers, unsigned int k);
list<my_vector>* update2(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k);

list<my_vector>* update1(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k, lsh_vector* lsh_model);
list<my_curve>* update1_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k, lsh_curve* lsh_model);//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus

int main(){
  unsigned int k=4;

  cout<<"arxisame!!\n";

  list <my_vector>* data=read_vector_file("./test.txt");

  lsh_vector *lsh_model=new lsh_vector(data->front().get_dimentions(),5,6,4,4);
  lsh_model->train(data);

  #if DEBUG
  cout<<"data\n";
  for(auto i : *data)
    i.print_vec();
  #endif

  list <my_vector>* centers=initialization1(data,k);//-----------initialization1

  list<my_vector*> *old_clusters;
  list<my_vector*> *clusters;

  clusters=assigment1(data,centers,k);//------------------------------assigment1

  unsigned int iteration=0;
  while(1){
    cout<<"iteration="<<iteration<<endl;

  // --------------------------------------------------------------sin8iki break
    bool changed=false;
    for(unsigned int i=0;i<k;i++)
      if(old_clusters[i]!=clusters[i]){
        changed=true;
        break;
      }
    if(!changed){
      cout<<"success!!\n";
      cout<<"end centers\n";
      for(auto i : *centers)
        i.print_vec();
      cout<<"clusters\n";
      unsigned int ii=0;
      for(unsigned int i=0;i<k;i++){
        cout<<"\tcluster"<<++ii<<endl;
        for(my_vector* j : clusters[i])
          j->print_vec();
      }
      break;
    }

    centers=update1(data, centers, clusters, k,lsh_model);//-----------------------update2

    //TODO delete old
    old_clusters=clusters;
    clusters=assigment1(data,centers,k);//----------------------------assigment1

    iteration++;
  }

  cout<<"\n\nEND!!\n\n";

  return 0;
}

my_vector get_vector_mean(unsigned int dimentions, list<my_vector*> &cluster){
  my_vector average(dimentions);

  for(unsigned int i=0;i<dimentions;i++)
   average.coordinates[i]=0;

  for(my_vector* j : cluster)
    for(unsigned int i=0;i<dimentions;i++)
      average.coordinates[i]+=j->coordinates[i];
  //FIXEME mpori na ine 0 to clusters[i].size()
  for(unsigned int i=0;i<dimentions;i++)
    average.coordinates[i]/=cluster.size();

  return average;
}

template<class T>
list<T>* initialization1(list <T>* data, unsigned int k){
  // --------------------------------------------------------------------initialization 1 ( nomizo afto lei)
    srand (time(NULL));
    //TODO mpori na gini beltistopiisi. prota epilogi ton kentron kai meta ena perasma tis listas
    list <T>* centers=new list <T>;
    for(unsigned int i=0;i<k;i++)
      centers->push_back(*next(data->begin(), rand()%data->size()));

    #if DEBUG
    cout<<"centers\n";
    for(auto i : *centers)
      i.print_vec();
    #endif
  return centers;
}

list<my_vector*>* assigment1(list<my_vector>* data, list<my_vector>* centers, unsigned int k){
// --------------------------------------------------------------------asigment 1 (nomizo enoi brute force)
  list<my_vector*> *clusters=new list<my_vector*>[k];
  //TODO chech vector typed for more apropriate
  unsigned int cluster_no=0, best_cluster=0;
  double minn=DBL_MAX,tmp=0.0;
  for(auto &i : *data){
    minn=DBL_MAX;
    cluster_no=0;
    best_cluster=0;
    for(auto j : *centers){
      cluster_no++;
      tmp=manhattan_distance(i, j);
      if(tmp<minn){
        best_cluster=cluster_no;
        minn=tmp;
      }
    }
    if(best_cluster==0){
      cerr<<"\n\n!! asigment 1 ERROR cluster not found !!\n\n";
      exit(1);
    }
    clusters[best_cluster-1].push_back(&i);
  }
  //FIXEME mpor;i ena cluster na ine keno!!
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()==0){
      cout<<"cluster of size 0\n";
      exit(1);
    }

  #if DEBUG
  cout<<"clusters\n";
  unsigned int ii=0;
  for(unsigned int i=0;i<k;i++){
    cout<<"\tcluster"<<++ii<<endl;
    for(my_vector* j : clusters[i])
      j->print_vec();
  }
  #endif

  return clusters;
}

list<my_vector>* update2(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k){
// --------------------------------------------------------------------update 2 (nomizo enoi brute force)
  centers->clear();

  for(unsigned int i=0;i<k;i++)
    centers->push_back(get_vector_mean(data->front().get_dimentions(),clusters[i]));

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif

  return centers;
}

void init1_ass1_update2(list <my_vector>* data, unsigned int k){
  list <my_vector>* centers=initialization1(data,k);//-----------initialization1

  list<my_vector*> *old_clusters;
  list<my_vector*> *clusters;

  clusters=assigment1(data,centers,k);//------------------------------assigment1

  unsigned int iteration=0;
  while(1){
    cout<<"iteration="<<iteration<<endl;

  // --------------------------------------------------------------sin8iki break
    bool changed=false;
    for(unsigned int i=0;i<k;i++)
      if(old_clusters[i]!=clusters[i]){
        changed=true;
        break;
      }
    if(!changed){
      cout<<"success!!\n";
      cout<<"end centers\n";
      for(auto i : *centers)
        i.print_vec();
      cout<<"clusters\n";
      unsigned int ii=0;
      for(unsigned int i=0;i<k;i++){
        cout<<"\tcluster"<<++ii<<endl;
        for(my_vector* j : clusters[i])
          j->print_vec();
      }
      break;
    }

    centers=update2(data, centers, clusters, k);//-----------------------update2

    //TODO delete old
    old_clusters=clusters;
    clusters=assigment1(data,centers,k);//----------------------------assigment1

    iteration++;
  }

  cout<<"\n\nEND!!\n\n";
}

list<my_vector>* update1(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k, lsh_vector* lsh_model){//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus
// --------------------------------------------------------------------update 1
  centers->clear();

  for(unsigned int i=0;i<k;i++){
    my_vector average(get_vector_mean(data->front().get_dimentions(),clusters[i]));
    centers->push_back(*lsh_model->find_NN(average, manhattan_distance).first);
  }

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif

  return centers;
}

my_curve get_curve_mean(unsigned int dimentions, list<my_curve*> &cluster){
  //bres to l
  double average=0.0;
  for(my_curve* i: cluster)
    average+=i->numofvectors;
  average/=cluster.size();//FIXEME mpori na ine 0
  unsigned int l=(unsigned int)floor(average);

  //arxikopiise to C
  my_curve *c;
  for(my_curve* it: cluster)//FIXEME perno to proto me length==l eno prepi na perno tixeo
    if(it->numofvectors>=average){
      c=new my_curve(l,dimentions);
      unsigned int jj=rand()%(it->numofvectors-l+1);//ksekina apo mia 8esi
      for(unsigned int j=0;j<l;j++)//ke pigene mexri jj+l
        for(unsigned int i=0;i<dimentions;i++)//gia ka8e diastasi
          c->vectors[j][i]=it->vectors[jj+j][i];//antegrapse sto C
      break;
    }


  my_curve *c2;
  list<my_vector*>* a=new list<my_vector*>[l];
  list<pair<unsigned int,unsigned int>>* ipairs;
  while(1){
    c2=c;
    for(unsigned int i=0;i<l;i++)
      a[i].clear();
    for(my_curve* it: cluster){
      //TODO add best traversal ipairs=best_traversal(c,it);
      for(pair<unsigned int,unsigned int> ij: *ipairs)
        a[ij.first].push_back(it->vectors[ij.second]);
    }
    for(unsigned int i=0;i<l;i++)
      *c->vectors[i]=get_vector_mean(dimentions,a[i]);
    //TODO add == operator if(*c2==*c)
    //   break;

  }

  return *c;
}

list<my_curve>* update1_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k, lsh_curve* lsh_model){//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus
// --------------------------------------------------------------------update 1
  centers->clear();

  for(unsigned int i=0;i<k;i++){
    my_curve average(get_curve_mean(data->front().numofvectors,clusters[i]));
    centers->push_back(*lsh_model->find_NN(average, Dtw, manhattan_distance).first);
  }

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif

  return centers;
}
