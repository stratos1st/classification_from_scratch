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
#include "init_k++.hpp"

using namespace std;

#define DEBUG 1

void init1_ass1_update2(list <my_vector>* data,unsigned int k);
template<class T>
list<T>* initialization1(list <T>* data, unsigned int k);
template list<my_vector>* initialization1<>(list <my_vector>* data, unsigned int k);
template list<my_curve>* initialization1<>(list <my_curve>* data, unsigned int k);
list<my_vector*>* assigment1(list<my_vector>* data, list<my_vector>* centers, unsigned int k);
list<my_vector>* update2(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k);

list<my_vector>* update1(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k, lsh_vector* lsh_model);
list<my_curve>* update1_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k, lsh_curve* lsh_model);//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus

//filippos was here
void vectorizeandsilhouette(list<my_vector*>* lclusters,list<my_vector>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
void silhouette(vector<my_vector*>* clusters,vector<my_vector>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);

int main(){
  unsigned int k=20;

  cout<<"arxisame!!\n";

  list <my_vector>* data=read_vector_file("./Input/test.txt");

  lsh_vector *lsh_model=new lsh_vector(data->front().get_dimentions(),5,6,4,4);
  lsh_model->train(data);

  #if DEBUG
  cout<<"data\n";
  for(auto i : *data)
    i.print_vec();
  #endif

  list <my_vector>* centers=initialization2(data,k);//-----------initialization1

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
  vectorizeandsilhouette(clusters,centers,centers->size(),data->size());

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



void silhouette(vector<my_vector*>* clusters,vector<my_vector>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)){
  double *a= new double[n];
  double *b= new double[n];
  unsigned int clustersize = 0;
  unsigned int i = 0;//this is for a[i]
  unsigned int r = 0;//this is for b[r]
  double *clustercoef = new double[k];
  double **distarray;
  for (unsigned int x  = 0; x < n; x++) {
    b[x] = 0.0;
    a[x] = 0.0;
  }
  for (unsigned int x  = 0; x < k; x++) {
    clustercoef[x] = 0.0;
  }

  //NearestC & 2ndNearestC
  for(unsigned int j=0; j<k; j++){//take a center and its cluster(they sould be sync)

    clustersize = clusters[j].size();
    #if DEBUG
      std::cout << "j="<<j<<"k="<<k <<"i,r = "<<i<<","<<r <<"size="<<clustersize<<"with center"<<(*centers)[j].id<< '\n';
    #endif

    //create lower triangular matrix
    distarray = new double* [clustersize-1];//create distance array:for cluster
    if(distarray == NULL)
      std::cout << "bad allocate" << '\n';
    for (unsigned int v = 0; v<clustersize-1; v++) {//z:0->clustersize-2 , lenght clustersize-1
      distarray[v]=new double[v+1];
      if(distarray[v] == NULL)
        std::cout << "bad allocate1" << '\n';
      for (unsigned int l = 0; l < v+1; l++) {//fill all z blocks of mem with doubles
        distarray[v][l] = distance_metric(*clusters[j][v],*clusters[j][l+1]);
      }
    }//done with assigment and creation of dist array
    //---------------------------------------------------compute a[i] for each i in j
    for (unsigned int l = 0; l < clustersize; l++) {
      for (unsigned int m= 0; m < clustersize; m++) {
        if(l>m)
        a[i] += distarray[l-1][m];
        else if(l<m)
        a[i] += distarray[m-1][l];
      }
      if (clustersize>1)
        a[i]=a[i]/(double)(clustersize-1);
      else
        a[i]=0;
      i++;
    }
    //delete dist array
    for (size_t l = 0; l < clustersize-1; l++) {
      delete distarray[l];
    }
    delete distarray;
    //--------------compute b[i] for each i in j-------------------------------
    for (unsigned int v = 0; v < clustersize; v++) {//take a vector of j cluster
      unsigned int secondcenter = k ;// initailized out of bounds on perpuse
      double min=DBL_MAX;
      for (unsigned int c = 0; c < k; c++) {// take a center
        if (c!=j ) {//exclude the center we are in  //maybe add z!=j
          if(distance_metric(*clusters[j][v],(*centers)[c])<min){//find the second best for that v
            min = distance_metric(*clusters[j][v],(*centers)[c]);
            secondcenter = c;
            #if DEBUG
              std::cout << "found the secondcenter of "<<clusters[j][v]->id<<" == "<<(*centers)[c].id<< '\n';
            #endif
          }
        }
      }
      //found second center
      if(secondcenter < k ){//find b[r]
        unsigned int secondcentersize = clusters[secondcenter].size();
        for (unsigned int v2= 0; v2 < secondcentersize; v2++) {//take a vector from the 2nd best cluster
          b[r]+=distance_metric(*clusters[j][v],*clusters[secondcenter][v2]);//distance from v to every other v on the second best cluster
        }
        b[r]=b[r]/(double)secondcentersize;
        r++;
      }
      else
        std::cout << "error did not find 2nd center of vector "<<clusters[j][v]->id << '\n';
    }
  }

  // ---------------------------------------cumpute clustercoef[j]
  clustersize = 0;
  i = 0;   // i goes from 0->n because sum(clustersize[j]) over j:1->k == n
  unsigned int sum = 0;
  for (unsigned int j = 0; j < k; j++) {
    clustersize = clusters[j].size();
    for (; i < sum+clustersize; i++) {
      // this is s(i)
      if(a[i]<b[i]){
        #if DEBUG
          std::cout << b[i]<<"+"<<a[i] << '\n';
          std::cout << 1-(a[i]/b[i]) << '\n';
        #endif
        clustercoef[j]+=1-(a[i]/b[i]);
      }
      else if(a[i]>b[i])
      {
        #if DEBUG
          std::cout << b[i]<<"-"<<a[i] << '\n';
          std::cout << (b[i]/a[i])-1 << '\n';
        #endif
        clustercoef[j]+=(b[i]/a[i])-1;
      }
      else
        clustercoef[j]+=0;
    }
    sum+=clustersize;
    clustercoef[j]=clustercoef[j]/(double)(clustersize);
  }
  #if DEBUG
    for (unsigned int x = 0; x < n; x++) {
      cout<<x+1<<". "<<b[x]<<endl;
    }
    for (unsigned int x = 0; x < n; x++) {
      cout<<x+1<<". "<<a[x]<<endl;
    }
  #endif
  std::cout << "\nCluster Coefficients" << '\n';
  std::cout << "--------------------" << '\n';
  for (unsigned int x = 0; x < k; x++) {
    cout<<x+1<<". "<<clustercoef[x]<<endl;
  }
  //TODO delete a ,b ,clustercoef
}

void vectorizeandsilhouette(list<my_vector*>* lclusters,list<my_vector>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)){
  vector<my_vector> centers(lcenters->begin(), lcenters->end());
  vector<my_vector*>* clusters = new vector<my_vector*> [k];
  for (unsigned int j = 0; j < k; j++) {
    for (my_vector* c: lclusters[j]) {//naive approach
  	   clusters[j].push_back(c);
  	}
  }

  //calls silhouette
  silhouette(clusters,&centers,k,n,distance_metric);
  // TODO delete
}
