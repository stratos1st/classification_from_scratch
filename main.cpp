#include <iostream>
#include <list>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <float.h>
#include <algorithm>
#include <random>

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
list<my_curve*>* assigment1_curve(list<my_curve>* data, list<my_curve>* centers, unsigned int k);
list<my_vector>* update2(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k);
my_curve get_curve_mean(unsigned int dimentions, list<my_curve*> &cluster);
list<my_curve>* update2_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k);

list<my_vector>* update1(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k, lsh_vector* lsh_model);
list<my_curve>* update1_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k, lsh_curve* lsh_model);//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus

double list_diff(list<my_vector> *list1, list<my_vector> *list2);

list<my_vector*>* assigment2(list<my_vector>* data, list<my_vector>* centers, unsigned int k, lsh_vector* lsh_model){
// --------------------------------------------------------------------asigment 2
  #if DEBUG
  cout<<"assigment2\n";
  #endif
  list<my_vector*> *clusters=new list<my_vector*>[k];
  hash<my_vector*> hasher;

  unordered_map<unsigned int, my_vector*> **close_points=new unordered_map<unsigned int, my_vector*>*[k];

  for(unsigned int i=0;i<k;i++){
    close_points[i]=lsh_model->find_bucket(*next(centers->begin(), i), manhattan_distance);
    #if DEBUG
    cout<<"lsh buckets\n-----------"<<i<<endl;
    for (auto v : *close_points[i])
      v.second->print_vec();
    #endif
  }


  /*create uper triangular table with no diagonal
    the table is filled like the example
    ex k=5
      0 1 2 3 4
     -----------
    0|
    1|.
    2|. .
    3|. . .
    4|. . . .

  *///make_pair(kk(it->second),it->second)
  unordered_map<unsigned int, my_vector*> **intersection=new unordered_map<unsigned int, my_vector*>*[k-1];
  for(unsigned int i=0;i<k-1;i++)
    intersection[i]= new unordered_map<unsigned int, my_vector*>[i+1];
  //find intersection of close_points
  for(unsigned int ii=0;ii<k-1;ii++)
    for(unsigned int j=0;j<=ii;j++){//for every . in intersection
      unordered_map<unsigned int, my_vector*>* big=close_points[ii+1];
      unordered_map<unsigned int, my_vector*>* small=close_points[j];
      if(big->size()<small->size())
        swap(big,small);//iterate throw the smallest set
      for(unordered_map<unsigned int, my_vector*>::iterator i = small->begin(); i != small->end(); i++)//for every element in the set
        if(big->find(hasher(i->second)) != big->end())//if the element exists in the other set
          intersection[ii][j].insert(*i);//insert it in intersection
    }

  #if DEBUG
  cout<<"intersection\n";
  for(unsigned int i=0;i<k-1;i++)
    for(unsigned int j=0;j<=i;j++){
      cout<<"\t"<<i+1<<" "<<j<<endl;
      for (auto v : intersection[i][j])
            v.second->print_vec();
    }
  #endif

  unordered_map<unsigned int, my_vector*> left_to_classify;
  for(my_vector& i : *data)//FIXME perini ora. mpori ta data na prepi na ine se set
    left_to_classify.insert(make_pair(hasher(&i),&i));

  double min_dist,tmp1,tmp2;
  unsigned int closest_center;

  cout<<"aaaaaaa\n";
  //classify all intersection elements
  for(unsigned int i=0;i<k-1;i++){
    for(unsigned int j=0;j<=i;j++){
      for(unordered_map<unsigned int, my_vector*>::iterator v=intersection[i][j].begin();v!=intersection[i][j].end();){//for every element in intersection
        if(left_to_classify.find(hasher(v->second)) != left_to_classify.end()){//if the element has yet to be classified
          min_dist=DBL_MAX;
          for(unsigned int ik=0;ik<k-1;ik++){
            for(unsigned int jk=0;jk<=i;jk++){//search all the intersections
              if( intersection[ik][jk].find(hasher(v->second)) !=  intersection[ik][jk].end()){//if it exists in any of them find the min dist and the closest_center
                tmp1=manhattan_distance(*next(centers->begin(), ik), *v->second);
                tmp2=manhattan_distance(*next(centers->begin(), jk), *v->second);
                if(tmp1<min_dist){
                  min_dist=tmp1;
                  closest_center=ik+1;
                }
                if(tmp2<min_dist){
                  min_dist=tmp2;
                  closest_center=jk;
                }
              }
            }
          }
          (v->second)->print_vec();
          clusters[closest_center].push_back(v->second);//push the element in the closest cluster
          left_to_classify.erase(v++);//erase it from left_to_classify
        }
        else
          ++v;
      }
    }
  }

  //classify all close_points elements (exept the ones classified in intersection)
  for(unsigned int i=0;i<k;i++)
    for (auto v : *close_points[i])//for all the close points
      if(left_to_classify.find(hasher(v.second)) != left_to_classify.end()){//if the element has yet to be classified
        clusters[k].push_back(v.second);//push the element in the closest cluster
        left_to_classify.erase(hasher(v.second));//erase it from left_to_classify
      }

  //classify all the rest
  for (unordered_map<unsigned int, my_vector*>::iterator v=left_to_classify.begin();v!=left_to_classify.end();){//for every element in left_to_classify
    min_dist=DBL_MAX;
    for(unsigned int i=0;i<k;i++){//check all the centers
      tmp1=manhattan_distance(*v->second, *next(centers->begin(),i));
      if(min_dist>tmp1){//find the cosest
        min_dist=tmp1;
        closest_center=i;
      }
    }
    clusters[closest_center].push_back(v->second);//push the element in the closest cluster
    left_to_classify.erase(v++);//erase it from left_to_classify
  }

  if(left_to_classify.size()!=0){
    cerr<<"\n\nleft_to_classify.size()!=0\n\n";
    for (auto v : left_to_classify)
      v.second->print_vec();
    exit(1);
  }

  //mpori ena cluster na ine keno!!
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()==0){
      cout<<"cluster of size 0\n";
      // exit(1);
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
  exit(0);
  return clusters;
}

int main(){
  srand (time(NULL));
  unsigned int k=4, max_iterations=10;
  double center_tol=0.01;

  cout<<"arxisame!!\n";

  // list <my_curve>* data=read_curve_file("./test_curves.txt");
  //
  // lsh_curve *lsh_model=new lsh_curve(data->front().vectordimentions,5,5,1000,4,4);
  // lsh_model->train(data);
  //
  // #if DEBUG
  // cout<<"data\n";
  // for(auto i : *data)
  //   i.print_vec();
  // #endif
  //
  // list <my_curve>* centers=initialization1(data,k);//-----------initialization1
  //
  // list<my_curve*> *old_clusters;
  // list<my_curve*> *clusters;


  list <my_vector>* data=read_vector_file("./test.txt");

  lsh_vector *lsh_model=new lsh_vector(data->front().get_dimentions(),5,5,4,4);
  lsh_model->train(data);

  #if DEBUG
  cout<<"data\n";
  for(auto i : *data)
    i.print_vec();
  #endif

  list<my_vector> *centers=initialization1(data,k);//-----------initialization1
  list<my_vector> *old_centers;

  list<my_vector*> *old_clusters;
  list<my_vector*> *clusters;

  //me tis epiloges gia initialization pou dini ipoti8ete oti ta cluster 8a exoun toulaxiston 1 simio
  //kanonika 8a prepi na to ele3oume prin kalesoume assigment
  clusters=assigment2(data,centers,k,lsh_model);//------------------------------assigment1

  unsigned int iteration=0;
  while(max_iterations--!=0){
    cout<<"\niteration="<<iteration<<endl;

  // --------------------------------------------------------------check if cluster==0
  // list<my_curve>::iterator it = next(centers->begin(), i);
  // //------------------------------------------replace center
  // while(next(centers->begin(), i))




    old_centers=centers;
    centers=update2(data, centers, clusters, k);//-----------------------update2

    //sin8iki break:: nea kentra konta sta palia
    // bool changed=false;
    // if(list_diff(old_centers,centers)>center_tol)//prepi na oriso sinartisi gia to - (curves kai vectors)
    //   changed=true;
    // if(!changed){
    //   cout<<"success!!\n";
    //   cout<<"end centers "<<list_diff(old_centers,centers)<<"\n";
    //   for(auto i : *centers)
    //     i.print_vec();
    //   cout<<"clusters\n";
    //   unsigned int ii=0;
    //   for(unsigned int i=0;i<k;i++){
    //     cout<<"\tcluster"<<++ii<<endl;
    //     for(auto* j : clusters[i])
    //       j->print_vec();
    //   }
    //   break;
    // }


    //TODO delete old
    old_clusters=clusters;
    clusters=assigment2(data,centers,k,lsh_model);//----------------------------assigment1

    //sin8iki break:: nea cluster idia me palia
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
        for(auto* j : clusters[i])
          j->print_vec();
      }
      break;
    }


    iteration++;
  }

  cout<<"\n\nEND!!\n\n";

  return 0;
}

my_vector get_vector_mean(unsigned int dimentions, list<my_vector*> &cluster){//DEN prepi to cluster na ine 0!!
  #if DEBUG
  cout<<"get_vector_mean\n";
  #endif
  my_vector average(dimentions);

  for(unsigned int i=0;i<dimentions;i++)
   average.coordinates[i]=0;

  for(my_vector* j : cluster)
    for(unsigned int i=0;i<dimentions;i++)
      average.coordinates[i]+=j->coordinates[i];
  //DEN mpori na ine 0 to clusters[i].size()
  for(unsigned int i=0;i<dimentions;i++)
    average.coordinates[i]/=cluster.size();

  return average;
}

template<class T>
list<T>* initialization1(list <T>* data, unsigned int k){//8eli srand (time(NULL)); apo tin main na kalesti
  // --------------------------------------------------------------------initialization 1 ( nomizo afto lei)
  #if DEBUG
  cout<<"initialization1\n";
  #endif

  //TODO mpori na iparxi kalieros tropos
  std::vector<int> random_numbers;
  for(unsigned int i=0; i<data->size(); i++)
      random_numbers.push_back(i);
  unsigned seed = time(NULL);
  shuffle(random_numbers.begin(), random_numbers.end(), std::default_random_engine(seed));

  //TODO mpori na gini beltistopiisi. prota epilogi ton kentron kai meta ena perasma tis listas
  list <T>* centers=new list <T>;
  for(unsigned int i=0;i<k;i++)
    centers->push_back(*next(data->begin(), random_numbers[i]%data->size()));

  #if DEBUG
  cout<<"centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif
  return centers;
}

list<my_vector*>* assigment1(list<my_vector>* data, list<my_vector>* centers, unsigned int k){
// --------------------------------------------------------------------asigment 1 (nomizo enoi brute force)
  #if DEBUG
  cout<<"assigment1\n";
  #endif
  list<my_vector*> *clusters=new list<my_vector*>[k];
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
  //mpori ena cluster na ine keno!!
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()==0){
      cout<<"cluster of size 0\n";
      // exit(1);
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

list<my_curve*>* assigment1_curve(list<my_curve>* data, list<my_curve>* centers, unsigned int k){//poli mikri alagi apo to assigment1
// --------------------------------------------------------------------asigment 1 (nomizo enoi brute force)
  #if DEBUG
  cout<<"assigment1_curve\n";
  #endif
  list<my_curve*> *clusters=new list<my_curve*>[k];
  unsigned int cluster_no=0, best_cluster=0;
  double minn=DBL_MAX,tmp=0.0;
  for(auto &i : *data){
    minn=DBL_MAX;
    cluster_no=0;
    best_cluster=0;
    for(auto j : *centers){
      cluster_no++;
      tmp=Dtw(i, j);
      if(tmp<minn){
        best_cluster=cluster_no;
        minn=tmp;
      }
    }
    // if(best_cluster==0){
    //   cerr<<"\n\n!! asigment 1 ERROR cluster not found !!\n\n";
    //   exit(1);
    // }
    clusters[best_cluster-1].push_back(&i);
  }
  //mpor;i ena cluster na ine keno!!
  // for(unsigned int i=0;i<k;i++)
  //   if(clusters[i].size()==0){
  //     cout<<"cluster of size 0\n";
  //     exit(1);
  //   }

  #if DEBUG
  cout<<"clusters\n";
  unsigned int ii=0;
  for(unsigned int i=0;i<k;i++){
    cout<<"\tcluster"<<++ii<<endl;
    for(my_curve* j : clusters[i])
      j->print_vec();
  }
  #endif

  return clusters;
}

list<my_vector>* update2(list<my_vector>* data, list<my_vector>* centers, list<my_vector*> *clusters, unsigned int k){
// --------------------------------------------------------------------update 2 (nomizo enoi brute force)
  list<my_vector>* new_centers=new list<my_vector>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[k].size()!=0)
      new_centers->push_back(get_vector_mean(data->front().get_dimentions(),clusters[i]));
    else
      new_centers->push_back(*next(centers->begin(), k));
      //TODO make it funcion

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  //TODO delete centers

  return new_centers;
}

list<my_curve>* update2_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k){
// --------------------------------------------------------------------update 2 (nomizo enoi brute force)
  #if DEBUG
  cout<<"update2_curves\n";
  #endif

  list<my_curve>* new_centers=new list<my_curve>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[k].size()!=0)
      new_centers->push_back(get_curve_mean(data->front().vectordimentions,clusters[i]));
    else
      new_centers->push_back(*next(centers->begin(), k));
      //TODO make it funcion

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
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
  #if DEBUG
  cout<<"update1\n";
  #endif
  list<my_vector>* new_centers=new list<my_vector>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[k].size()!=0){
      my_vector average(get_vector_mean(data->front().get_dimentions(),clusters[i]));
      pair<my_vector*,double> tmp=lsh_model->find_NN(average, manhattan_distance);
      if(tmp.second==DBL_MAX){//TODO what to do if lsh cant find options: brute force or throw errors
        cerr<<"wrong lsh arguments!! no NN found!!\n\n";
        exit(1);
      }
      new_centers->push_back(*tmp.first);
    }
    else
      new_centers->push_back(*next(centers->begin(), k));
      //TODO make it funcion

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
}

my_curve get_curve_mean(unsigned int dimentions, list<my_curve*> &cluster){//DEN prepi to cluster na ine 0!!
  #if DEBUG
  cout<<"get_curve_mean\n";
  #endif
  //bres to l
  double average=0.0;
  for(my_curve* i: cluster)
    average+=i->numofvectors;
  average/=cluster.size();//DEN mpori na ine 0!!
  unsigned int l=(unsigned int)floor(average);
  //arxikopiise to C
  my_curve *c;
  for(my_curve* it: cluster)//FIXME perno to proto me length==l eno prepi na perno tixeo
    if(it->numofvectors>=average){
      c=new my_curve(l,dimentions);
      unsigned int jj=rand()%(it->numofvectors-l+1);//ksekina apo mia 8esi
      for(unsigned int j=0;j<l;j++)//ke pigene mexri jj+l
        for(unsigned int i=0;i<dimentions;i++)//gia ka8e diastasi
          c->vectors[j]->coordinates[i]=it->vectors[jj+j]->coordinates[i];//antegrapse sto C
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
      ipairs=MinMatching(*c,*it);
      for(pair<unsigned int,unsigned int> ij: *ipairs)
        a[ij.first-1].push_back(it->vectors[ij.second-1]);
    }
    for(unsigned int i=0;i<l;i++)
      *c->vectors[i]=get_vector_mean(dimentions,a[i]);
    if(*c2==*c)
      break;

  }

  return *c;
}

list<my_curve>* update1_curves(list<my_curve>* data, list<my_curve>* centers, list<my_curve*> *clusters, unsigned int k, lsh_curve* lsh_model){//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus
// --------------------------------------------------------------------update 1
  #if DEBUG
  cout<<"update1_curves\n";
  #endif
  list<my_curve>* new_centers=new list<my_curve>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[k].size()!=0){
      my_curve average(get_curve_mean(data->front().vectordimentions,clusters[i]));
      pair<my_curve*,double> tmp=lsh_model->find_NN(average, Dtw, manhattan_distance);
      if(tmp.second==DBL_MAX){//TODO what to do if lsh cant find options: brute force or throw errors
        cerr<<"wrong lsh arguments!! no NN found!!\n\n";
        exit(1);
      }
      new_centers->push_back(*tmp.first);
    }
    else
      new_centers->push_back(*next(centers->begin(), k));
      //TODO make it funcion

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
}

double list_diff(list<my_vector> *list1, list<my_vector> *list2){
  double maxx=0;
  list<my_vector>::iterator i=list1->begin();

  for(my_vector j : *list2){
    for(unsigned int k=0;k<i->get_dimentions();k++)
      maxx=max(abs(i->coordinates[k]-j.coordinates[k]),maxx);
    i++;
  }
  return maxx;
}
