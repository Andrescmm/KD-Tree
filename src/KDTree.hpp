//
//  KDTree.hpp
//  KD-Tree Andres
//
//  Created by Andres  on 30/09/21.
//

// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"
#include "KdTreeNode.hpp"

using namespace std;

template <size_t N, typename ElemType>
class KDTree {
    
 public:
  typedef std::pair<Point<N>, ElemType> value_type;
  KDTreeNode<N,ElemType> *root;
  KDTree();
  ~KDTree();
  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);
  size_t dimension() const;
  size_t size() const;
  bool empty() const;
  bool contains(const Point<N> &pt) const;
  void insert(const Point<N> &pt, const ElemType &value=ElemType());
  KDTreeNode<N,ElemType>* busqueda( KDTreeNode<N,ElemType>* temp,const Point<N>& pt) const;
  KDTreeNode<N,ElemType>* niveles(typename vector<pair<Point<N>, ElemType>>::iterator start,typename vector<pair<Point<N>, ElemType>>::iterator end, int currLevel);
  ElemType &operator[](const Point<N> &pt);// valor del punto
  ElemType &at(const Point<N> &pt); // parecido al de arriba
  const ElemType &at(const Point<N> &pt) const; // cuando el KD-tree es de tipo constante llamara a la funcioarriba
   ElemType knn_value(const Point<N> &key, size_t k) const; // el valor de knn mas comun de los k mas cercanos, tipo de dato mas comun, si es empate, cuaquiera
  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;// consulta, regresa el vector con los k mas cercanos
 private:
  ElemType temporal;
  size_t dimension_;
  size_t size_;
    
};


// Constructor
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    dimension_=N;
    size_=0;
    root= nullptr;
}

//Destructor
template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    
}

//Constructor Copia
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    this->dimension_=rhs.dimension();
    this->size_=rhs.size();
    this->root= rhs.root;
}

//


template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    if (this != &rhs) {
        freeResource(root);
        this->root = rhs.root;
        this->size = rhs.size;
    }
    return *this;
}

// Dimenssiones
template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

// Cantidad de elementos insertados 
template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    //cout<<"size= "<<size_<<endl;
  return size_;
}

// Si esta vacio
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    if (size_==0){
        return true;
    }
  return false;
}

template <size_t N, typename ElemType>
KDTreeNode<N,ElemType>* niveles(typename vector<pair<Point<N>, ElemType>>::iterator start,typename vector<pair<Point<N>, ElemType>>::iterator end, int currLevel){
    if (start >= end) return NULL;

    int axis = currLevel % N;
    auto cmp = [axis](const std::pair<Point<N>, ElemType>& p1, const std::pair<Point<N>, ElemType>& p2) {
        return p1.first[axis] < p2.first[axis];
    };
    std::size_t len = end - start;
    auto mid = start + len / 2;
    std::nth_element(start, mid, end, cmp);

    while (mid > start && (mid - 1)->first[axis] == mid->first[axis]) {
        --mid;
    }

    KDTreeNode<N,ElemType>* newNode = new KDTreeNode<N,ElemType>(mid->first, currLevel, mid->second);
    newNode->left = buildTree(start, mid, currLevel + 1);
    newNode->right = buildTree(mid + 1, end, currLevel + 1);
    return newNode;
}

//Busqueda
template <size_t N, typename ElemType>
KDTreeNode<N,ElemType>* KDTree<N, ElemType>::busqueda(KDTreeNode<N,ElemType>* temp, const Point<N>& pt) const{
    if (temp ==NULL || temp->p == pt) {
        return temp;
    }
    
    //int iter = 0;
    
    //while(temp){
    const Point<N>& point1 = temp->p;
    int level = temp ->level;
    if (pt[level % N] < point1[level % N]) {
        return temp ->leftNode == NULL ? temp : busqueda(temp->leftNode ,pt);
    } else {
        return temp->rightNode == NULL ? temp : busqueda(temp->rightNode, pt);
      }
       level++;
    //}
}



template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    auto node = busqueda(root, pt);
    return node != NULL && node->p == pt;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
    
    auto targetNode = busqueda(root, pt);
    if (targetNode == NULL) {
        root = new KDTreeNode<N,ElemType>(pt, 0, value);
        size_ = 1;
    } else {
        if (targetNode->p == pt) {
            targetNode->value = value;
        } else {
            int currLevel = targetNode->level;
            KDTreeNode<N,ElemType>* newNode = new KDTreeNode<N,ElemType>(pt, currLevel + 1, value);
            if (pt[currLevel%N] < targetNode->p[currLevel%N]) {
                targetNode->leftNode = newNode;
            } else {
                targetNode->rightNode = newNode;
            }
            ++size_;
        }
    }
}


template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    auto node = busqueda(root, pt);
    if (node != NULL && node->p == pt) { // Si pt estuviera en el arbol
        return node->value;
    } else { 
        insert(pt);
        if (node == NULL) return root->value;
        else return (node->leftNode != NULL && node->leftNode->p == pt) ? node->leftNode->value: node->rightNode->value;
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    const KDTree<N, ElemType>& constThis = *this;
    return const_cast<ElemType&>(constThis.at(pt));
}


template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    auto node = busqueda(root, pt);
    if (node == NULL || node->p != pt) {
        throw std::out_of_range("Point not found in the KD-Tree");
    } else {
        return node->value;
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,size_t k) const {
    
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_

