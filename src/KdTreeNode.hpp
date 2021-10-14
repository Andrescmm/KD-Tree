#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
//#include <pair>
#include "Point.hpp"


template <size_t N, typename ElemType>
class KDTreeNode {
  public:
    typedef std::pair<Point<N>, ElemType> value_type;
    KDTreeNode *leftNode;
    KDTreeNode *rightNode;
    ElemType value;
    int level;
    Point<N> p;
    KDTreeNode();
    KDTreeNode(value_type &value,int _level);
    KDTreeNode(const Point<N>& _pt, int _level, const ElemType& _value=ElemType()):
        p(_pt), leftNode(NULL), rightNode(NULL), level(_level), value(_value) {}
    
   
};

template <size_t N, typename ElemType>
KDTreeNode<N, ElemType>::KDTreeNode() {
  leftNode=nullptr;
  rightNode=nullptr;
  this-> level = level;
}
