#ifndef COLL_TREE_H
#define COLL_TREE_H
class Atoms;
class CollisionTree {
    public:
        CollisionTree(void);
        ~CollisionTree(void);
        double init(Atoms *);
        int compute(Atoms *, double);
        void updatePointers(void);
        void print(void);
    private:
        double _center[3];      //!<\brief Node center coordinates.
        double _size;           //!<\brief Node size.
        CollisionTree *_child;  //!<\brief Node children array.
        CollisionTree *_next;   //!<\brief Pointer for tree walking.
        CollisionTree *_skip;   //!<\brief Pointer for tree walking.
        int _n;                 //!<\brief Node occupation number.
        int _i;                 //!<\brief Index of the first atom.
        int _j;                 //!<\brief Index of the second atom.
};
#endif //COLL_TREE_H
/* coltree.h */
