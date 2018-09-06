#include <iostream>
#include <arrayfire.h>
#include <vector>
#include <memory> // for using smart pointers
#include "MatrixData.hpp"

class FMM2DBox 
{
public:
    double c_x, c_y, r_x, r_y; // centers and radii for this box
    bool is_root, is_leaf, is_empty, charge_computed;    
    array scaled_nodes_x, scaled_nodes_y;
    array inds_in_box;

    int N;          // total number of points in this box
    int N_level;    // level in the tree
    int N_assigned; // node number assigned

    // By using smart pointers, memory doesn't have to be manually deleted:
    std::shared_ptr<FMM2DBox> parent;
    std::vector<std::shared_ptr<FMM2DBox>> children(4);
    std::vector<std::shared_ptr<FMM2DBox>> neighbor(8);
    std::vector<std::shared_ptr<FMM2DBox>> interaction(27);

    int N_interaction_allocated; // number of points which aren't nullptr in interaction
    int N_neighbor_allocated; // number of points which aren't nullptr in neighbor    
    double charge;
    double potential;
    array node_potential;
    array L2L;

    // Constructors for the class
    FMM2DBox(){};
    // Destructor for the class:
    ~FMM2DBox(){};
};
