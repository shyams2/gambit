#include <iostream>
#include <arrayfire.h>
#include <vector>;

class Node 
{
  
    private:
    int box_number;
    int parent_number;
    int children_numbers[4];
    int neighbor_numbers[8];
    int inner_numbers[16];
    int outer_numbers[24];

    Eigen::VectorXd multipoles;
    Eigen::VectorXd locals;

    std::vector cheb_nodes_x, cheb_nodes_y;

    public:
    // Constructor for the class
    Node();
};

// Assigning all values to -1 upon initialization:
Node::Node()
{
    this->box_number    = -1;
    this->parent_number = -1;

    for (unsigned i = 0; i < 4; i++) 
    {
        this->children_numbers[i] = -1;
    }

    for (unsigned i = 0; i < 8; i++) 
    {
        this->neighbor_numbers[i] = -1;
    }

    for (unsigned i = 0; i < 16; i++) 
    {
        this->inner_numbers[i] = -1;
    }

    for (unsigned i = 0; i < 24; i++) 
    {
        this->outer_numbers[i] = -1;
    }
}
