FMM2DTree::getTransferParentToChild(int N_child, array &L2L)
{
    // Dividing the domain [-1, 1] to [-1, 0] and [0, 1]:
    // Child Nodes which are < 0:
    array child_nodes_1 = 0.5 * (this->standard_nodes - 1);
    // Child Nodes which are > 0:
    array child_nodes_2 = 0.5 * (this->standard_nodes + 1);

    // Getting nodes_x and nodes_y depending upon the quadrant:
    array child_nodes_x, child_nodes_y;

    if(N_child == 1)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_1, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_1.T(), this->N_nodes));
    }

    else if(N_child == 2)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_2, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_1.T(), this->N_nodes));
    }

    else if(N_child == 3)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_2, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_2.T(), this->N_nodes));
    }

    else if(N_child == 4)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_1, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_2.T(), this->N_nodes));
    }

    else
    {
        cout << "INVALID!!!" << endl;
        exit(1);
    }

    getL2L2D(child_nodes_x, child_nodes_y, this->standard_nodes, L2L);
}