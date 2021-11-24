#include "mesh.h"
#include "poisson.h"
#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
    Mesh* mesh = new Mesh("infile/mesh.in");
    Poisson* poisson = new Poisson(mesh);

    cout << mesh->num_cells() << endl;
    for(int i = 0; i < 10000; ++i)
        poisson->showinstantphi();
    
    return 0;
}