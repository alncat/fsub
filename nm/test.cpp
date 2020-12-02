#include "nm.h"
#include <iostream>
#include <vector>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/vec3.h>

using namespace std;
namespace array = scitbx::af;
using scitbx::vec3;
class read_nm_file{
private:
    typedef array::shared<double> normal_modes;
public:
    vec3<normal_modes> r;
    array::shared<vec3<normal_modes> > modes;
    read_nm_file(std::size_t n){
        normal_modes nm1(n, array::init_functor_null<double>()), nm2(n, array::init_functor_null<double>()), nm3(n, array::init_functor_null<double>());
        for(std::size_t i = 0; i < n; i++){
            nm1[i] = i;
            nm2[i] = i;
            nm3[i] = i;
        }
        r[0] = nm1;
        r[1] = nm2;
        r[2] = nm3;
        modes.push_back(r);
   }
};

int main(int argc, char **argv) {
    array::tiny<double, 3> aa(1,2,3);
    array::tiny<double, 3> at(aa);
    array::c_grid<2> grid(11,20);
    array::versa<double, array::c_grid<2> > ab(grid,1);
//    typedef array::tiny<double, 50> normal_modes;
//    normal_modes nm1, nm2, nm3;
//    for(std::size_t i = 0; i < 50; i++){
//        nm1[i] = i;
//        nm2[i] = i;
//        nm3[i] = i;
//    }
//    vec3<normal_modes> r;
//    r[0] = nm1;
//    r[1] = nm2;
//    r[2] = nm3;
//    array::shared<vec3<normal_modes> > modes(1, array::init_functor_null<vec3<normal_modes> >());
//    modes[0] = r;
    read_nm_file file1(50);
    cout << file1.r[0][1] << endl;
    cout << file1.modes[0][0][1] << endl;
    cout << at.size() <<endl;
    cout << at.elems[2] << endl;
    cout << at[2] << endl;
    cout << "grid[1]" << grid[1] <<endl;
    cout << "grid size is "<< grid.size_1d() <<endl;
    cout << "grid index is " << grid(3,4) <<endl;
    cout << "grid index value is " << ab(3,4)<<endl;
    cout << voea::tiny_to_cpp( array::tiny<double, 2>(1.1, 2.2) ) <<endl;
//    cout << modes[0][0][48] <<endl;
    return 0;
}

