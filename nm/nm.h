#ifndef MMTBX_NM_H
#define MMTBX_NM_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/versa_matrix.h>
#include <vector>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/targets.h>
#include <scitbx/matrix/outer_product.h>
#include <iotbx/pdb/hierarchy.h>
#include <fstream>

using namespace std;
namespace mmtbx { namespace nm {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;
typedef af::shared<vec3<double> > normal_mode;

class nm_init {
private:
    af::shared<normal_mode> modes;
    af::shared<normal_mode> zero_modes;
    af::shared<int> nm_index;
    std::size_t count_zero_modes;
    af::shared<af::shared<bool> > selections;

public:
    void is_na(char* res){
        if(strcmp(res, " DA") == 0 or strcmp(res, " DC") == 0 or strcmp(res, " DG") == 0 or strcmp(res, " DI") == 0
                or strcmp(res, " DT") == 0 or strcmp(res, " DU") == 0){
           res[1] = '*';
        }
        else if(strcmp(res, "  A") == 0 or strcmp(res, "  C") == 0 or strcmp(res, "  G") == 0 or strcmp(res, "  I") == 0
                or strcmp(res, "  T") == 0 or strcmp(res, "  U") == 0){
           res[1] = '*';
        }
    }
    nm_init(
                const char* filename,
                std::size_t n_modes,
                af::shared<iotbx::pdb::hierarchy::atom> const& atoms,
                bool zero_mode_input_flag,
                bool zero_mode_flag)
    {
//      modes(atoms.size(), af::init_functor_null<vec3<normal_modes> >());
	count_zero_modes = 0;
        vec3<double> null_mode_i(0.0, 0.0, 0.0);
        normal_mode null_mode;
        nm_index.resize(atoms.size(), -1);
//        for(std::size_t i = 0; i < atoms.size(); i++){
//            null_mode.push_back(null_mode_i);
//        }
//        vec3<normal_modes> null_mode;
//        for(std::size_t i = 0; i < 3; i++){
//           null_mode[i] = null_mode_i;
//       }
        std::size_t nmode_start = 0;
        if(zero_mode_flag and (not zero_mode_input_flag)){
            nmode_start = 6;
        }
        else {
            nmode_start = 0;
        }
        ifstream file(filename, ios::in | ios::binary);
        std::size_t total_evc = 0;
        std::size_t align_mark = 0;
        std::size_t resi_count = 0;
        if( file.is_open() ){
            std::size_t i = 0;
//            for(std::size_t i = 0; i <= atoms.size(); i++){
            while(1){
                char serial[9] = "";
                file.read(serial, 8);
                file.seekg(5, ios::cur);
//                char * resname = new char[8];
                char resname[4] = "";
                file.read(resname, 3);
                file.seekg(4, ios::cur);
                char atmname[5] = "";
//                char * atmname = new char[8];
                file.read(atmname, 4);
//                cout << serial << resname << atmname << endl;
                for(std::size_t j = align_mark; j < atoms.size(); j++){
                    char res_pdb[4] = "";
//                    char * res_pdb = atoms[j].parent()->data->resname.elems;
                    int charp = 0;
                    while(atoms[j].parent()->data->resname.elems[charp] != '\0'){
                        res_pdb[charp] = atoms[j].parent()->data->resname.elems[charp];
                        charp++;
                    }
                    is_na(res_pdb);
//                    cout << j <<  atoms[j].data->name.elems << "," << atmname << endl;
//                    cout << atoms[j].parent()->data->resname.elems << "," << resname << endl;
                    if( atoms[j].parent()->data->altloc != "" and atoms[j].parent()->data->altloc != "A" ){
                        std::size_t n_atoms = atoms[j].parent()->atoms_size();
                        align_mark += n_atoms;
                        j += n_atoms - 1;
                    }
//                    else if( atoms[j].parent()->data->resname == resname ){
                    else if(strcmp(res_pdb, resname) == 0){
                        std::size_t n_atoms = atoms[j].parent()->atoms_size();
                        if(n_atoms == 3 and (strcmp(resname, "ACE") != 0)){
                            n_atoms += atoms[j+3].parent()->atoms_size();
                        }
                        else if(atoms[j].parent()->data->altloc == "A"){
                            n_atoms += 3;
                        }
                        for(std::size_t k = 0; k < n_atoms; k++){
                            if( atoms[k+align_mark].data->name == atmname ){
//                                cout << atoms[k+align_mark].data->name.elems << "," << atmname << endl;
 //                               cout << atoms[k+align_mark].parent()->data->resname.elems << "," << resname << endl;
                                nm_index[k+align_mark] = i;
                                resi_count++;
//                                cout << resi_count << "," << n_atoms << endl;
                                if(n_atoms == resi_count){
                                    align_mark += n_atoms;
                                    resi_count = 0;
                                }
                                break;
                            }
                        }
                        break;
                    }
                    else{
                        std::size_t n_atoms = atoms[j].parent()->atoms_size();
                        j += n_atoms - 1;
                        align_mark += n_atoms;
                    }
                }
                if(strncmp(serial, "END", 3) == 0){
                    total_evc = i;
                    break;
                }
                i++;
            }
//            for(std::size_t i = 0; i < nm_index.size(); i++){
//                if(nm_index[i] == -1){
//                    cout << i << "," << nm_index[i] << endl;
//                }
//            }
//            for(std::size_t i = 0; i < atoms.size(); i++){
//                if(atoms[i].parent()->atoms_size() == 3 ){
//                    cout << atoms[i].data->name.elems << atoms[i].parent()->data->resname.elems << nm_index[i] << endl;
//                }
//                else if(atoms[i].parent()->data->altloc == "A"){
//                    cout << atoms[i].data->name.elems << atoms[i].parent()->data->altloc.elems << atoms[i].parent()->atoms_size() << nm_index[i] << endl;
//                }
//                else if(atoms[i].parent()->data->altloc == "B"){
//                    cout << atoms[i].data->name.elems << atoms[i].parent()->data->altloc.elems << atoms[i].parent()->atoms_size() << nm_index[i] << endl;
//                }
//            }

            cout << "number of atoms in evecin: " << total_evc << "," <<"in pdb: " <<  atoms.size() << endl;
            for(std::size_t i = 0; i < nmode_start; i++){
                normal_mode null_mode(atoms.size(), null_mode_i);
                modes.push_back(null_mode);
            }
            for(std::size_t i = nmode_start; i < n_modes; i++){
//                normal_mode tmp_mode(atoms.size(), af::init_functor_null<vec3<double> >());
                normal_mode tmp_mode;
                tmp_mode.resize(atoms.size(), null_mode_i);
                for(std::size_t j = 0; j < total_evc; j++){
                    double vx, vy, vz;
                    file.read( (char *) &vx, 8);
                    file.read( (char *) &vy, 8);
                    file.read( (char *) &vz, 8);
                    for(std::size_t k = 0; k < atoms.size() ; k++){
                        if( nm_index[k] == j ){
                            tmp_mode[k][0] = vx;
                            tmp_mode[k][1] = vy;
                            tmp_mode[k][2] = vz;
                            break;
                        }
                    }
                }
//                cout << modes[i-1][0][0] << modes[i-1][0][1] << modes[i-1][0][2] <<endl;
                double vx, vy, vz;
                file.read((char *) &vx, 8);
                file.read((char *) &vy, 8);
                file.read((char *) &vz, 8);
                if( vx!=111111111. and vy!=-222222222. and vz!=333333333. ){
                    cout << "error in reading eigenvector file" <<endl;
                    cout << "please check consistency!!!!!!" << endl;
                    cout << vx << vy << vz <<endl;
                    exit(EXIT_FAILURE);
                }
                else{
                    cout << "finish reading mode " << i <<endl;
                    modes.push_back(tmp_mode);
                }

//                for( std::size_t w = 0; w <= i; w++ ){
//                    cout<< modes[w][0][0] << modes[w][0][1] << modes[w][0][2] << endl;
//                }
            }//read non-zero mode
//            if( nmode_start != 0){
//                gen_zero_modes(modes);
//            } move to python
        }
        else {
            cout << "cannot open eigenvector file" << filename << endl;
            exit(EXIT_FAILURE);
        }
        for(std::size_t i = 0; i < nm_index.size(); i++){
            if( nm_index[i]== -1 ){
//                cout << i << nm_index[i] << endl;
//               cout << "cannot find " << atoms[i].data->name.elems << atoms[i].parent()->data->resname.elems << atoms[i].parent()->data->altloc.elems << atoms[i].parent()->parent()->data->resseq.elems << atoms[i].parent()->parent()->parent()->data->id << " in eigenvector file, it is removed from your selection!" << endl;
	    }
        }
    }
    void schmidt(af::shared<normal_mode> & pall)
    {
//        cout << pall.size() << endl;
        MMTBX_ASSERT(pall.size() == 6*(count_zero_modes+1));
        double anorm = 0.0;
        std::size_t padd = count_zero_modes*6;
//        cout << padd << endl;
//        cout << pall[0+padd].size() << endl;
        for(std::size_t i = 0; i < pall[0+padd].size() ; i++){
            //anorm += pall[padd][i][0]*pall[padd][i][0];
            anorm += pall[padd][i].length_sq();
        }
        MMTBX_ASSERT(anorm != 0.0);
        anorm = 1./sqrt(anorm);

        for(std::size_t i = 0; i < pall[padd].size() ; i++){
            //pall[padd][i][0] *= anorm;
            pall[padd][i] *= anorm;
        }
        double rec[6][6];
        for(std::size_t i = 1; i < pall.size()/(count_zero_modes+1) ; i++){
            for(std::size_t j = 0; j < i; j++){
                rec[j][i] = 0.0;
                for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                    rec[j][i] += pall[j+padd][k][0]*pall[i+padd][k][0] + pall[j+padd][k][1]*pall[i+padd][k][1] + pall[j+padd][k][2]*pall[i+padd][k][2];
                }
            }
            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                double aaa = 0.0;
                for(std::size_t j = 0; j < i; j++){
                    aaa += (pall[j+padd][k][0] + pall[j+padd][k][1] + pall[j+padd][k][2])*rec[j][i];
                }
                pall[i+padd][k][0] -= aaa;
                pall[i+padd][k][1] -= aaa;
                pall[i+padd][k][2] -= aaa;
            }

            anorm = 0.0;
            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                anorm += pall[i+padd][k].length_sq();
            }
            MMTBX_ASSERT(anorm != 0);
            anorm = 1./sqrt(anorm);

            for(std::size_t k = 0; k < pall[i+padd].size(); k++){
                pall[i+padd][k] *= anorm;
            }
        }
    }
    void gen_zero_modes(af::shared<vec3<double> > const& sites_cart,
                        af::shared<double> const& weights,
                        af::shared<bool> const& selection)
    {
        MMTBX_ASSERT(sites_cart.size() == weights.size() );
        MMTBX_ASSERT(selection.size() == nm_index.size() );
        double xcm = 0.0;
        double ycm = 0.0;
        double zcm = 0.0;
        double tmass = 0.0;
        af::shared<vec3<double> > sites_cart_new(sites_cart.size(), af::init_functor_null<vec3<double> >());
        af::shared<bool> selection_new(nm_index.size(), true);
        for( std::size_t i = 0; i < sites_cart.size() ; i++ ){
            double weight = weights[i];
            //double weight = 1.;
            vec3<double> site = sites_cart[i];
            sites_cart_new[i] = sites_cart[i];
            xcm += site[0] * weight;
            ycm += site[1] * weight;
            zcm += site[2] * weight;
            tmass += weight;
        }
        MMTBX_ASSERT(tmass != 0);
        xcm /= tmass;
        ycm /= tmass;
        zcm /= tmass;
        for( std::size_t i = 0; i < sites_cart_new.size() ; i++){
            sites_cart_new[i][0] -= xcm;
            sites_cart_new[i][1] -= ycm;
            sites_cart_new[i][2] -= zcm;
        }
        for( std::size_t i = 0; i < 6 ; i++ ){
            af::shared<vec3<double> > null_mode( sites_cart.size(), af::init_functor_null<vec3<double> >());
            zero_modes.push_back(null_mode);
        }
        std::size_t padd = count_zero_modes*6;
        for( std::size_t i = 0; i < 6 ; i++ ){
            for( std::size_t j = 0 ; j < sites_cart.size(); j++ ){
                zero_modes[i+padd][j][0] = 0.0;
                zero_modes[i+padd][j][1] = 0.0;
                zero_modes[i+padd][j][2] = 0.0;
            }
        }
        std::size_t i = 0;
        std::size_t i_skip = 0;
        for( std::size_t j = 0; j < nm_index.size() ; j++ ){
            if( selection[j] == true and nm_index[j] != -1 ){
                double sqrt_weight = sqrt(weights[i]);
                //double sqrt_weight = 1.0;
                zero_modes[0+padd][i][0] = sqrt_weight;
                zero_modes[1+padd][i][1] = sqrt_weight;
                zero_modes[2+padd][i][2] = sqrt_weight;
                zero_modes[3+padd][i][1] = sqrt_weight*sites_cart_new[i][2];
                zero_modes[3+padd][i][2] = -sqrt_weight*sites_cart_new[i][1];
                zero_modes[4+padd][i][0] = -sqrt_weight*sites_cart_new[i][2];
                zero_modes[4+padd][i][2] = sqrt_weight*sites_cart_new[i][0];
                zero_modes[5+padd][i][0] = sqrt_weight*sites_cart_new[i][1];
                zero_modes[5+padd][i][1] = -sqrt_weight*sites_cart_new[i][0];
//                zero_modes[0+padd][i][0] = sqrt_weight;
//                zero_modes[1+padd][i][1] = sqrt_weight;
//                zero_modes[2+padd][i][2] = sqrt_weight;
//                zero_modes[3+padd][i][1] = -sqrt_weight*sites_cart_new[i][2];
//                zero_modes[3+padd][i][2] = sqrt_weight*sites_cart_new[i][1];
//                zero_modes[4+padd][i][0] = sqrt_weight*sites_cart_new[i][2];
//                zero_modes[4+padd][i][2] = -sqrt_weight*sites_cart_new[i][0];
//                zero_modes[5+padd][i][0] = -sqrt_weight*sites_cart_new[i][1];
//                zero_modes[5+padd][i][1] = sqrt_weight*sites_cart_new[i][0];
                i++;
                selection_new[j] = true;
            }
            else if( selection[j] == true and nm_index[j] == -1 ){
                i++;
                i_skip++;
                selection_new[j] = false;
            }
            else {
                selection_new[j] = false;
            }
        }
        //cout << "new selection" << i - i_skip << "," << selection.size() << endl;
        schmidt(zero_modes);
        count_zero_modes++;
        selections.push_back(selection_new);
//        cout << count_zero_modes << endl;
    }

    void normalize(std::size_t n_modes){
        MMTBX_ASSERT( n_modes = modes.size() );
        for(std::size_t i = 0; i < n_modes; i++){
            double norm_mode = 0.0;
            for(std::size_t j = 0; j < modes[i].size(); j++){
                norm_mode += modes[i][j].length_sq();
            }
            double sqrt_norm = sqrt(norm_mode);
            for(std::size_t j = 0; j < modes[i].size(); j++){
                modes[i][j] /= sqrt_norm;
            }
        }
    }

    af::shared<vec3<double> > return_modes( std::size_t i ) { return modes[i]; }
    af::shared<vec3<double> > return_zero_modes( std::size_t i ) { return zero_modes[i]; }
    af::shared<bool> return_new_selection( std::size_t i ) { return selections[i]; }
    void print_eigenvector( std::size_t i )
    {
        for( std::size_t j = 0; j < modes[i].size(); j++ ){
            cout<< modes[i][j][0] << modes[i][j][1] << modes[i][j][2] << endl;
        }
    }
};

void normalize_mode(af::shared<vec3<double> > & mode){
    double norm_mode = 0.0;
    for(std::size_t j = 0; j < mode.size(); j++){
        norm_mode += mode[j].length_sq();
    }
    MMTBX_ASSERT(norm_mode != 0);
    double sqrt_norm = sqrt(norm_mode);
    for(std::size_t j = 0; j < mode.size(); j++){
        mode[j] /= sqrt_norm;
    }
}

af::shared<sym_mat3<double> > init_nm_adp(af::shared<vec3<double> > const& modes,
                                            af::shared<double> const& weights,
                                            std::size_t n_modes,
                                            bool zero_mode_flag)
{
//    af::shared<sym_mat3<double> > adp_nma;
    std::size_t unit_len = weights.size();
    std::size_t total_len = unit_len*n_modes;
//    cout << unit_len << endl;
//    cout << modes.size() << endl;
    MMTBX_ASSERT(modes.size() == total_len);
    af::shared<sym_mat3<double> > adp_nma(n_modes*total_len, af::init_functor_null<sym_mat3<double> >());
//    cout << "all right" << endl;
//    af::shared<sym_mat3<double> > adp_nma_i(modes[0].size(), af::init_functor_null<sym_mat3<double> >());
    for(std::size_t n = 0; n < n_modes ; n++){
        for(std::size_t m = n; m < n_modes ; m++){
//            af::shared<sym_mat3<double> > adp_nma_i(modes[n].size(), af::init_functor_null<sym_mat3<double> >());
            if(zero_mode_flag){
                if( n < 6 and m > 7 )
                    continue;
            }
            for( std::size_t i = 0; i < weights.size() ; i++){
                sym_mat3<double> adp;
                adp[0] = modes[n*unit_len + i][0]*modes[m*unit_len + i][0];
                adp[1] = modes[n*unit_len + i][1]*modes[m*unit_len + i][1];
                adp[2] = modes[n*unit_len + i][2]*modes[m*unit_len + i][2];
                adp[3] = (modes[n*unit_len + i][0]*modes[m*unit_len + i][1] + modes[n*unit_len + i][1]*modes[m*unit_len + i][0])/2.0;
                adp[4] = (modes[n*unit_len + i][0]*modes[m*unit_len + i][2] + modes[n*unit_len + i][2]*modes[m*unit_len + i][0])/2.0;
                adp[5] = (modes[n*unit_len + i][1]*modes[m*unit_len + i][2] + modes[n*unit_len + i][2]*modes[m*unit_len + i][1])/2.0;
                adp = adp/weights[i];
                adp_nma[i + m*unit_len + n_modes*n*unit_len] = adp;
//                adp_nma.push_back(adp);
            }
        }
    }
    return adp_nma;
}

af::versa<double, af::c_grid<2> > unpack_x(af::shared<double> const& x, std::size_t n_modes, bool zero_mode_flag)
{
//    MMTBX_ASSERT(x.size() == n_modes*n_modes);
    af::versa<double, af::c_grid<2> > s(af::c_grid<2>(n_modes, n_modes), 0.0);
    if(zero_mode_flag){
        std::size_t nx = 0;
        for(std::size_t i = 0; i < 6; i++){
            for(std::size_t j = 0; j <= i; j++){
                s(i, j) = x[nx];
                nx++;
            }
        }
        for(std::size_t i = 6; i < n_modes; i++){
            for(std::size_t j = 6; j <= i; j++){
                s(i, j) = x[nx];
                nx++;
            }
        }
    }
    else{
        std::size_t nx = 0;
        for(std::size_t i = 0; i < n_modes; i++){
            for(std::size_t j = 0; j <= i; j++){
                s(i, j) = x[nx];
                nx++;
            }
        }
    }

    return s;
}

af::shared<double> scale_x(af::shared<double> const& x,
                           af::shared<sym_mat3<double> > const& uanisos,
                           af::shared<sym_mat3<double> > const& adp_all,
                           std::size_t n_modes,
                           bool zero_mode_flag)
{
    MMTBX_ASSERT(uanisos.size() == adp_all.size());
    double iso_all = 0.0;
    double aniso_all = 0.0;
    af::shared<sym_mat3<double> > adp_nma;
    for(std::size_t i = 0; i < uanisos.size(); i++){
//        if(uanisos[i][2] == 0.0){
//            iso_all += uanisos[i][0]*3.;
//        }
//        else{
          for(std::size_t j = 0; j < 3; j++){
                iso_all += uanisos[i][j];
            }
//        }
        for(std::size_t j = 0; j < 3; j++){
            aniso_all += adp_all[i][j];
        }
    }
    aniso_all = aniso_all/3.0/uanisos.size();
    iso_all = iso_all/3.0/uanisos.size();
// ios_all is slight incorrect, but it's ok
    double scaling = sqrt(iso_all/aniso_all);
    af::shared<double> scaled_x;
    scaled_x = x*scaling;
    cout << "x is scaled by " << scaling << endl;
//    for(std::size_t i = 0; i < x.size(); i++){
//        scaled_x[i] = x[i]*scaling;
//    }
    return scaled_x;
}

class uaniso_from_s {
private:
    af::versa<double, af::c_grid<2> > sigma;
    af::shared<sym_mat3<double> > adp_all;
    af::versa<double, af::c_grid<2> > s;
public:
    void s2sigma(af::versa<double, af::c_grid<2> > const& s, bool zero_mode_flag)
    {
        MMTBX_ASSERT(s.accessor().is_square());
        std::size_t n_modes = s.accessor()[0];
        sigma.resize(af::c_grid<2>(n_modes, n_modes), 0.0);
        if(zero_mode_flag){
            for( std::size_t i = 0; i < 6; i++ ){
                for( std::size_t j = i; j < 6; j++ ){
                    for( std::size_t k = 0; k <= i; k++ ){
                        sigma(i, j) += s(i, k)*s(j, k);
                    }
                }
            }
            for( std::size_t i = 6; i < n_modes ; i++ ){
                for( std::size_t j = i; j < n_modes ; j++ ){
                    for( std::size_t k = 0; k <= i ; k++ ){
                        sigma(i, j) += s(i, k)*s(j, k);
                    }
                }
            }
        }
        else{
            for( std::size_t i = 0; i < n_modes ; i++ ){
                for( std::size_t j = i; j < n_modes ; j++ ){
                    for( std::size_t k = 0; k <= i ; k++ ){
                        sigma(i, j) +=  s(i, k)*s(j, k);
                    }
                }
            }
        }
    }

    uaniso_from_s(af::shared<double> const& x,
                  af::shared<sym_mat3<double> > const& adp_nma,
                  std::size_t n_modes,
                  bool zero_mode_flag)
    {
//        MMTBX_ASSERT(adp_nma.size() == n_modes*n_modes*weights.size());
        s = unpack_x(x, n_modes, zero_mode_flag);
//        sigma = af::matrix_multiply(s, s);
        s2sigma(s, zero_mode_flag);
//        af::shared<sym_mat3<double> > adp_nma;
//        adp_nma = init_nm_adp(modes1d, weights, n_modes, zero_mode_flag);
        sym_mat3<double> nul_sym_mat3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        std::size_t total_atoms = adp_nma.size()/(n_modes*n_modes);
        adp_all.resize(total_atoms, nul_sym_mat3);
        for(std::size_t n = 0; n < total_atoms ; n++){
            for(std::size_t i = 0; i < n_modes ; i++){
                adp_all[n] += sigma(i, i)*adp_nma[n + i*n_modes*total_atoms + i*total_atoms];
            }
        }
        if(zero_mode_flag){
            for(std::size_t n = 0; n < total_atoms; n++){
                for(std::size_t i = 0; i < 5; i++){
                    for(std::size_t j = i + 1; j < 6; j++){
                        adp_all[n] += 2*sigma(i, j)*adp_nma[n + i*n_modes*total_atoms + j*total_atoms];
                    }
                }
                for(std::size_t i = 6; i < n_modes - 1; i++){
                    for(std::size_t j = i + 1; j < n_modes; j++){
                        adp_all[n] += 2*sigma(i, j)*adp_nma[n + i*n_modes*total_atoms + j*total_atoms];
                    }
                }
            }
        }
        else{
            for(std::size_t n = 0; n < total_atoms; n++){
                for(std::size_t i = 0; i < n_modes - 1 ; i++){
                    for(std::size_t j = i + 1; j < n_modes ; j++){
                        adp_all[n] += 2*sigma(i, j)*adp_nma[n + i*n_modes*total_atoms + j*total_atoms];
                    }
                }
            }
        }
    }
    af::shared<sym_mat3<double> > u_cart() { return adp_all; }
};

class d_uaniso_d_nm {
public:
    d_uaniso_d_nm(af::shared<sym_mat3<double> > const& adp_nma,
                  af::versa<double, af::c_grid<2> > const& s,
                  std::size_t atm_num,
                  std::size_t n_modes,
                  std::size_t unit_len,
                  bool zero_mode_flag)
    {
        MMTBX_ASSERT(s.accessor().is_square());
        MMTBX_ASSERT(s.accessor()[0] == n_modes);
        MMTBX_ASSERT(adp_nma.size() == n_modes*unit_len*n_modes);
        if(zero_mode_flag){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j <= i; j++){
                    sym_mat3<double> dpart(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    for(int k = i; k < 6; k++){
                        dpart += s(k, j)*adp_nma[atm_num + i*n_modes*unit_len + k*unit_len];
                    }
                    for(int k = j; k < i; k++){
//                        cout << "all right! k, i, j" << k <<"," << i << "," << j << endl;
                        dpart += s(k, j)*adp_nma[atm_num + k*n_modes*unit_len + i*unit_len];
                    }
                    dpart *= 2.0;
                    dudnm_.push_back(dpart);
                }
            }
            for(int i = 6; i < n_modes; i++){
                for(int j = 6; j <= i; j++){
                    sym_mat3<double> dpart(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    for(int k = i; k < n_modes; k++){
                        dpart += s(k, j)*adp_nma[atm_num + i*n_modes*unit_len + k*unit_len];
                    }
                    for(int k = j; k < i; k++){
                        dpart += s(k, j)*adp_nma[atm_num + k*n_modes*unit_len + i*unit_len];
                    }
                    dpart *= 2.0;
                    dudnm_.push_back(dpart);
                }
            }
            MMTBX_ASSERT(dudnm_.size() == 21 + (n_modes - 5)*(n_modes - 6)/2);
        }
        else{
            for(int i = 0; i < n_modes; i++){
                for(int j = 0; j <= i; j++){
                    sym_mat3<double> dpart(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    for(int k = i; k < n_modes; k++){
                        dpart += s(k, j)*adp_nma[atm_num + i*n_modes*unit_len + k*unit_len];
                    }
                    for(int k = j; k < i; k++){
                        dpart += s(k, j)*adp_nma[atm_num + k*n_modes*unit_len + i*unit_len];
                    }
                    dpart *= 2.0;
                    dudnm_.push_back(dpart);
                }
            }
            MMTBX_ASSERT(dudnm_.size() == (n_modes + 1)*n_modes/2);
        }
    }
    af::shared<sym_mat3<double> > d_u_d_nm() { return dudnm_; }
private:
        af::shared<sym_mat3<double> > dudnm_;
};

class d_target_d_nm {
public:
    d_target_d_nm(af::shared<sym_mat3<double> > const& adp_nma,
                  af::shared<sym_mat3<double> > const& d_target_d_uaniso,
                  af::shared<double> const& x,
                  std::size_t n_modes,
                  bool zero_mode_flag)
    {
        std::size_t unit_len = d_target_d_uaniso.size();
        MMTBX_ASSERT( adp_nma.size() == n_modes*unit_len*n_modes );
//        MMTBX_ASSERT( weights.size() == unit_len );
        std::size_t n_nmpars = x.size();
//        cout << "n_nmpars: " << n_nmpars << endl;
        if(zero_mode_flag){
            MMTBX_ASSERT( n_nmpars == 21 + (n_modes - 5)*(n_modes - 6)/2 );
        }
        else{
            MMTBX_ASSERT( n_nmpars == (n_modes + 1)*n_modes/2 );
        }
        s = unpack_x(x, n_modes, zero_mode_flag);
        gNM.resize(n_nmpars, 0.0);
//        nm_adp = init_nm_adp(modes1d, weights, n_modes, zero_mode_flag);
//        cout << "init nm adp successfully!" << endl;
        for(std::size_t i = 0; i < unit_len; i++){
            d_uaniso_d_nm d_uaniso_d_nm_manager(adp_nma, s, i, n_modes, unit_len, zero_mode_flag);
            af::shared<sym_mat3<double> > d_u_d_nm = d_uaniso_d_nm_manager.d_u_d_nm();
            MMTBX_ASSERT(d_u_d_nm.size() == n_nmpars);
//            cout << "generate derivates successfully!" << endl;
            for(std::size_t k = 0; k < n_nmpars; k++){
                for(std::size_t m = 0; m < 6; m++){
                    gNM[k] += d_target_d_uaniso[i][m]*d_u_d_nm[k][m];
                }
            }
        }
    }

    af::shared<double> grad_nm() { return gNM; }
private:
    af::shared<double> gNM;
    af::versa<double, af::c_grid<2> > s;
//    af::shared<sym_mat3<double> > nm_adp;
};

double grad_resi(af::shared<double> const& gradients){
    double resi = 0.0;
    for(std::size_t i = 0; i < gradients.size(); i++){
        resi += gradients[i]*gradients[i];
    }
    return resi;
}

af::shared<double> add_grads(af::shared<double> const& gradients1,
                     af::shared<double> const& gradients2,
                     double weight){
    MMTBX_ASSERT(gradients1.size()==gradients2.size());
    af::shared<double> gradients(gradients1.size(), af::init_functor_null<double>());
    for(std::size_t i = 0; i < gradients.size(); i++){
        gradients[i] = gradients1[i] + gradients2[i]*weight;
    }
    return gradients;
}

class nm_from_uaniso_target_and_grads {
public:
    nm_from_uaniso_target_and_grads(
                                   af::shared<double> const& x,
                                   af::shared<sym_mat3<double> > const& adp_nma,
                                   af::shared<sym_mat3<double> > const& uanisos,
                                   std::size_t n_modes,
                                   bool zero_mode_flag)
    {
        tg = 0.0;
        uaniso_from_s uaniso_from_s_manager(x, adp_nma, n_modes, zero_mode_flag);
        unm = uaniso_from_s_manager.u_cart();
        for(std::size_t i=0; i < uanisos.size(); i++){
            sym_mat3<double> diff = unm[i] - uanisos[i];
            for(std::size_t k=0; k < diff.size(); k++){
                tg += diff[k]*diff[k];
            }
            diffs.push_back(diff*2.);
        }
        d_target_d_nm d_target_d_nm_manager(adp_nma, diffs, x, n_modes, zero_mode_flag);
        gNM = d_target_d_nm_manager.grad_nm();
    }
    double target() const { return tg; }
    af::shared<double> grad_nm() { return gNM; }
private:
    double tg;
    af::shared<double> gNM;
    af::shared<sym_mat3<double> > diffs;
    af::shared<sym_mat3<double> > unm;
};

}}//namespace mmtbx::nm
#endif //MMTBX_NM_H
