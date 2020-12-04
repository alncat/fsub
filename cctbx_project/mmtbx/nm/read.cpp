#include <iostream>
#include <fstream>
using namespace std;
int main(){
    ifstream file ("4ki8_evec.dat", ios::in | ios::binary);
    std::size_t total;
    if (file.is_open()){
        for(std::size_t i = 0; i < 100000; i++){
            char serial[9];
            file.read(serial, 8);
            file.seekg(5, ios::cur);
            char resname[4];
            file.read(resname, 3);
            file.seekg(4, ios::cur);
            char atmname[5];
            file.read(atmname, 4);
            if(serial[0]=='E' and serial[1]=='N' and serial[2]=='D'){
                total = i;
                break;
            }
            if( strncmp(resname, "VAL", 3) ){
//                cout << "good" << endl;
            }
            cout << serial << resname << atmname <<endl;
        }
        for(std::size_t i = 0; i < 50; i++){
            for(std::size_t j = 0; j < total; j++){
                double vx, vy, vz;
                file.read((char *) &vx, 8);
                file.read((char *) &vy, 8);
                file.read((char *) &vz, 8);
                cout << vx << vy << vz<< endl;
            }
            double vx, vy, vz;
            file.read((char *) &vx, 8);
            file.read((char *) &vy, 8);
            file.read((char *) &vz, 8);
            if( vx==111111111. and vy==-222222222. and vz==333333333. )
                cout <<"good"<< endl;
        }
    }
    else cout << "unable to open file";
    return 0;
}
