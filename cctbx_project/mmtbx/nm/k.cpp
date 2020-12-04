#include <iostream>
using namespace std;
int main(){
    for(std::size_t i = 0; i < 6; i++){
        for(std::size_t j =0; j <= i; j++){
            for(int k = i; k < 6; k++){
                cout<< "k =" << k << "i=" << i << "j=" << j;
            }
            cout << endl;
            for(int k = j; k <= i - 1; k++){
                cout<< "k =" << k << "i=" << i << "j=" << j;
            }
            cout << endl;
        }
    }
    return 1;
}
