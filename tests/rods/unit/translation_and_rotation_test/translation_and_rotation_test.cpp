#include "rod_structure.h"
#include <string>
#include <cmath>

int main(){
    rod::Rod test_rod("realistic_rod.rodtraj", 0);
    test_rod.load_header("realistic_rod.rodtraj");
    test_rod.load_contents("realistic_rod.rodtraj");
    test_rod.set_units();
    std::array<float, 30> ref_current_r;
    for (int i=0; i<30; i++){
        ref_current_r[i] = test_rod.current_r[i];
    }
    rod::float3 euler_angles = {3.141592654*2, 3.141592654*2, 3.141592654*2};
    test_rod.rotate_rod(euler_angles);
    for (int i=0; i<30; i++){
        if (std::abs(ref_current_r[i] - test_rod.current_r[i]) > 0.01){
            return 1;
        }
    }
    return 0;
} 
