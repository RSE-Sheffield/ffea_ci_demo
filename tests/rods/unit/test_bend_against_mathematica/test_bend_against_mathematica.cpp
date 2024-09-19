#include "rod_math_v9.h"

int main(){
    float mathematica_result = 059.6893*2;
    rod::float4 bbar = {1,0,0,1};
    rod::float3 eim1 = {0,0, 0.008};
    rod::float3 ei = {0, 0.00841471, 0.00540302};
    rod::float3 eim1bar = {0,0,0.01};
    rod::float3 eibar = {0,0,0.01};
    rod::float3 mim12 = {1,0,0};
    rod::float3 mim1 = {0, 1, 0};
    rod::float3 mim12bar = {1,0,0};
    rod::float3 mim1bar = {0,1,0};
    rod::float3 mi2 = {1,0,0};
    rod::float3 mi = {0, -0.540302, -0.841471};
    rod::float3 mi2bar = {1,0,0};
    rod::float3 mibar = {0,1,0};
    //rod::float3 kbi;
    //rod::float3 kbibar;
    
    float computed_energy = rod::get_bend_energy_from_p(
        eim1,
        ei,
        eim1bar,
        eibar,
        mim12,
        mim1,
        mim12bar,
        mim1bar,
        mi2,
        mi,
        mi2bar,
        mibar,
        bbar,
        bbar);
        
    if (computed_energy > mathematica_result - 0.01 && computed_energy < mathematica_result + 0.01){
        return 0;
    }
    else{
        std::cout << "Computed energy: " << computed_energy << "\n";
        std::cout << "Analytical energy: " << mathematica_result << "\n";
        return 1;
    }
}
