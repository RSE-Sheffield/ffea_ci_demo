#include "rod_math_v9.h"

int main(){
    // testing rodriguez
    rod::float3 e_i = {1,0,0};
    rod::float3 m_i = {0,1,0};
    rod::float3 rotated_mi;
    float theta = 1.570796327*2;
    rod::rodrigues_rotation(m_i, e_i, theta, rotated_mi);
    if (rotated_mi[1] < -0.99 && rotated_mi[1] > -1.01 && rotated_mi[2] < 0.001 && rotated_mi[2] > -0.001){
        return 0;
    }
    return 1;
} 
