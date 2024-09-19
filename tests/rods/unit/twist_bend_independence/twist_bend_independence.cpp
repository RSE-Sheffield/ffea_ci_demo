#include "rod_math_v9.h"

int main(){
    rod::float3 m_i = {0,1,0};
    rod::float3 m_im1 = {0,1,0};
    rod::float3 m_i_equil = {0,1,0};
    rod::float3 m_im1_equil = {0,1,0};
    rod::float3 e_im1 = {1,0,0};
    rod::float3 e_i = {1,0,0};
    rod::float3 e_im1_equil = {1,0,0};
    rod::float3 e_i_equil = {1,0,0};
    rod::rodrigues_rotation(m_i, e_i, 4.2, m_i);
    float twist = rod::get_twist_energy(1, m_i, m_im1, m_i_equil, m_im1_equil, e_im1, e_i, e_im1_equil, e_i_equil);
    
    rod::float3 m_i_2;
    rod::cross_product(m_i, e_i, m_i_2);

    rod::float3 m_i_2_equil;
    rod::cross_product(m_i_equil, e_i_equil, m_i_2_equil);
    
    rod::float3 m_im1_2;
    rod::cross_product(m_im1, e_im1, m_im1_2);

    rod::float3 m_im1_2_equil;
    rod::cross_product(m_im1_equil, e_im1_equil, m_im1_2_equil);
    
    rod::float4 B_equil = {1,0,0,1};
    
    float bend = rod::get_bend_energy_from_p(e_im1, e_i, e_im1_equil, e_i_equil, m_im1_2, m_im1, m_im1_2_equil,
        m_im1_equil,
        m_i_2,
        m_i,
        m_i_2_equil,
        m_i_equil,
        B_equil,
        B_equil);

    if (twist > 0 && bend > -0.001 && bend < 0.001){
        return 0;
    }

//    std::cout << "m_i rotated = " << m_i[0] << ", " << m_i[1] << ", " << m_i[2] << "\n";
//    std::cout << "twist emergy (should be >0)" << twist << "\n";
//    std::cout << "bend energy (should be ~0)" << bend << "\n";

    return 1;
} 
