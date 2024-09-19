#include "rod_math_v9.h"

int main(){
    rod::float9 nodes = {0,0,0, 1,0,0, 2,0,0};
    rod::float3 current_x = {nodes[0], nodes[1], nodes[2]};
    rod::float3 current_xp1 = {nodes[3], nodes[4], nodes[5]};
    rod::float3 current_xp2 = {nodes[6], nodes[7], nodes[8]};
    rod::float3 e_i_equil;
    rod::float3 e_ip1_equil;
    rod::get_p_i(current_x, current_xp1, e_i_equil);
    rod::get_p_i(current_xp1, current_xp2, e_ip1_equil);
    
    rod::float3 m_i = {0,1,0};
    rod::float3 m_ip1 = {0,1,0};
    
    nodes[5] += 0.5;
    
    rod::float3 new_x = {nodes[0], nodes[1], nodes[2]};
    rod::float3 new_xp1 = {nodes[3], nodes[4], nodes[5]};
    rod::float3 new_xp2 = {nodes[6], nodes[7], nodes[8]};
    rod::float3 e_i;
    rod::float3 e_ip1;
    rod::get_p_i(new_x, new_xp1, e_i);
    rod::get_p_i(new_xp1, new_xp2, e_ip1);
    
    rod::float3 m_i_prime;
    rod::float3 m_ip1_prime;
    
    rod::float3 t_i;
    rod::float3 t_ip1;
    rod::float3 t_i_equil;
    rod::float3 t_ip1_equil;
    rod::normalize(e_i, t_i);
    rod::normalize(e_i_equil, t_i_equil);
    rod::normalize(e_ip1_equil, t_ip1_equil);
    rod::normalize(e_ip1, t_ip1);
    
    rod::update_m1_matrix( m_i, e_i_equil, e_i, m_i_prime );
    rod::update_m1_matrix( m_ip1, e_ip1_equil, e_ip1, m_ip1_prime );
       
    float twist = rod::get_twist_energy(1, m_ip1_prime, m_i_prime, m_ip1, m_i, e_i, e_ip1, e_i_equil, e_ip1_equil);
    
    if (twist<0.001){
        return 0;
    }
    //std::cout << twist << "\n";
    return 1;
    // test: doing this material frame update (check on wolfram?) SHOULD NOT MAKE A TWIST ENERGY    
}
