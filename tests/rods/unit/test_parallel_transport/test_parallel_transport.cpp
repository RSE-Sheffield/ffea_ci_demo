#include "mat_vec_fns_II.h"
#include "rod_math_v9.h"

int main(){
    rod::float3 m_im1 = {-8.78125100e-11, 1.16555800e-10, -8.71992700e-11};
    const rod::float3 e_im1 = {1.34365400e-08, 6.62241900e-09, -4.68984700e-09};
    const rod::float3 e_i = {6.55867000e-09, 2.03361000e-09, 5.91157700e-09};
    rod::float3 m_prime;
    rod::parallel_transport(m_im1, m_prime, e_im1, e_i);
    const float dotprod = dot(m_prime, e_i);
    if (dotprod < 0.01){
        return 0;
    }
    else{
        return 1;
    }
}
