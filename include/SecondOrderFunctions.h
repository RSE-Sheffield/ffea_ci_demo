#ifndef SECONDORDERFUNCTIONS_H_INCLUDED
#define SECONDORDERFUNCTIONS_H_INCLUDED

#define GRAD_PSI_1 0
#define GRAD_PSI_2 1
#define GRAD_PSI_3 2
#define GRAD_PSI_4 3
#define GRAD_PSI_5 4
#define GRAD_PSI_6 5
#define GRAD_PSI_7 6
#define GRAD_PSI_8 7
#define GRAD_PSI_9 8
#define GRAD_PSI_10 9

#define DS_BY_DX 0
#define DT_BY_DX 1
#define DU_BY_DX 2
#define DS_BY_DY 3
#define DT_BY_DY 4
#define DU_BY_DY 5
#define DS_BY_DZ 6
#define DT_BY_DZ 7
#define DU_BY_DZ 8


#include "mesh_node.h"

namespace SecondOrderFunctions {

    struct stu {
        scalar s, t, u;
    };

    struct abcd {
        scalar a, b, c, d;
    };

     static stu stu_lookup[10] ={
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {.5, 0, 0},
        {0, .5, 0},
        {0, 0, .5},
        {.5, .5, 0},
        {.5, 0, .5},
        {0, .5, .5}
    }; 

    static void calc_psi(scalar psi[10], scalar s, scalar t, scalar u) {
        scalar r = 1 - s - t - u;
        psi[0] = (2 * r - 1) * r;
        psi[1] = (2 * s - 1) * s;
        psi[2] = (2 * t - 1) * t;
        psi[3] = (2 * u - 1) * u;
        psi[4] = 4 * r * s;
        psi[5] = 4 * r * t;
        psi[6] = 4 * r * u;
        psi[7] = 4 * s * t;
        psi[8] = 4 * s * u;
        psi[9] = 4 * t * u;
    } 

     static void calc_grad_psi(vector3 grad_psi[10], scalar s, scalar t, scalar u, scalar J_inv[9]) {
        scalar dpsi1_by_dstu = 4 * (s + t + u) - 3;
        grad_psi[GRAD_PSI_1].x = dpsi1_by_dstu * (J_inv[DS_BY_DX] + J_inv[DT_BY_DX] + J_inv[DU_BY_DX]);
        grad_psi[GRAD_PSI_1].y = dpsi1_by_dstu * (J_inv[DS_BY_DY] + J_inv[DT_BY_DY] + J_inv[DU_BY_DY]);
        grad_psi[GRAD_PSI_1].z = dpsi1_by_dstu * (J_inv[DS_BY_DZ] + J_inv[DT_BY_DZ] + J_inv[DU_BY_DZ]);

        scalar dpsi2_by_ds = 4 * s - 1;
        grad_psi[GRAD_PSI_2].x = dpsi2_by_ds * J_inv[DS_BY_DX];
        grad_psi[GRAD_PSI_2].y = dpsi2_by_ds * J_inv[DS_BY_DY];
        grad_psi[GRAD_PSI_2].z = dpsi2_by_ds * J_inv[DS_BY_DZ];

        scalar dpsi3_by_dt = 4 * t - 1;
        grad_psi[GRAD_PSI_3].x = dpsi3_by_dt * J_inv[DT_BY_DX];
        grad_psi[GRAD_PSI_3].y = dpsi3_by_dt * J_inv[DT_BY_DY];
        grad_psi[GRAD_PSI_3].z = dpsi3_by_dt * J_inv[DT_BY_DZ];

        scalar dpsi4_by_du = 4 * u - 1;
        grad_psi[GRAD_PSI_4].x = dpsi4_by_du * J_inv[DU_BY_DX];
        grad_psi[GRAD_PSI_4].y = dpsi4_by_du * J_inv[DU_BY_DY];
        grad_psi[GRAD_PSI_4].z = dpsi4_by_du * J_inv[DU_BY_DZ];

        scalar dpsi5_by_ds = 4 * (1 - 2 * s - t - u), dpsi5_by_dtu = -4 * s;
        grad_psi[GRAD_PSI_5].x = dpsi5_by_ds * J_inv[DS_BY_DX] + dpsi5_by_dtu * (J_inv[DT_BY_DX] + J_inv[DU_BY_DX]);
        grad_psi[GRAD_PSI_5].y = dpsi5_by_ds * J_inv[DS_BY_DY] + dpsi5_by_dtu * (J_inv[DT_BY_DY] + J_inv[DU_BY_DY]);
        grad_psi[GRAD_PSI_5].z = dpsi5_by_ds * J_inv[DS_BY_DZ] + dpsi5_by_dtu * (J_inv[DT_BY_DZ] + J_inv[DU_BY_DZ]);

        scalar dpsi6_by_dt = 4 * (1 - s - 2 * t - u), dpsi6_by_dsu = -4 * t;
        grad_psi[GRAD_PSI_6].x = dpsi6_by_dt * J_inv[DT_BY_DX] + dpsi6_by_dsu * (J_inv[DS_BY_DX] + J_inv[DU_BY_DX]);
        grad_psi[GRAD_PSI_6].y = dpsi6_by_dt * J_inv[DT_BY_DY] + dpsi6_by_dsu * (J_inv[DS_BY_DY] + J_inv[DU_BY_DY]);
        grad_psi[GRAD_PSI_6].z = dpsi6_by_dt * J_inv[DT_BY_DZ] + dpsi6_by_dsu * (J_inv[DS_BY_DZ] + J_inv[DU_BY_DZ]);

        scalar dpsi7_by_du = 4 * (1 - s - t - 2 * u), dpsi7_by_dst = -4 * u;
        grad_psi[GRAD_PSI_7].x = dpsi7_by_du * J_inv[DU_BY_DX] + dpsi7_by_dst * (J_inv[DS_BY_DX] + J_inv[DT_BY_DX]);
        grad_psi[GRAD_PSI_7].y = dpsi7_by_du * J_inv[DU_BY_DY] + dpsi7_by_dst * (J_inv[DS_BY_DY] + J_inv[DT_BY_DY]);
        grad_psi[GRAD_PSI_7].z = dpsi7_by_du * J_inv[DU_BY_DZ] + dpsi7_by_dst * (J_inv[DS_BY_DZ] + J_inv[DT_BY_DZ]);

        scalar dpsi8_by_ds = 4 * t, dpsi8_by_dt = 4 * s;
        grad_psi[GRAD_PSI_8].x = dpsi8_by_ds * J_inv[DS_BY_DX] + dpsi8_by_dt * J_inv[DT_BY_DX];
        grad_psi[GRAD_PSI_8].y = dpsi8_by_ds * J_inv[DS_BY_DY] + dpsi8_by_dt * J_inv[DT_BY_DY];
        grad_psi[GRAD_PSI_8].z = dpsi8_by_ds * J_inv[DS_BY_DZ] + dpsi8_by_dt * J_inv[DT_BY_DZ];

        scalar dpsi9_by_ds = 4 * u, dpsi9_by_du = 4 * s;
        grad_psi[GRAD_PSI_9].x = dpsi9_by_ds * J_inv[DS_BY_DX] + dpsi9_by_du * J_inv[DU_BY_DX];
        grad_psi[GRAD_PSI_9].y = dpsi9_by_ds * J_inv[DS_BY_DY] + dpsi9_by_du * J_inv[DU_BY_DY];
        grad_psi[GRAD_PSI_9].z = dpsi9_by_ds * J_inv[DS_BY_DZ] + dpsi9_by_du * J_inv[DU_BY_DZ];

        scalar dpsi10_by_dt = 4 * u, dpsi10_by_du = 4 * t;
        grad_psi[GRAD_PSI_10].x = dpsi10_by_dt * J_inv[DT_BY_DX] + dpsi10_by_du * J_inv[DU_BY_DX];
        grad_psi[GRAD_PSI_10].y = dpsi10_by_dt * J_inv[DT_BY_DY] + dpsi10_by_du * J_inv[DU_BY_DY];
        grad_psi[GRAD_PSI_10].z = dpsi10_by_dt * J_inv[DT_BY_DZ] + dpsi10_by_du * J_inv[DU_BY_DZ];

        //		for(int i = 0; i < 10; i++) {
        //			printf("dpsi%d/dx = %+e\n", i + 1, grad_psi[i].x);
        //			printf("dpsi%d/dy = %+e\n", i + 1, grad_psi[i].y);
        //			printf("dpsi%d/dz = %+e\n", i + 1, grad_psi[i].z);
        //		}
    } 

  /** Construct jacobian column coefficients */ 
  static void calc_jacobian_column_coefficients(mesh_node *n[10], abcd J_coeff[3][3]) {
        // dx/ds
        J_coeff[0][0].a = 4 * (n[0]->pos.x + n[1]->pos.x - 2 * n[4]->pos.x);
        J_coeff[0][0].b = 4 * (n[0]->pos.x - n[4]->pos.x - n[5]->pos.x + n[7]->pos.x);
        J_coeff[0][0].c = 4 * (n[0]->pos.x - n[4]->pos.x - n[6]->pos.x + n[8]->pos.x);
        J_coeff[0][0].d = 4 * n[4]->pos.x - 3 * n[0]->pos.x - n[1]->pos.x;

        // dy/ds
        J_coeff[1][0].a = 4 * (n[0]->pos.y + n[1]->pos.y - 2 * n[4]->pos.y);
        J_coeff[1][0].b = 4 * (n[0]->pos.y - n[4]->pos.y - n[5]->pos.y + n[7]->pos.y);
        J_coeff[1][0].c = 4 * (n[0]->pos.y - n[4]->pos.y - n[6]->pos.y + n[8]->pos.y);
        J_coeff[1][0].d = 4 * n[4]->pos.y - 3 * n[0]->pos.y - n[1]->pos.y;

        // dz/ds
        J_coeff[2][0].a = 4 * (n[0]->pos.z + n[1]->pos.z - 2 * n[4]->pos.z);
        J_coeff[2][0].b = 4 * (n[0]->pos.z - n[4]->pos.z - n[5]->pos.z + n[7]->pos.z);
        J_coeff[2][0].c = 4 * (n[0]->pos.z - n[4]->pos.z - n[6]->pos.z + n[8]->pos.z);
        J_coeff[2][0].d = 4 * n[4]->pos.z - 3 * n[0]->pos.z - n[1]->pos.z;


        // dx/dt
        J_coeff[0][1].a = 4 * (n[0]->pos.x - n[4]->pos.x - n[5]->pos.x + n[7]->pos.x);
        J_coeff[0][1].b = 4 * (n[0]->pos.x + n[2]->pos.x - 2 * n[5]->pos.x);
        J_coeff[0][1].c = 4 * (n[0]->pos.x - n[5]->pos.x - n[6]->pos.x + n[9]->pos.x);
        J_coeff[0][1].d = 4 * n[5]->pos.x - 3 * n[0]->pos.x - n[2]->pos.x;

        // dy/dt
        J_coeff[1][1].a = 4 * (n[0]->pos.y - n[4]->pos.y - n[5]->pos.y + n[7]->pos.y);
        J_coeff[1][1].b = 4 * (n[0]->pos.y + n[2]->pos.y - 2 * n[5]->pos.y);
        J_coeff[1][1].c = 4 * (n[0]->pos.y - n[5]->pos.y - n[6]->pos.y + n[9]->pos.y);
        J_coeff[1][1].d = 4 * n[5]->pos.y - 3 * n[0]->pos.y - n[2]->pos.y;

        // dz/dt
        J_coeff[2][1].a = 4 * (n[0]->pos.z - n[4]->pos.z - n[5]->pos.z + n[7]->pos.z);
        J_coeff[2][1].b = 4 * (n[0]->pos.z + n[2]->pos.z - 2 * n[5]->pos.z);
        J_coeff[2][1].c = 4 * (n[0]->pos.z - n[5]->pos.z - n[6]->pos.z + n[9]->pos.z);
        J_coeff[2][1].d = 4 * n[5]->pos.z - 3 * n[0]->pos.z - n[2]->pos.z;


        // dx/du
        J_coeff[0][2].a = 4 * (n[0]->pos.x - n[4]->pos.x - n[6]->pos.x + n[8]->pos.x);
        J_coeff[0][2].b = 4 * (n[0]->pos.x - n[5]->pos.x - n[6]->pos.x + n[9]->pos.x);
        J_coeff[0][2].c = 4 * (n[0]->pos.x + n[3]->pos.x - 2 * n[6]->pos.x);
        J_coeff[0][2].d = 4 * n[6]->pos.x - 3 * n[0]->pos.x - n[3]->pos.x;

        // dy/du
        J_coeff[1][2].a = 4 * (n[0]->pos.y - n[4]->pos.y - n[6]->pos.y + n[8]->pos.y);
        J_coeff[1][2].b = 4 * (n[0]->pos.y - n[5]->pos.y - n[6]->pos.y + n[9]->pos.y);
        J_coeff[1][2].c = 4 * (n[0]->pos.y + n[3]->pos.y - 2 * n[6]->pos.y);
        J_coeff[1][2].d = 4 * n[6]->pos.y - 3 * n[0]->pos.y - n[3]->pos.y;

        // dz/du
        J_coeff[2][2].a = 4 * (n[0]->pos.z - n[4]->pos.z - n[6]->pos.z + n[8]->pos.z);
        J_coeff[2][2].b = 4 * (n[0]->pos.z - n[5]->pos.z - n[6]->pos.z + n[9]->pos.z);
        J_coeff[2][2].c = 4 * (n[0]->pos.z + n[3]->pos.z - 2 * n[6]->pos.z);
        J_coeff[2][2].d = 4 * n[6]->pos.z - 3 * n[0]->pos.z - n[3]->pos.z;
    }

   static scalar calc_det_J(abcd J_coeff[3][3], scalar s, scalar t, scalar u, scalar J_inv[9]) {
        scalar J[3][3];

        // Construct Jacobian for the part of the element at (s, t, u)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                J[i][j] = J_coeff[i][j].a * s + J_coeff[i][j].b * t + J_coeff[i][j].c * u + J_coeff[i][j].d;
            }
        }

        // Invert the Jacobian and transpose it to obtain the derivatives of local coordinates wrt cartesian coords
        J_inv[DS_BY_DX] = J[2][2] * J[1][1] - J[2][1] * J[1][2];
        J_inv[DS_BY_DY] = J[2][1] * J[0][2] - J[2][2] * J[0][1];
        J_inv[DS_BY_DZ] = J[1][2] * J[0][1] - J[1][1] * J[0][2];

        J_inv[DT_BY_DX] = J[2][0] * J[1][2] - J[2][2] * J[1][0];
        J_inv[DT_BY_DY] = J[2][2] * J[0][0] - J[2][0] * J[0][2];
        J_inv[DT_BY_DZ] = J[1][0] * J[0][2] - J[1][2] * J[0][0];

        J_inv[DU_BY_DX] = J[2][1] * J[1][0] - J[2][0] * J[1][1];
        J_inv[DU_BY_DY] = J[2][0] * J[0][1] - J[2][1] * J[0][0];
        J_inv[DU_BY_DZ] = J[1][1] * J[0][0] - J[1][0] * J[0][1];

        // Calculate determinant
        scalar det = J[0][0] * J_inv[DS_BY_DX] + J[1][0] * J_inv[DS_BY_DY] + J[2][0] * J_inv[DS_BY_DZ];

        for (int i = 0; i < 9; i++) {
            J_inv[i] /= det;
        }

        
         //               printf("J:\n");
         //               printf("%+e %+e %+e\n", J[0][0], J[0][1], J[0][2]);
         //               printf("%+e %+e %+e\n", J[1][0], J[1][1], J[1][2]);
         //               printf("%+e %+e %+e\n", J[2][0], J[2][1], J[2][2]);

         //               printf("J_inv:\n");
         //               printf("%+e %+e %+e\n", J_inv[0], J_inv[1], J_inv[2]);
         //               printf("%+e %+e %+e\n", J_inv[3], J_inv[4], J_inv[5]);
         //               printf("%+e %+e %+e\n", J_inv[6], J_inv[7], J_inv[8]);
         

        return fabs(det);
    }

}

#endif
