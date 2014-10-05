#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include <math.h>

#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SecondOrderFunctions.h"

class Face
{
	public:
		Face()
		{
			n[0] = NULL;
			n[1] = NULL;
			n[2] = NULL;
			e = NULL;
			vdw_interaction_type = -1;
			area_0 = 0;
			zero_force();
			vdw_xz_interaction_flag = false;
			daddy_blob = NULL;
		}

		~Face()
		{
			n[0] = NULL;
			n[1] = NULL;
			n[2] = NULL;
			e = NULL;
			vdw_interaction_type = -1;
			area_0 = 0;
			zero_force();
			vdw_xz_interaction_flag = false;
			daddy_blob = NULL;
		}

		void init(tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, SecondOrderFunctions::stu centroid_stu, Blob *daddy_blob)
		{
			this->e = e;
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;

			calc_area_normal_centroid();
			area_0 = area;

			this->centroid_stu.s = centroid_stu.s;
			this->centroid_stu.t = centroid_stu.t;
			this->centroid_stu.u = centroid_stu.u;

			this->daddy_blob = daddy_blob;
		}

		void init(mesh_node *n0, mesh_node *n1, mesh_node *n2, Blob *daddy_blob)
		{
			this->e = NULL;
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;

			calc_area_normal_centroid();
			area_0 = area;

			this->centroid_stu.s = 0;
			this->centroid_stu.t = 0;
			this->centroid_stu.u = 0;

			this->daddy_blob = daddy_blob;
		}

		void set_vdw_interaction_type(int vdw_interaction_type)
		{
			this->vdw_interaction_type = vdw_interaction_type;
		}

		void calc_area_normal_centroid()
		{
			// (1/2) * |a x b|
			vector3 a = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
				b = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
			normal.x = a.y * b.z - a.z * b.y;
			normal.y = a.z * b.x - a.x * b.z;
			normal.z = a.x * b.y - a.y * b.x;

			scalar normal_mag = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

			area = .5 * normal_mag;

			// Normalise normal
			normal.x /= normal_mag;
			normal.y /= normal_mag;
			normal.z /= normal_mag;

			// Find centroid
			centroid.x = (1.0/3.0) * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x);
			centroid.y = (1.0/3.0) * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y);
			centroid.z = (1.0/3.0) * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z);
		}

		/* Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3 */
		void barycentric_calc_point(scalar b1, scalar b2, scalar b3, vector3 *p)
		{
			p->x = b1 * n[0]->pos.x + b2 * n[1]->pos.x + b3 * n[2]->pos.x;
			p->y = b1 * n[0]->pos.y + b2 * n[1]->pos.y + b3 * n[2]->pos.y;
			p->z = b1 * n[0]->pos.z + b2 * n[1]->pos.z + b3 * n[2]->pos.z;
		}

		/* Returns the average electrostatic potential of this face */
		scalar average_phi()
		{
			return (1.0/3.0) * (n[0]->phi + n[1]->phi + n[2]->phi);
		}

		scalar get_normal_flux()
		{
			vector3 dphi;
			e->get_grad_phi_at_stu(&dphi, centroid_stu.s, centroid_stu.t, centroid_stu.u);
			return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
		}

/*
                scalar get_normal_flux()
                {
                       	vector3 dphi =  {
                                        e->n[0]->phi * e->dpsi[0] + e->n[1]->phi * e->dpsi[1] + e->n[2]->phi * e->dpsi[2] + e->n[3]->phi * e->dpsi[3],
                                        e->n[0]->phi * e->dpsi[4] + e->n[1]->phi * e->dpsi[5] + e->n[2]->phi * e->dpsi[6] + e->n[3]->phi * e->dpsi[7],
                                       	e->n[0]->phi * e->dpsi[8] + e->n[1]->phi * e->dpsi[9] + e->n[2]->phi * e->dpsi[10] + e->n[3]->phi * e->dpsi[11]
                                        };
                        return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
                }
*/

//		void apply_normal_force(scalar i)
//		{
//			normal_force.x = i * normal.x;
//			normal_force.y = i * normal.y;
//			normal_force.z = i * normal.z;
//		}

//		void zero_normal_force()
//		{
//			normal_force.x = 0;
//			normal_force.y = 0;
//			normal_force.z = 0;
//		}

		void add_force_to_node(int i, vector3 *f)
		{
			force[i].x += f->x;
			force[i].y += f->y;
			force[i].z += f->z;
		}

		void add_force_to_node_atomic(int i, vector3 *f)
		{
			force[i].x += f->x;
			force[i].y += f->y;
			force[i].z += f->z;
		}

		void zero_force()
		{
			for(int i = 0; i < 3; i++) {
				force[i].x = 0;
				force[i].y = 0;
				force[i].z = 0;
			}
		}

		void set_vdw_xz_interaction_flag(bool state)
		{
			vdw_xz_interaction_flag = state;
		}
		
		/*scalar get_vdw_energy_with_face(Face *f2, scalar vdw_r_eq, scalar vdw_eps) 
		{
			if(dot(normal, &f2->normal) < 0) {
				return 0;
			}
		
			const int num_tri_gauss_quad_points = 3;

			struct tri_gauss_point
			{
				scalar W;
				scalar eta[3];
			};

			const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] =
				{
					// Weight, eta1, eta2, eta3
				        {0.333333333333333,     {0.666666666666667, 0.166666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.666666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.166666666666667, 0.666666666666667}}
				};

			vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_points];
	
			// Convert all area coordinate gauss points to cartesian
			for(int i = 0; i < num_tri_gauss_quad_points; i++) {
				barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
				f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);		
			}

			// Get the energy by integrating U(p, q), where p and q are the gauss points on each face
			scalar vdw_E = 0;
			scalar ApAq = area * f2->area;

			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
				for(int l = k; l < num_tri_gauss_quad_points; l++) {
					vector3 r;
					r.x = p[k].x - q[l].x;
					r.y = p[k].y - q[l].y;
					r.z = p[k].z - q[l].z;
						
					scalar mag_r = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);

					vdw_E += pow(vdw_r_eq, 6) * vdw_eps * (pow(vdw_r_eq, 6) * pow(mag_r, -12) - 2 * pow(mag_r, -6));
				}
			}

			vdw_E *= ApAq;

			return vdw_E;
		}*/

		scalar get_vdw_energy_with_xz_plane(scalar dim_y, scalar vdw_r_eq, scalar vdw_eps)
		{
			// Check which wall face is facing
			scalar y_wall;
			if(normal.y < 0) {
				y_wall = 0;
			} else {
				y_wall = dim_y;
			}

			const int num_tri_gauss_quad_points = 3;

			struct tri_gauss_point
			{
				scalar W;
				scalar eta[3];
			};

			const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] =
				{
					// Weight, eta1, eta2, eta3
				        {0.333333333333333,     {0.666666666666667, 0.166666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.666666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.166666666666667, 0.666666666666667}}
				};

			vector3 p[num_tri_gauss_quad_points];

			// Convert all area coordinate gauss points to cartesian
			for(int i = 0; i < num_tri_gauss_quad_points; i++) {
				barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
			}

			// Get the energy by integrating U(p, q), where p and q are the gauss points on each face
			scalar vdw_E = 0;
			scalar Asq = area * area;
			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
				for(int l = k; l < num_tri_gauss_quad_points; l++) {
					scalar mag_r = p[k].y - y_wall;

					vdw_E += pow(vdw_r_eq, 6) * vdw_eps * (pow(vdw_r_eq, 6) * pow(mag_r, -12) - 2 * pow(mag_r, -6));
				}
			}
			vdw_E *= Asq;

			return vdw_E;
		}

		bool is_vdw_active()
		{
			if(vdw_interaction_type == -1) {
				return false;
			} else {
				return true;
			}
		}

		/* Pointers to the surface nodes that define this face */
		mesh_node *n[3];

		/* Pointer to the element this is a face of */
		tetra_element_linear *e;

		/* Van der Waals interaction type */
		int vdw_interaction_type;

		/* Initial, equilibrium area of this face */
		scalar area_0;

		/* Last calculated area (by most recent call to calc_area_normal_centroid()) */
		scalar area;

		/* Last calculated normal (by most recent call to calc_area_normal_centroid()) */
		vector3 normal;

		/* Last calculated centroid (by most recent call to calc_area_normal_centroid()) */
		vector3 centroid;

		/* Stores the current force applied to this face */
		vector3 force[3];

		
		/* Stores the natural (shape function) coords of the centroid of this face in the parent element */
		SecondOrderFunctions::stu centroid_stu;

		bool vdw_xz_interaction_flag;

		Blob *daddy_blob;
};

#endif
