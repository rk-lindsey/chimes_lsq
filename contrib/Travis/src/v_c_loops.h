// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file c_loops.hh
 * \brief Header file for the loop classes. */

#ifndef VOROPP_C_LOOPS_HH
#define VOROPP_C_LOOPS_HH

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "v_config.h"


/** \brief Base class for looping over particles in a container.
 *
 * This class forms the base of all classes that can loop over a subset of
 * particles in a contaner in some order. When initialized, it stores constants
 * about the corresponding container geometry. It also contains a number of
 * routines for interrogating which particle currently being considered by the
 * loop, which are common between all of the derived classes. */
class c_loop_base {
	public:
		/** The number of blocks in the x direction. */
		const int nx;
		/** The number of blocks in the y direction. */
		const int ny;
		/** The number of blocks in the z direction. */
		const int nz;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines that step through blocks in
		 * sequence. */
		const int nxy;
		/** A constant, set to the value of nx*ny*nz, which is used in
		 * the routines that step through blocks in sequence. */
		const int nxyz;
		/** The number of floating point numbers per particle in the
		 * associated container data structure. */
		const int ps;
		/** A pointer to the particle position information in the
		 * associated container data structure. */
		double **p;
		/** A pointer to the particle ID information in the associated
		 * container data structure. */
		int **id;
		/** A pointer to the particle counts in the associated
		 * container data structure. */
		int *co;
		/** The current x-index of the block under consideration by the
		 * loop. */
		int i;
		/** The current y-index of the block under consideration by the
		 * loop. */
		int j;
		/** The current z-index of the block under consideration by the
		 * loop. */
		int k;
		/** The current index of the block under consideration by the
		 * loop. */
		int ijk;
		/** The index of the particle under consideration within the current
		 * block. */
		int q;
		/** The constructor copies several necessary constants from the
		 * base container class.
		 * \param[in] con the container class to use. */
		template<class c_class>
		c_loop_base(c_class &con) : nx(con.nx), ny(con.ny), nz(con.nz),
					    nxy(con.nxy), nxyz(con.nxyz), ps(con.ps),
					    p(con.p), id(con.id), co(con.co) {}
		/** Returns the position vector of the particle currently being
		 * considered by the loop.
		 * \param[out] (x,y,z) the position vector of the particle. */
		inline void pos(double &x,double &y,double &z) {
			double *pp=p[ijk]+ps*q;
			x=*(pp++);y=*(pp++);z=*pp;
		}
		/** Returns the ID, position vector, and radius of the particle
		 * currently being considered by the loop.
		 * \param[out] pid the particle ID.
		 * \param[out] (x,y,z) the position vector of the particle.
		 * \param[out] r the radius of the particle. If no radius
		 * 		 information is available the default radius
		 * 		 value is returned. */
		inline void pos(int &pid,double &x,double &y,double &z,double &r) {
			pid=id[ijk][q];
			double *pp=p[ijk]+ps*q;
			x=*(pp++);y=*(pp++);z=*pp;
			r=ps==3?default_radius:*(++pp);
		}
		/** Returns the x position of the particle currently being
		 * considered by the loop. */
		inline double x() {return p[ijk][ps*q];}
		/** Returns the y position of the particle currently being
		 * considered by the loop. */
		inline double y() {return p[ijk][ps*q+1];}
		/** Returns the z position of the particle currently being
		 * considered by the loop. */
		inline double z() {return p[ijk][ps*q+2];}
		/** Returns the ID of the particle currently being considered
		 * by the loop. */
		inline int pid() {return id[ijk][q];}
};


/** \brief A class for looping over all particles in a container_periodic or
 * container_periodic_poly class.
 *
 * Since the container_periodic and container_periodic_poly classes have a
 * fundamentally different memory organization, the regular loop classes cannot
 * be used with them. */
class c_loop_all_periodic : public c_loop_base {
	public:
		/** The constructor copies several necessary constants from the
		 * base periodic container class.
		 * \param[in] con the periodic container class to use. */
		template<class c_class>
		c_loop_all_periodic(c_class &con) : c_loop_base(con), ey(con.ey), ez(con.ez), wy(con.wy), wz(con.wz),
			ijk0(nx*(ey+con.oy*ez)), inc2(2*nx*con.ey+1) {}
		/** Sets the class to consider the first particle.
		 * \return True if there is any particle to consider, false
		 * otherwise. */
		inline bool start() {
			i=0;
			j=ey;
			k=ez;
			ijk=ijk0;
			q=0;
			while(co[ijk]==0) if(!next_block()) return false;
			return true;
		}
		/** Finds the next particle to test.
		 * \return True if there is another particle, false if no more
		 * particles are available. */
		inline bool inc() {
			q++;
			if(q>=co[ijk]) {
				q=0;
				do {
					if(!next_block()) return false;
				} while(co[ijk]==0);
			}
			return true;
		}
	private:
		/** The lower y index (inclusive) of the primary domain within
		 * the block structure. */
		int ey;
		/** The lower y index (inclusive) of the primary domain within
		 * the block structure. */
		int ez;
		/** The upper y index (exclusive) of the primary domain within
		 * the block structure. */
		int wy;
		/** The upper z index (exclusive) of the primary domain within
		 * the block structure. */
		int wz;
		/** The index of the (0,0,0) block within the block structure.
		 */
		int ijk0;
		/** A value to increase ijk by when the z index is increased.
		 */
		int inc2;
		/** Updates the internal variables to find the next
		 * computational block with any particles.
		 * \return True if another block is found, false if there are
		 * no more blocks. */
		inline bool next_block() {
			i++;
			if(i==nx) {
				i=0;j++;
				if(j==wy) {
					j=ey;k++;
					if(k==wz) return false;
					ijk+=inc2;
				} else ijk++;
			} else ijk++;
			return true;
		}
};


#endif


