/*
 * d3db.h
 *
 *  Created on: Dec 9, 2013
 *      Author: bylaska
 */

#ifndef D3DB_H_
#define D3DB_H_

extern int d3db_nfft3d();
extern int d3db_nfft3d_map();
extern int d3db_n2ft3d();
extern int d3db_n2ft3d_map();
extern int d3db_nq();
extern int d3db_nx();
extern int d3db_ny();
extern int d3db_nz();
extern int d3db_zplane_size();

extern void d3db_end();
extern void d3db_init(int, int, int, int);
extern void d3db_cr_fft3b(REAL*, REAL*, REAL*);

#endif /* D3DB_H_ */
