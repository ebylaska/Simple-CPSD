/*
 * Parallel.h
 *
 *  Created on: Dec 9, 2013
 *      Author: bylaska
 */

#ifndef PARALLEL_H_
#define PARALLEL_H_

extern int Parallel_np();
extern int Parallel_taskid();

extern void Parallel_init(int*, char***);
extern void Parallel_end();



#endif /* PARALLEL_H_ */
