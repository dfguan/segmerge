/*
 * =====================================================================================
 *
 *       Filename:  opt.h
 *
 *    Description:  options 
 *
 *        Version:  1.0
 *        Created:  02/05/2018 19:47:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _OPT_H
#define _OPT_H


#include <stdint.h>

typedef struct {
	//int32_t max_oh; //maximum overhanged length 
	//float ratio; //minimum mapped ratio
	//char *gfa_fn; // not really alloated
	char *paf_fn; //not really allocated
	//char *out;//output file
	//int32_t min_ovlp:31, fmt:1;//minimum overlap
	int32_t min_bl;
	int32_t max_gs;
}opt;
int parse_args(int argc, char *argv[], opt *o);
#endif
