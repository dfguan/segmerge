/*
 * =====================================================================================
 * *       Filename:  opt.c
 *
 *    Description:  opt.c
 *
 *        Version:  1.0
 *        Created:  02/05/2018 19:51:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "opt.h"

int help()
{
	fprintf(stderr, "\nUsage: update [options] <PAF>\n");
	fprintf(stderr, "Options:\n");	
	/*fprintf(stderr, "         -r    FLOAT    minimum overlap ratio for an alignment [0.8]\n");	*/
	fprintf(stderr, "         -m    INT      minimum alignment block length [3K]\n");
	fprintf(stderr, "         -M    INT      maximum gap size for chaining [20K]\n")	;
	fprintf(stderr, "         -r             read to reference alignment [NO]\n")	;
	/*fprintf(stderr, "         -O    STR      output file format: GFA, FASTA [GFA]\n");*/
	/*fprintf(stderr, "         -o    FILE     output file [stdout]\n");*/

	fprintf(stderr, "         -h             help\n")	;
	return 0;
}

int parse_args(int argc, char *argv[], opt *o)
{
	o->max_gs = 20000;
	o->min_bl = 3000;
	o->s2s = 1;
	int c;
	while ((c = getopt(argc, argv, "m:M:rh")) != -1) {
		switch (c) {
			/*case 'r':*/
				/*o->ratio = strtof(optarg, NULL);*/
				/*break;*/
			case 'M':
				o->max_gs = atoi(optarg);
				break;
			case 'm':
				o->min_bl = atoi(optarg);
				break;
			case 'r':
				o->s2s = 0;
				break;
			/*case 'O':*/
				/*if (!strcmp("FASTA", optarg))*/
					/*o->fmt = 0;*/
				/*break;*/
			/*case 'o':*/
				/*o->out = optarg;*/
			case 'h':
				help();
				return 1;
			default:
				fprintf(stderr,"[E::%s] undefined option %c\n", __func__, c);
				help();
				return 1;
		}
	}	
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] paf file can be omitted!\n", __func__);
		help();
		return 1;
	} 
	o->paf_fn = argv[optind];
	return 0;
}

