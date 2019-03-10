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
	fprintf(stderr, "\nUsage: update [options] <GFA> <PAF>\n");
	fprintf(stderr, "Options:\n");	
	fprintf(stderr, "         -r    FLOAT    minimum overlap ratio for an alignment [0.8]\n");	
	fprintf(stderr, "         -m    INT      maximum overhangs [1000]\n");
	fprintf(stderr, "         -M    INT      minimum mapped length [5000]\n")	;
	fprintf(stderr, "         -O    STR      output file format: GFA, FASTA [GFA]\n");
	fprintf(stderr, "         -o    FILE     output file [stdout]\n");

	fprintf(stderr, "         -h             help\n")	;
	return 0;
}

int parse_args(int argc, char *argv[], opt *o)
{
	o->max_oh = 1000;
	o->ratio = 0.8;
	o->min_ovlp = 5000;
	o->out = NULL;
	o->fmt = 1;
	int c;
	while ((c = getopt(argc, argv, "r:m:M:o:O:h")) != -1) {
		switch (c) {
			case 'r':
				o->ratio = strtof(optarg, NULL);
				break;
			case 'M':
				o->max_oh = atoi(optarg);
				break;
			case 'm':
				o->min_ovlp = atoi(optarg);
				break;
			case 'O':
				if (!strcmp("FASTA", optarg))
					o->fmt = 0;
				break;
			case 'o':
				o->out = optarg;
			case 'h':
				help();
				return 1;
			default:
				fprintf(stderr,"[E::%s] undefined option %c\n", __func__, c);
				help();
				return 1;
		}
	}	
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] neither gfa nor paf file can be omitted!\n", __func__);
		help();
		return 1;
	} else {
		o->gfa_fn = argv[optind++];
		o->paf_fn = argv[optind];
	}
	return 0;
}

/*opt *init_opt() */
/*{*/
	/*opt *o = (opt *)calloc(sizeof())*/

/*}*/
/*int destroy_opt(opt *o)*/
/*{*/
	/*if (o->)*/


/*}*/

