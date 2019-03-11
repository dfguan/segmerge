/*
 * =====================================================================================
 *
 *       Filename:  eg.c
 *
 *    Description:  enrich GFA by using miniasm assembly
 *
 *        Version:  1.0
 *        Created:  02/05/2018 19:40:57
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
/*#include <omp.h>*/

#include "kvec.h"
#include "paf.h"
#include "sdict.h"
#include "opt.h"

#define eg_idx_qn(a, b) (((a).qns>>32) < ((b).qns>>32))
#define eg_idx_tn(a, b) (((a).tns >> 32) < (((b).tns >> 32)))
/*#define INT_MATCH 1*/
/*#define OVLP_MATCH 2*/
/*#define CON_MATCH 3*/
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))


typedef struct {
	uint64_t qns, tns;
	uint32_t qe, te;
	uint32_t ql, tl;
	uint32_t ml:30, rev:1, del:1;
	uint32_t bl:30, tail:1, con:1; //mapped at the end of reference, contained by another query
}eg_hit_t;

typedef struct {size_t n, m; eg_hit_t *a;} eg_hit_v;
typedef struct {size_t n; eg_hit_t *rht; uint64_t *idx;} eg_hits_t;

/**
 * @func	print_hit
 * @brief	print single hit
 *
 */


/*int print_hit(eg_hit_t *r, gfa_t *gf)*/
/*{*/
	/*size_t i = 0; eg_hit_t * rht = r; fprintf(stderr, "%s\t%u\t%u\t%u\t%llu\t%u\t%u\t%u\t%d\t%d\n", gf->seg[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].tns>>32 , rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);*/
	/*return 0;*/
/*}*/

//print a single hit
int print_hit(eg_hit_t *rht, sdict_t *tn)
{
	size_t i = 0;
	/*fprintf(stderr, "Order: QID\tQLEN\tQS\tQE\tTID\tTLEN\tTS\tTE\n");*/
	if (!rht[i].del) fprintf(stdout, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%d\t%d\n", tn->seq[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].rev?'-':'+', tn->seq[rht[i].tns>>32].name, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);
	return 0;
}
/**
 * @func   print_hits
 * @brief  print
 */

int print_hits(eg_hit_t *rht, size_t s, size_t e, sdict_t *tn)
{
	size_t i;
	/*fprintf(stderr, "Order: QID\tQLEN\tQS\tQE\tTID\tTLEN\tTS\tTE\n");*/
	for (i = s; i < e; ++i) { 
		/*if (!rht[i].del) fprintf(stdout, "%u: %s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%d\t%d\n", i, tn->seq[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].rev?'-':'+', tn->seq[rht[i].tns>>32].name, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);*/
		if (!rht[i].del) fprintf(stdout, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%d\t%d\n", tn->seq[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].rev?'-':'+', tn->seq[rht[i].tns>>32].name, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].bl, rht[i].ml);
	}
	return 0;
}


/**
 * @func   f_qn
 * @brief  return qn 
 *
 *
 */


uint32_t f_qn (eg_hit_t *r)
{
	return r->qns >> 32;
}
/**
 * @func   f_tn
 * @brief  return tn 
 */
uint32_t f_tn (eg_hit_t *r)
{
	return r->tns >> 32;
}	

/**
 * @func  index
 * @brief index alignments according to fun
 */

uint64_t *hit_index(eg_hit_t *rht, size_t n_rht, size_t n_idx, uint32_t (*f)(eg_hit_t *))
{
	uint64_t *idx = (uint64_t *)calloc(n_idx, sizeof(uint64_t)); //don't use malloc here cause some of the n_idx doesn't have value or pass **idx here, 
	if (!idx) {
		fprintf(stderr, "[E::%s] failed to allocate memory space, required %lu\n", __func__, n_idx * sizeof(uint64_t));
		exit(1);
	}
	
	/*fprintf(stderr, "%p\n",idx);*/
	size_t i, last;
	n_idx = 0;
	for (i=1, last = 0; i <= n_rht; ++i) {
		if (i == n_rht || f(rht+i) != f(rht + last)) {
			/*fprintf(stderr, "%u\t%u\n", i, last);*/
			idx[f(rht+last)] = (uint64_t)last << 32 | (i - last); //don't forget the left bracket here !!!
			/*fprintf(stderr, "%u\n", f(rht+i-1));	*/
			last = i;
			/*++*n_idx;*/
		}
	}
	/*for(i = 0; i <*n_idx; ++i)	fprintf(stderr, "this:%u\t%u\t%u\n", i, idx[i]>>32, (uint32_t)idx[i]);*/
	return idx;
}

/**
 * @func    cmp_qtn
 * @brief   compare query and target's names then query start position
 *
 */


int cmp_qtn (const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	uint32_t pn = p->qns >> 32;
	uint32_t qn = q->qns >> 32;
	if (pn == qn) {
		uint32_t ptn = p->tns >> 32;
		uint32_t qtn = q->tns >> 32;	
		if (ptn > qtn) return -1;
		else if (ptn == qtn) {
			if (p->qns > q->qns) return 1;
			else if (p->qns == q->qns) return 0;
			else return -1;
		} 
		else return 1;	
	} else if (pn > qn) return -1;
	else return 1;
}
/**
 * @func    cmp_q
 * @brief   compare query's 
 *
 */


int cmp_q (const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	if (p->qns == q->qns) {
		if (p->qe > q->qe) return -1;
		else if (p->qe == q->qe) return 0;
		else return 1;
	} else if (p->qns > q->qns) return 1;
	else return -1;
}

/**
 * @func    cmp_t
 * @brief   compare target's 
 *
 */
//tn   
int cmp_t(const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	if (p->tns == q->tns) {
		if (p->te > q->te) return -1;
		else if (p->te == q->te) return 0;
		else return 1;
	} else if (p->tns > q->tns) return 1;
	else return -1;
	/*fprintf(stderr, "%d\t",z);*/
}

/**
 * @func  mt_sort
 * @brief mulitple threads sort according to tn,ts,te using openmp
 */

int mt_sort(eg_hit_t *rht, uint64_t *idx, size_t n_ind, int (*cmp) (const void *r, const void *s))
{
	size_t i;
	/*#pragma parallel for number_threads(10) //10 threads*/
	for (i=0; i < n_ind; ++i) qsort(&rht[idx[i]>>32], (uint32_t)idx[i], sizeof(eg_hit_t), cmp);		
	return 0;	
}


/**
 * @func  set_cov 
 * @brief set covered alignments 
 *
 */

/*int set_cov(rg_hit_t *rht)*/
/*{*/
	/*return 0;*/
/*}*/

int merge_seg_core2(eg_hit_t *rht, size_t s, size_t e, uint32_t gs, sdict_t *tn)
{
    /*
     qry -----------------
	 ref
	 |   \  \
	 |   \ \
	 |   
	 |
	 |
	 */ 
	size_t i, j, k;
 	//score: ml backtrace: bl 
	for ( i = s; i < e; ++i) {
		rht[i].tail = 1;
		rht[i].bl = -1;
		uint32_t qml = rht[i].ml;
		rht[i].del = 0;
		/*fprintf(stdout, "(%u, %u)\n", i, __LINE__ );*/
		/*print_hit(rht+i,tn);*/
		for ( j = i - 1; j + 1 > s; --j) {
			if (rht[i].rev == rht[j].rev) {
				// j u < i v
				if (rht[i].rev) {
					/*fprintf(stdout, "(%u, %u, %u)\n", j, i, __LINE__ );*/
					/*print_hit(rht+j,tn);*/
					/*print_hit(rht+i, tn);*/
					if (rht[j].qns < rht[i].qns && rht[j].qe < rht[i].qe && rht[j].te > rht[i].te && rht[j].tns > rht[i].tns) {
						if ((uint32_t)rht[i].qns > rht[j].qe && (uint32_t)rht[i].qns - rht[j].qe > gs) continue;
					/*fprintf(stdout, "pass1 (%u, %u, %u)\n", j, i, __LINE__ );*/
						if ((uint32_t)rht[j].tns > rht[i].te && (uint32_t)rht[j].tns - rht[i].te > gs) continue;	
					/*fprintf(stdout, "pass (%u, %u, %u)\n", j, i, __LINE__ );*/
						//otherwise
					} else 
						continue;
				} else {
					if (rht[j].qns < rht[i].qns && rht[j].qe < rht[i].qe && rht[j].te < rht[i].te && rht[j].tns < rht[i].tns) {
						if ((uint32_t)rht[i].qns > rht[j].qe && (uint32_t)rht[i].qns - rht[j].qe > gs) continue;
						if ((uint32_t)rht[i].tns > rht[j].te && (uint32_t) rht[i].tns - rht[j].te > gs) continue;	
					} else 
						continue;	
				}		
				rht[j].tail = 0;
				rht[j].del = 1;
				if (rht[i].ml < rht[j].ml + qml) {
					rht[i].ml = rht[j].ml + qml;
					rht[i].bl = j;
				}	
			}		
		}
	}
	//update coordinate
	uint32_t max = 0;	
	uint32_t max_idx = -1;
	uint32_t sent = 0X3FFFFFFF;	
	for (i = s; i < e; ++i) {
		if (rht[i].tail) {
			//backtrace to find start 
			if (rht[i].ml > max) {
				max = rht[i].ml;
				max_idx = i;
			} 
			for (k = i; rht[k].bl != sent; k = rht[k].bl); 
			rht[i].qns = rht[k].qns;
			if (rht[i].rev) rht[i].te = rht[k].te;
			else rht[i].tns = rht[k].tns;	
		}
	}	
	return max_idx;
}


int merge_seg_core(eg_hit_t *rht, size_t s, size_t e, uint32_t max_gs, sdict_t *tn)
{
	size_t i, j;
	uint32_t rtn;
			/*print_hits(rht, s,e, tn);*/
	for ( i = j = s; i <= e; ++i) {
		if (i == e || (rht[j].tns >> 32) != (rht[i].tns >> 32)) {
			/*fprintf(stdout, "index %d\t%d\n", j, i);*/

			/*print_hits(rht, j,i, tn);*/
			rtn = merge_seg_core2(rht, j, i, max_gs, tn);
			/*if (~rtn) print_hit(rht+rtn, tn);*/
			j = i;	
		} 
	}
	return 0;
}

/**
 * @func    merge_segs
 * @brief   merge segment if they are aligned closed to each other 
 * @alg     construct alignment graph, and find paths.  
 */
int merge_segs(eg_hit_t *rht, uint64_t *idx, size_t n_idx, uint32_t max_gs, sdict_t *tn)
{
	size_t j;
	/*#pragma parallel for number_threads(4)*/
	for (j = 0; j < n_idx; ++j) {
		merge_seg_core(rht , idx[j] >> 32, (idx[j] >> 32) + (uint32_t)idx[j], max_gs, tn); 
	}
	return 0;
}
/**
 * @func   rm_neml
 * @brief  remove alignments without enough mapped length or internal match (not sure if appropriate)
 *
 */
#define MAX_LEFT 5000 //maybe use fraction
int flt_alns(eg_hit_t *rht, size_t n_rht, int min_ovlp, float int_fract)
{
	uint64_t i = 0;
	for (i = 0; i < n_rht; ++i) {
		//not deleted and not at the end of reference  
		if (!rht[i].del) {
				int32_t ovlp = rht[i].qe - (int32_t)rht[i].qns;
				if (ovlp < min_ovlp || ovlp < rht[i].ql * int_fract ) {
					if ((uint32_t)rht[i].tns > MAX_LEFT && rht[i].tl - rht[i].te > MAX_LEFT) 
						rht[i].tail = 0, rht[i].del = 1;
					else 
						rht[i].tail = 1;
				} else rht[i].tail = 0; 
		}
	}
	return 0;
}

/** 
 * @func	is_cont
 * @brief	if aln q is contained in aln p 
 * @notice	not helpful don't use
 */


int is_cont(eg_hit_t *p, eg_hit_t *q) 
{
	uint32_t max_ts = max((uint32_t) p->tns, (uint32_t)q->tns);
	uint32_t min_te = min(p->te, q->te);
	uint32_t p_qs, p_qe, q_qs, q_qe;
	int32_t p_ext5, p_ext3, q_ext5, q_ext3;
	p_qs = (uint32_t) p->qns + max_ts - (uint32_t) p->tns;
	p_qe = p->qe + min_te - p->te;
	q_qs = (uint32_t) q->qns + max_ts - (uint32_t) q->tns;
	q_qe = q->qe + min_te - q->te;
	if (p->rev) 
		p_ext5 = p->ql - p->qe, p_ext3 = p_qs;
	else
		p_ext3 = p->ql - p->qe, p_ext5 = p_qs;
	if (q->rev) 
		q_ext5 = q->ql - q->qe, q_ext3 = q_qs;
	else
		q_ext3 = q->ql - q->qe, q_ext5 = q_qs;
	return !(p_ext5-q_ext5) || !(p_ext3-q_ext3) || ((p_ext5 - q_ext5) ^ (p_ext3 - q_ext3)) >= 0; 
}


/**
 * @func   set_cont
 * @brief  mark contained queries
 */

int set_cont(eg_hit_t *rht, uint64_t *idx, size_t n_idx)
{
	size_t i = 0, last;
		
	for (i = 0; i < n_idx; ++i) {
		for (last = idx[i]>>32, i = last + 1; i < (idx[i]>>32) + (uint32_t)idx[i]; ++i) {
			/*if (is_cont(rht + i, rht + i-1))*/
			if (rht[i].te < rht[i-1].te)
				rht[i].con = 1;
			else {
				rht[i].con = 0;	
				last = i;
			}
		}
	}
	return 0;	
}

/**
 * @func   cleanup_alns 
 * @brief  clean up alns after 
 */

size_t cleanup_alns(eg_hit_t *rht, size_t n_rht)
{
	size_t i, j;
	size_t n_del = 0;
	for (i = 0, j = 0; i < n_rht; ++i) 
		if (!rht[i].del) 
			rht[j++] = rht[i], ++n_del;
	return n_rht - n_del;
}

/**
 * @func     update_cords
 * @brief    update coordinates for alns 
 * @notice   keep ml bl unchanged
 */
int update_cords(eg_hit_t *rht, size_t n_rht)
{
	size_t i;
	for (i = 0; i < n_rht; ++i) {
		if (rht[i].tail) 
			continue;
		int32_t ext5, ext3; 
		if (rht[i].rev) 
			ext3 = rht[i].ql - rht[i].qe, ext5 = (int32_t)rht[i].qns;
		else
			ext5 = rht[i].ql - rht[i].qe, ext3 = (int32_t)rht[i].qns;
		rht[i].qns = rht[i].qns >> 32 << 32;
		rht[i].qe = rht[i].ql;
		rht[i].tns = (rht[i].tns >> 32 << 32) | ((int32_t) rht[i].tns > ext5 ? (int32_t) rht[i].tns - ext5 : 0);
		rht[i].te =  rht[i].te + ext3; 	
	}
	return 0;
}


/**
 * @func    eg_hit_read
 * @brief   read alns from paf file *
 */
eg_hit_t *eg_hit_read(char *paf_fn, sdict_t* tn, size_t *n, uint32_t min_match)
{
	
	paf_file_t *fp;
	paf_rec_t r;
	eg_hit_v h = {0, 0, 0};
	fp = paf_open(paf_fn);
	if (!fp) {
		fprintf(stderr, "[E::%s] can not open PAF file %s\n",__func__, paf_fn);
		exit(1);	
	}	
	while (paf_read(fp, &r) >= 0) {
		/*int32_t n2id = gfa_name2id(gf, r.qn);*/
		/*if (n2id < 0) {*/
			/*fprintf(stderr, "[W::%s] query name %s not found in GFA file, won't be saved\n", __func__, r.qn);*/
			/*continue;*/
		/*}*/
		if (r.ml <= min_match || strcmp(r.qn, r.tn) >= 0) continue; //self alignment
		eg_hit_t *p;
		kv_pushp(eg_hit_t, h, &p);
		p->qns = (uint64_t)sd_put(tn, r.qn, r.ql)<<32 | r.qs;	
		p->ql = r.ql; p->qe = r.qe;
		p->tns = (uint64_t)sd_put(tn, r.tn, r.tl) << 32 | r.ts;
		p->tl = r.tl; p->te = r.te;
		p->rev = r.rev; p->ml = r.ml; p->bl = r.bl; 	
		p->con = 0; p->del = 0;	
	}
	paf_close(fp);
	*n = h.n;
	return h.a;
}	



eg_hits_t *eg_init()
{
	return (eg_hits_t *)calloc(1, sizeof(eg_hits_t));
}


int eg_destroy(eg_hits_t *r)
{
	if (r) {
		if (r->idx) free(r->idx);
		if (r->rht) free(r->rht);
		free(r);
	}
	return 0;

}


int main(int argc, char *argv[])
{
	opt opts;
	if (parse_args(argc, argv, &opts))    return 1;		
	
	/*fprintf(stderr,"[M::%s] parsing gfa...\n", __func__);*/
	//read gfa from gfa file
	/*gfa_t *gf = gfa_read(opts.gfa_fn);*/
	/*gfa_print(gf, stdout, 1);*/
	//read alns from paf file
	sdict_t *rn = sd_init();
	eg_hits_t *rhts = eg_init();

	/*fprintf(stderr,"[M::%s] parsing paf...\n", __func__);*/
	rhts->rht = eg_hit_read(opts.paf_fn, rn, &rhts->n, opts.min_bl);		
	/*print_hits(rhts->rht, 0, rhts->n, rn);*/
	/*fprintf(stderr, "%llu\n", rhts->n);	*/
	/*fprintf(stderr,"[M::%s] indexing query...\n", __func__);*/
	size_t n_ind = rn->n_seq;
	rhts->idx = hit_index(rhts->rht, rhts->n, n_ind, f_qn);
	/*fprintf(stderr, "%p\n", rhts->idx);*/
	/*size_t i = 0;*/
	/*for ( i = 0; i < n_ind; ++i) { fprintf(stderr, "%u : %u\n", rhts->idx[i]>>32, (uint32_t)rhts->idx[i]); }	*/
	/*fprintf(stderr,"[M::%s] sorting target...\n", __func__);*/
	//sort according to tn,ts,te
	mt_sort(rhts->rht, rhts->idx, n_ind, cmp_qtn);
	//merge near alns
	/*fprintf(stderr,"[M::%s] merging alignments...\n", __func__);*/
	/*print_hits(rhts->rht, 0, rhts->n, rn);*/
	merge_segs(rhts->rht, rhts->idx, n_ind, opts.max_gs, rn);
	print_hits(rhts->rht, 0, rhts->n, rn);
	/*print_hits(rhts->rht, rhts->n, rn);*/
	//rm not engouh mapped length and ! OVLP and internal matching, update coordination
	/*fprintf(stderr,"[M::%s] filtering alignments...\n", __func__);*/
	/*flt_alns(rhts->rht, rhts->n, opts.min_ovlp, opts.ratio);*/
	//print kept alns probably only one or two alignment was kept 
	//clean up	
		
	/*rhts->n = cleanup_alns(rhts->rht, rhts->n);*/
	//sort according to tn,ts,te
	/*fprintf(stderr, "%lu\t%lu\n", rhts->n, sizeof(eg_hit_t));*/
	/*qsort(rhts->rht, rhts->n, sizeof(eg_hit_t), cmp_t);	*/
	//index target here
	/*if (rhts->idx) free(rhts->idx);*/
	/*rhts->idx = hit_index(rhts->rht, rhts->n, &n_ind, f_tn);*/
	/*rhts->idx = hit_index(rhts->rht, rhts->n, &n_ind, f_qn);*/
	//set contained
	/*set_cont(rhts->rht, rhts->idx, n_ind);*/
	//adjust coordinate
	/*print_hits(rhts->rht, rhts->n, gf, tn);*/
	/*update_cords(rhts->rht, rhts->n);*/
	//update gfa
	/*fprintf(stderr,"[M::%s] updating gfa...\n", __func__);*/
	/*print_hits(rhts->rht, rhts->n, rn);*/
	/*update_gfa(rhts->rht, rhts->idx, n_ind, gf);*/

	/*gfa_t *ug = gfa_ug_gen(gf);	*/
	/*FILE *fp_out = opts.out ? fopen(opts.out, "w") : stdout;*/
	/*if (opts.fmt) */
		/*gfa_print(ug, fp_out, 1);*/
	/*else */
		/*fa_print(ug, fp_out);*/
	
	/*if (opts.out) fclose(fp_out);*/
	if (rhts) eg_destroy(rhts);	
	sd_destroy(rn);	
	/*if (gf) gfa_destroy(gf);*/
	/*if (ug) gfa_destroy(ug);*/
	return 0;	
}




