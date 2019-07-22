#include "convhull_3d.h"

/************
 * INTERNAL:
 ***********/
 /*
	 Redesigned, integrated, and improved by Rinat Mukhometzianov, 2019
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define CV_STRNCPY(a,b,c) strncpy_s(a,c+1,b,c);
#define CV_STRCAT(a,b) strcat_s(a,sizeof(b),b);
#else
#define CV_STRNCPY(a,b,c) strncpy(a,b,c);
#define CV_STRCAT(a,b) strcat(a,b);
#endif

#ifdef CONVHULL_3D_USE_FLOAT_PRECISION
#define CH_FLT_MIN FLT_MIN
#define CH_FLT_MAX FLT_MAX
#define CH_NOISE_VAL 0.00001f
#else
#define CH_FLT_MIN DBL_MIN
#define CH_FLT_MAX DBL_MAX
#define CH_NOISE_VAL 0.0000001
#endif

#ifndef MIN
#define MIN(a,b) (( (a) < (b) ) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )
#endif
#define CH_MAX_NUM_FACES 50000

 /* structs for qsort */
typedef struct float_w_idx {
	CH_FLOAT val;
	unsigned long idx;
}float_w_idx;

typedef struct int_w_idx {
	unsigned long val;
	unsigned long idx;
}int_w_idx;

/* internal functions prototypes: */
static int cmp_asc_float(const void*, const void*);
static int cmp_desc_float(const void*, const void*);
static int cmp_asc_int(const void*, const void*);
static int cmp_desc_int(const void*, const void*);
static void sort_float(CH_FLOAT*, CH_FLOAT*, unsigned long*, unsigned long, unsigned long);
static void sort_int(unsigned long*, unsigned long*, unsigned long*, unsigned long, unsigned long);
static ch_vec3 cross(ch_vec3*, ch_vec3*);
static CH_FLOAT det_4x4(CH_FLOAT*);
static void plane_3d(CH_FLOAT*, CH_FLOAT*, CH_FLOAT*);
static void ismember(unsigned long*, unsigned long*, unsigned long*, unsigned long, unsigned long);

/* internal functions definitions: */
static int cmp_asc_float(const void* a, const void* b) {
	struct float_w_idx* a1 = (struct float_w_idx*)a;
	struct float_w_idx* a2 = (struct float_w_idx*)b;
	if ((*a1).val < (*a2).val)return -1;
	else if ((*a1).val > (*a2).val)return 1;
	else return 0;
}

static int cmp_desc_float(const void* a, const void* b) {
	struct float_w_idx* a1 = (struct float_w_idx*)a;
	struct float_w_idx* a2 = (struct float_w_idx*)b;
	if ((*a1).val > (*a2).val)return -1;
	else if ((*a1).val < (*a2).val)return 1;
	else return 0;
}

static int cmp_asc_int(const void* a, const void* b) {
	struct int_w_idx* a1 = (struct int_w_idx*)a;
	struct int_w_idx* a2 = (struct int_w_idx*)b;
	if ((*a1).val < (*a2).val)return -1;
	else if ((*a1).val > (*a2).val)return 1;
	else return 0;
}

static int cmp_desc_int(const void* a, const void* b) {
	struct int_w_idx* a1 = (struct int_w_idx*)a;
	struct int_w_idx* a2 = (struct int_w_idx*)b;
	if ((*a1).val > (*a2).val)return -1;
	else if ((*a1).val < (*a2).val)return 1;
	else return 0;
}

static void sort_float
(
	CH_FLOAT* in_vec,  /* vector[len] to be sorted */
	CH_FLOAT* out_vec, /* if NULL, then in_vec is sorted "in-place" */
	unsigned long* new_idices,   /* set to NULL if you don't need them */
	unsigned long len,           /* number of elements in vectors, must be consistent with the input data */
	unsigned long descendFLAG    /* !1:ascending, 1:descending */
)
{
	unsigned long i;
	struct float_w_idx* data;

	data = (float_w_idx*)malloc(len * sizeof(float_w_idx));
	if (data) {
		for (i = 0; i < len; i++) {
			data[i].val = in_vec[i];
			data[i].idx = i;
		}
		if (descendFLAG)
			qsort(data, len, sizeof(data[0]), cmp_desc_float);
		else
			qsort(data, len, sizeof(data[0]), cmp_asc_float);
		for (i = 0; i < len; i++) {
			if (out_vec != NULL)
				out_vec[i] = data[i].val;
			else
				in_vec[i] = data[i].val; /* overwrite input vector */
			if (new_idices != NULL)
				new_idices[i] = data[i].idx;
		}
	}
	else {
		printf("Can't use sort_float. Not enough memory\n");
		exit(0);
	}
	free(data);
}

static void sort_int
(
	unsigned long* in_vec,     /* vector[len] to be sorted */
	unsigned long* out_vec,    /* if NULL, then in_vec is sorted "in-place" */
	unsigned long* new_idices, /* set to NULL if you don't need them */
	unsigned long len,         /* number of elements in vectors, must be consistent with the input data */
	unsigned long descendFLAG  /* !1:ascending, 1:descending */
)
{
	unsigned long i;
	struct int_w_idx* data;

	data = (int_w_idx*)malloc(len * sizeof(int_w_idx));
	if (data) {
		for (i = 0; i < len; i++) {
			data[i].val = in_vec[i];
			data[i].idx = i;
		}
		if (descendFLAG)
			qsort(data, len, sizeof(data[0]), cmp_desc_int);
		else
			qsort(data, len, sizeof(data[0]), cmp_asc_int);
		for (i = 0; i < len; i++) {
			if (out_vec != NULL)
				out_vec[i] = data[i].val;
			else
				in_vec[i] = data[i].val; /* overwrite input vector */
			if (new_idices != NULL)
				new_idices[i] = data[i].idx;
		}
	}
	else {
		printf("Can't use sort_float. Not enough memory\n");
		exit(0);
	}
	free(data);
}

static ch_vec3 cross(ch_vec3* v1, ch_vec3* v2)
{
	ch_vec3 cross;
	cross.x = v1->y * v2->z - v1->z * v2->y;
	cross.y = v1->z * v2->x - v1->x * v2->z;
	cross.z = v1->x * v2->y - v1->y * v2->x;
	return cross;
}

/* calculates the determinent of a 4x4 matrix */
static CH_FLOAT det_4x4(CH_FLOAT* m) {
	return
		m[3] * m[6] * m[9] * m[12] - m[2] * m[7] * m[9] * m[12] -
		m[3] * m[5] * m[10] * m[12] + m[1] * m[7] * m[10] * m[12] +
		m[2] * m[5] * m[11] * m[12] - m[1] * m[6] * m[11] * m[12] -
		m[3] * m[6] * m[8] * m[13] + m[2] * m[7] * m[8] * m[13] +
		m[3] * m[4] * m[10] * m[13] - m[0] * m[7] * m[10] * m[13] -
		m[2] * m[4] * m[11] * m[13] + m[0] * m[6] * m[11] * m[13] +
		m[3] * m[5] * m[8] * m[14] - m[1] * m[7] * m[8] * m[14] -
		m[3] * m[4] * m[9] * m[14] + m[0] * m[7] * m[9] * m[14] +
		m[1] * m[4] * m[11] * m[14] - m[0] * m[5] * m[11] * m[14] -
		m[2] * m[5] * m[8] * m[15] + m[1] * m[6] * m[8] * m[15] +
		m[2] * m[4] * m[9] * m[15] - m[0] * m[6] * m[9] * m[15] -
		m[1] * m[4] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
}

/* Calculates the coefficients of the equation of a PLANE in 3D.
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 */
static void plane_3d
(
	CH_FLOAT* p,
	CH_FLOAT* c,
	CH_FLOAT* d
)
{
	unsigned long i, j, k, l;
	unsigned long r[3];
	CH_FLOAT sign, det, norm_c;
	CH_FLOAT pdiff[2][3], pdiff_s[2][2];

	for (i = 0; i < 2; i++)
		for (j = 0; j < 3; j++)
			pdiff[i][j] = p[(i + 1) * 3 + j] - p[i * 3 + j];
	memset(c, 0, 3 * sizeof(CH_FLOAT));
	sign = 1.0;
	for (i = 0; i < 3; i++)
		r[i] = i;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 2; j++) {
			for (k = 0, l = 0; k < 3; k++) {
				if (r[k] != i) {
					pdiff_s[j][l] = pdiff[j][k];
					l++;
				}
			}
		}
		det = pdiff_s[0][0] * pdiff_s[1][1] - pdiff_s[1][0] * pdiff_s[0][1];
		c[i] = sign * det;
		sign *= -1.0;
	}
	norm_c = (CH_FLOAT)0.0;
	for (i = 0; i < 3; i++)
		norm_c += (pow(c[i], 2.0));
	norm_c = sqrt(norm_c);
	for (i = 0; i < 3; i++)
		c[i] /= norm_c;
	(*d) = (CH_FLOAT)0.0;
	for (i = 0; i < 3; i++)
		(*d) += -p[i] * c[i];
}

static void ismember
(
	unsigned long* pLeft,          /* left vector; nLeftElements x 1 */
	unsigned long* pRight,         /* right vector; nRightElements x 1 */
	unsigned long* pOut,           /* 0, unless pRight elements are present in pLeft then 1; nLeftElements x 1 */
	unsigned long nLeftElements,   /* number of elements in pLeft */
	unsigned long nRightElements   /* number of elements in pRight */
)
{
	unsigned long i, j;
	memset(pOut, 0, nLeftElements * sizeof(unsigned long));
	for (i = 0; i < nLeftElements; i++)
		for (j = 0; j < nRightElements; j++)
			if (pLeft[i] == pRight[j])
				pOut[i] = 1;
}

/* A C version of the 3D quickhull matlab implementation from here:
 * https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox?focused=3851550&tab=example
 * (*out_faces) is returned as NULL, if triangulation fails *
 * Original Copyright (c) 2014, George Papazafeiropoulos
 * Distributed under the BSD (2-clause) license
 * Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford Barber, David P. Dobkin
 *             and Hannu Huhdanpaa, Geometry Center Technical Report GCG53, July 30, 1993"
 */
void convhull_3d_build
(
	ch_vertex* const in_vertices,
	const unsigned long nVert,
	unsigned long** out_faces,
	unsigned long* nOut_faces
)
{
	unsigned long i, j, k, l, h;
	unsigned long nFaces, p, d;
	unsigned long* aVec, * faces;
	CH_FLOAT dfi, v, max_p, min_p;
	CH_FLOAT* points, * cf, * cfi, * df, * p_s, * span;

	if (nVert < 3 || in_vertices == NULL) {
		(*out_faces) = NULL;
		(*nOut_faces) = 0;
		return;
	}

	/* 3 dimensions. The code should theoretically work for >=2 dimensions, but "plane_3d" and "det_4x4" are hardcoded for 3,
	 * so would need to be rewritten */
	d = 3;
	span = (CH_FLOAT*)malloc(d * sizeof(CH_FLOAT));
	for (j = 0; j < d; j++) {
		max_p = 2.23e-13; min_p = 2.23e+13;
		for (i = 0; i < nVert; i++) {
			max_p = MAX(max_p, in_vertices[i].v[j]);
			min_p = MIN(min_p, in_vertices[i].v[j]);
		}
		span[j] = max_p - min_p;
	}
	points = (CH_FLOAT*)malloc(nVert * (d + 1) * sizeof(CH_FLOAT));
	if (points) {
		for (i = 0; i < nVert; i++) {
			for (j = 0; j < d; j++)
				points[i * (d + 1) + j] = in_vertices[i].v[j] + CH_NOISE_VAL * rand() / (float)RAND_MAX; /* noise mitigates duplicates */
			points[i * (d + 1) + d] = 1.0f; /* add a last column of ones. Used only for determinant calculation */
		}
	}
	else {
		printf("Can't arrange enough memory for points variable. Can't build convex hull");
		exit(0);
	}

	/* The initial convex hull is a simplex with (d+1) facets, where d is the number of dimensions */
	nFaces = (d + 1);
	faces = (unsigned long*)calloc(nFaces * d, sizeof(unsigned long));
	aVec = (unsigned long*)malloc(nFaces * sizeof(unsigned long));
	if (aVec) {
		for (i = 0; i < nFaces; i++)
			aVec[i] = i;
	}
	else {
		printf("Can't initialize aVec");
		exit(0);
	}

	/* Each column of cf contains the coefficients of a plane */
	cf = (CH_FLOAT*)malloc(nFaces * d * sizeof(CH_FLOAT));
	cfi = (CH_FLOAT*)malloc(d * sizeof(CH_FLOAT));
	df = (CH_FLOAT*)malloc(nFaces * sizeof(CH_FLOAT));
	p_s = (CH_FLOAT*)malloc(d * d * sizeof(CH_FLOAT));
	for (i = 0; i < nFaces; i++) {
		/* Set the indices of the points defining the face  */
		if (faces) {
			for (j = 0, k = 0; j < (d + 1); j++) {
				if (aVec[j] != i) {
					faces[i * d + k] = aVec[j];
					k++;
				}
			}
		}
		else {
			printf("Can't arrange enough memory for faces variable. Can't build convex hull");
			exit(0);
		}

		/* Calculate and store the plane coefficients of the face */
		if (p_s) {
			for (j = 0; j < d; j++)
				for (k = 0; k < d; k++)
					p_s[j * d + k] = points[(faces[i * d + j]) * (d + 1) + k];
		}
		else {
			printf("Can't arrange enough memory for p_s variable. Can't build convex hull");
			exit(0);
		}

		/* Calculate and store the plane coefficients of the face */
		plane_3d(p_s, cfi, &dfi);
		if (cf && cfi) {
			for (j = 0; j < d; j++)
				cf[i * d + j] = cfi[j];
		}
		else {
			printf("Can't arrange enough memory for cf or cfi variable. Can't build convex hull");
			exit(0);
		}
		if (df && dfi) {
			df[i] = dfi;
		}
		else {
			printf("Can't arrange enough memory for df or dfi variable. Can't build convex hull");
			exit(0);
		}
	}
	CH_FLOAT* A;
	unsigned long* bVec, * fVec, * asfVec, * face_tmp;

	/* Check to make sure that faces are correctly oriented */
	bVec = (unsigned long*)malloc(4 * sizeof(unsigned long));
	if (bVec) {
		for (i = 0; i < d + 1; i++)
			bVec[i] = i;
	}
	else {
		printf("Can't arrange enough memory for bVec variable. Can't build convex hull");
		exit(0);
	}

	/* A contains the coordinates of the points forming a simplex */
	A = (CH_FLOAT*)calloc((d + 1) * (d + 1), sizeof(CH_FLOAT));
	face_tmp = (unsigned long*)malloc((d + 1) * sizeof(unsigned long));
	fVec = (unsigned long*)malloc((d + 1) * sizeof(unsigned long));
	asfVec = (unsigned long*)malloc((d + 1) * sizeof(unsigned long));
	if (fVec && asfVec) {
		for (k = 0; k < (d + 1); k++) {
			/* Get the point that is not on the current face (point p) */

			for (i = 0; i < d; i++)
				fVec[i] = faces[k * d + i];
			sort_int(fVec, NULL, NULL, d, 0); /* sort accending */

			p = k;
			for (i = 0; i < d; i++)
				for (j = 0; j < (d + 1); j++)
					A[i * (d + 1) + j] = points[(faces[k * d + i]) * (d + 1) + j];
			for (; i < (d + 1); i++)
				for (j = 0; j < (d + 1); j++)
					A[i * (d + 1) + j] = points[p * (d + 1) + j];

			/* det(A) determines the orientation of the face */
			v = det_4x4(A);

			/* Orient so that each point on the original simplex can't see the opposite face */
			if (v < 0) {
				/* Reverse the order of the last two vertices to change the volume */
				for (j = 0; j < d; j++)
					face_tmp[j] = faces[k * d + j];
				for (j = 0, l = d - 2; j < d - 1; j++, l++)
					faces[k * d + l] = face_tmp[d - j - 1];

				/* Modify the plane coefficients of the properly oriented faces */
				for (j = 0; j < d; j++)
					cf[k * d + j] = -cf[k * d + j];
				df[k] = -df[k];
				for (i = 0; i < d; i++)
					for (j = 0; j < (d + 1); j++)
						A[i * (d + 1) + j] = points[(faces[k * d + i]) * (d + 1) + j];
				for (; i < (d + 1); i++)
					for (j = 0; j < (d + 1); j++)
						A[i * (d + 1) + j] = points[p * (d + 1) + j];
			}
		}
	}
	else {
		printf("Can't arrange enough memory for fVec or asfVec variable. Can't build convex hull");
		exit(0);
	}
	/* Coordinates of the center of the point set */
	CH_FLOAT* meanp, * absdist, * reldist, * desReldist;
	meanp = (CH_FLOAT*)calloc(d, sizeof(CH_FLOAT));
	for (i = d + 1; i < nVert; i++)
		for (j = 0; j < d; j++)
			meanp[j] += points[i * (d + 1) + j];
	for (j = 0; j < d; j++)
		meanp[j] = meanp[j] / (CH_FLOAT)(nVert - d - 1);

	/* Absolute distance of points from the center */
	absdist = (CH_FLOAT*)malloc((nVert - d - 1) * d * sizeof(CH_FLOAT));
	for (i = d + 1, k = 0; i < nVert; i++, k++)
		for (j = 0; j < d; j++)
			absdist[k * d + j] = (points[i * (d + 1) + j] - meanp[j]) / span[j];

	/* Relative distance of points from the center */
	reldist = (CH_FLOAT*)calloc((nVert - d - 1), sizeof(CH_FLOAT));
	desReldist = (CH_FLOAT*)malloc((nVert - d - 1) * sizeof(CH_FLOAT));
	for (i = 0; i < (nVert - d - 1); i++)
		for (j = 0; j < d; j++)
			reldist[i] += pow(absdist[i * d + j], 2.0);

	/* Sort from maximum to minimum relative distance */
	unsigned long num_pleft, cnt;
	unsigned long* ind, * pleft;
	ind = (unsigned long*)malloc((nVert - d - 1) * sizeof(unsigned long));
	pleft = (unsigned long*)malloc((nVert - d - 1) * sizeof(unsigned long));
	sort_float(reldist, desReldist, ind, (nVert - d - 1), 1);

	/* Initialize the vector of points left. The points with the larger relative
	 distance from the center are scanned first. */
	num_pleft = (nVert - d - 1);
	for (i = 0; i < num_pleft; i++)
		pleft[i] = ind[i] + d + 1;

	/* Loop over all remaining points that are not deleted. Deletion of points
	 occurs every #iter2del# iterations of this while loop */
	memset(A, 0, (d + 1) * (d + 1) * sizeof(CH_FLOAT));

	/* cnt is equal to the points having been selected without deletion of
	 nonvisible points (i.e. points inside the current convex hull) */
	cnt = 0;

	/* The main loop for the quickhull algorithm */
	CH_FLOAT detA;
	CH_FLOAT* points_cf, * points_s, * tmp_ch_float;
	unsigned long* visible_ind, * visible, * nonvisible_faces, * f0, * face_s,
		* u, * gVec, * horizon, * hVec, * pp, * hVec_mem_face, * tmp_long;

	unsigned long num_visible_ind, num_nonvisible_faces, n_newfaces, count, vis;
	unsigned long f0_sum, u_len, start, num_p, index, horizon_size1;
	unsigned long FUCKED;
	FUCKED = 0;
	u = horizon = NULL;
	nFaces = d + 1;
	visible_ind = (unsigned long*)malloc(nFaces * sizeof(unsigned long));
	points_cf = (CH_FLOAT*)malloc(nFaces * sizeof(CH_FLOAT));
	points_s = (CH_FLOAT*)malloc(d * sizeof(CH_FLOAT));
	face_s = (unsigned long*)malloc(d * sizeof(unsigned long));
	gVec = (unsigned long*)malloc(d * sizeof(unsigned long));
	while ((num_pleft > 0)) {
		/* i is the first point of the points left */
		i = pleft[0];

		/* Delete the point selected */
		for (j = 0; j < num_pleft - 1; j++)
			pleft[j] = pleft[j + 1];
		num_pleft--;
		if (num_pleft == 0)
			free(pleft);
		else {
			tmp_long = (unsigned long*)realloc(pleft, num_pleft * sizeof(unsigned long));
			if (tmp_long != NULL)
			{
				pleft = tmp_long;
			}
		}

		/* Update point selection counter */
		cnt++;

		/* find visible faces */
		for (j = 0; j < d; j++)
			points_s[j] = points[i * (d + 1) + j];

		tmp_ch_float = (CH_FLOAT*)realloc(points_cf, nFaces * sizeof(CH_FLOAT));
		if (tmp_ch_float != NULL)
		{
			points_cf = tmp_ch_float;
		}

		tmp_long = (unsigned long*)realloc(visible_ind, nFaces * sizeof(unsigned long));
		if (tmp_long != NULL)
		{
			visible_ind = tmp_long;
		}
#ifdef CONVHULL_3D_USE_CBLAS
#ifdef CONVHULL_3D_USE_FLOAT_PRECISION
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, nFaces, d, 1.0f,
			points_s, d,
			cf, d, 0.0f,
			points_cf, nFaces);
#else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, nFaces, d, 1.0,
			points_s, d,
			cf, d, 0.0,
			points_cf, nFaces);
#endif
#else
		for (j = 0; j < nFaces; j++) {
			points_cf[j] = 0;
			for (k = 0; k < d; k++)
				points_cf[j] += points_s[k] * cf[j * d + k];
		}
#endif
		num_visible_ind = 0;
		for (j = 0; j < nFaces; j++) {
			if (points_cf[j] + df[j] > 0.0) {
				num_visible_ind++; /* will sum to 0 if none are visible */
				visible_ind[j] = 1;
			}
			else
				visible_ind[j] = 0;
		}
		num_nonvisible_faces = nFaces - num_visible_ind;

		/* proceed if there are any visible faces */
		if (num_visible_ind != 0) {
			/* Find visible face indices */
			visible = (unsigned long*)malloc(num_visible_ind * sizeof(unsigned long));
			for (j = 0, k = 0; j < nFaces; j++) {
				if (visible_ind[j] == 1) {
					visible[k] = j;
					k++;
				}
			}

			/* Find nonvisible faces */

			nonvisible_faces = (unsigned long*)malloc(num_nonvisible_faces * d * sizeof(unsigned long));
			f0 = (unsigned long*)malloc(num_nonvisible_faces * d * sizeof(unsigned long));
			for (j = 0, k = 0; j < nFaces; j++) {
				if (visible_ind[j] == 0) {
					for (l = 0; l < d; l++)
						nonvisible_faces[k * d + l] = faces[j * d + l];
					k++;
				}
			}

			/* Create horizon (count is the number of the edges of the horizon) */
			count = 0;
			for (j = 0; j < num_visible_ind; j++) {
				/* visible face */
				vis = visible[j];
				for (k = 0; k < d; k++)
					face_s[k] = faces[vis * d + k];
				sort_int(face_s, NULL, NULL, d, 0);
				ismember(nonvisible_faces, face_s, f0, num_nonvisible_faces * d, d);
				u_len = 0;

				/* u are the nonvisible faces connected to the face v, if any */
				for (k = 0; k < num_nonvisible_faces; k++) {
					f0_sum = 0;
					for (l = 0; l < d; l++)
						f0_sum += f0[k * d + l];
					if (f0_sum == d - 1) {
						u_len++;
						if (u_len == 1)
							u = (unsigned long*)malloc(u_len * sizeof(unsigned long));
						else {
							tmp_long = (unsigned long*)realloc(u, u_len * sizeof(unsigned long));
							if (tmp_long != NULL)
							{
								u = tmp_long;
							}
						}
						u[u_len - 1] = k;
					}
				}
				for (k = 0; k < u_len; k++) {
					/* The boundary between the visible face v and the k(th) nonvisible face connected to the face v forms part of the horizon */
					count++;
					if (count == 1)
						horizon = (unsigned long*)malloc(count * (d - 1) * sizeof(unsigned long));
					else {
						tmp_long = (unsigned long*)realloc(horizon, count * (d - 1) * sizeof(unsigned long));
						if (tmp_long != NULL)
						{
							horizon = tmp_long;
						}
					}
					for (l = 0; l < d; l++)
						gVec[l] = nonvisible_faces[u[k] * d + l];
					for (l = 0, h = 0; l < d; l++) {
						if (f0[u[k] * d + l]) {
							horizon[(count - 1) * (d - 1) + h] = gVec[l];
							h++;
						}
					}
				}
				if (u_len != 0)
					free(u);
			}
			horizon_size1 = count;
			for (j = 0, l = 0; j < nFaces; j++) {
				if (!visible_ind[j]) {
					/* Delete visible faces */
					for (k = 0; k < d; k++)
						faces[l * d + k] = faces[j * d + k];

					/* Delete the corresponding plane coefficients of the faces */
					for (k = 0; k < d; k++)
						cf[l * d + k] = cf[j * d + k];
					df[l] = df[j];
					l++;
				}
			}

			/* Update the number of faces */
			nFaces = nFaces - num_visible_ind;
			tmp_long = (unsigned long*)realloc(faces, nFaces * d * sizeof(unsigned long));
			if (tmp_long != NULL)
			{
				faces = tmp_long;
			}
			tmp_ch_float = (CH_FLOAT*)realloc(cf, nFaces * d * sizeof(CH_FLOAT));
			if (tmp_ch_float != NULL)
			{
				cf = tmp_ch_float;
			}
			tmp_ch_float = (CH_FLOAT*)realloc(df, nFaces * sizeof(CH_FLOAT));
			if (tmp_ch_float != NULL)
			{
				df = tmp_ch_float;
			}

			/* start is the first row of the new faces */
			start = nFaces;

			/* Add faces connecting horizon to the new point */
			n_newfaces = horizon_size1;
			for (j = 0; j < n_newfaces; j++) {
				nFaces++;
				tmp_long = (unsigned long*)realloc(faces, nFaces * d * sizeof(unsigned long));
				if (tmp_long != NULL)
				{
					faces = tmp_long;
				}
				tmp_ch_float = (CH_FLOAT*)realloc(cf, nFaces * d * sizeof(CH_FLOAT));
				if (tmp_ch_float != NULL)
				{
					cf = tmp_ch_float;
				}
				tmp_ch_float = (CH_FLOAT*)realloc(df, nFaces * sizeof(CH_FLOAT));
				if (tmp_ch_float != NULL)
				{
					df = tmp_ch_float;
				}
				for (k = 0; k < d - 1; k++)
					faces[(nFaces - 1) * d + k] = horizon[j * (d - 1) + k];
				faces[(nFaces - 1) * d + (d - 1)] = i;

				/* Calculate and store appropriately the plane coefficients of the faces */
				for (k = 0; k < d; k++)
					for (l = 0; l < d; l++)
						p_s[k * d + l] = points[(faces[(nFaces - 1) * d + k]) * (d + 1) + l];
				plane_3d(p_s, cfi, &dfi);
				for (k = 0; k < d; k++)
					cf[(nFaces - 1) * d + k] = cfi[k];
				df[(nFaces - 1)] = dfi;
				if (nFaces > CH_MAX_NUM_FACES) {
					FUCKED = 1;
					nFaces = 0;
					break;
				}
			}

			/* Orient each new face properly */
			hVec = (unsigned long*)malloc(nFaces * sizeof(unsigned long));
			hVec_mem_face = (unsigned long*)malloc(nFaces * sizeof(unsigned long));
			for (j = 0; j < nFaces; j++)
				hVec[j] = j;
			for (k = start; k < nFaces; k++) {
				for (j = 0; j < d; j++)
					face_s[j] = faces[k * d + j];
				sort_int(face_s, NULL, NULL, d, 0);
				ismember(hVec, face_s, hVec_mem_face, nFaces, d);
				num_p = 0;
				for (j = 0; j < nFaces; j++)
					if (!hVec_mem_face[j])
						num_p++;
				pp = (unsigned long*)malloc(num_p * sizeof(unsigned long));
				for (j = 0, l = 0; j < nFaces; j++) {
					if (!hVec_mem_face[j]) {
						pp[l] = hVec[j];
						l++;
					}
				}
				index = 0;
				detA = 0.0;

				/* While new point is coplanar, choose another point */
				while (detA == 0.0) {
					for (j = 0; j < d; j++)
						for (l = 0; l < d + 1; l++)
							A[j * (d + 1) + l] = points[(faces[k * d + j]) * (d + 1) + l];
					for (; j < d + 1; j++)
						for (l = 0; l < d + 1; l++)
							A[j * (d + 1) + l] = points[pp[index] * (d + 1) + l];
					index++;
					detA = det_4x4(A);
				}

				/* Orient faces so that each point on the original simplex can't see the opposite face */
				if (detA < 0.0) {
					/* If orientation is improper, reverse the order to change the volume sign */
					for (j = 0; j < d; j++)
						face_tmp[j] = faces[k * d + j];
					for (j = 0, l = d - 2; j < d - 1; j++, l++)
						faces[k * d + l] = face_tmp[d - j - 1];

					/* Modify the plane coefficients of the properly oriented faces */
					for (j = 0; j < d; j++)
						cf[k * d + j] = -cf[k * d + j];
					df[k] = -df[k];
					for (l = 0; l < d; l++)
						for (j = 0; j < d + 1; j++)
							A[l * (d + 1) + j] = points[(faces[k * d + l]) * (d + 1) + j];
					for (; l < d + 1; l++)
						for (j = 0; j < d + 1; j++)
							A[l * (d + 1) + j] = points[pp[index] * (d + 1) + j];
				}
				free(pp);
			}
			free(horizon);
			free(f0);
			free(nonvisible_faces);
			free(visible);
			free(hVec);
			free(hVec_mem_face);
		}
		if (FUCKED) {
			break;
		}
	}

	/* output */
	if (FUCKED) {
		(*out_faces) = NULL;
		(*nOut_faces) = 0;
	}
	else {
		(*out_faces) = (unsigned long*)malloc(nFaces * d * sizeof(unsigned long));
		memcpy((*out_faces), faces, nFaces * d * sizeof(unsigned long));
		(*nOut_faces) = nFaces;
	}

	/* clean-up */
	free(visible_ind);
	free(points_cf);
	free(points_s);
	free(face_s);
	free(gVec);
	free(meanp);
	free(absdist);
	free(reldist);
	free(desReldist);
	free(ind);
	free(span);
	free(points);
	free(faces);
	free(aVec);
	free(cf);
	free(cfi);
	free(df);
	free(p_s);
	free(face_tmp);
	free(fVec);
	free(asfVec);
	free(bVec);
	free(A);
}

void convhull_3d_export_obj
(
	ch_vertex* const vertices,
	const unsigned long nVert,
	unsigned long* const faces,
	const unsigned long nFaces,
	const unsigned long keepOnlyUsedVerticesFLAG,
	char* const obj_filename
)
{
	unsigned long i, j;
	char path[256] = "\0";
	CV_STRNCPY(path, obj_filename, strlen(obj_filename));
	path[strlen(obj_filename) - 1] = 0;
	FILE* obj_file;
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
	CV_STRCAT(path, ".obj");
	fopen_s(&obj_file, path, "wt");
#else
	obj_file = fopen(strcat(path, ".obj"), "wt");
#endif
	fprintf(obj_file, "o\n");
	CH_FLOAT scale;
	ch_vec3 v1, v2, normal;

	/* export vertices */
	if (keepOnlyUsedVerticesFLAG) {
		for (i = 0; i < nFaces; i++)
			for (j = 0; j < 3; j++)
				fprintf(obj_file, "v %f %f %f\n", vertices[faces[i * 3 + j]].x,
					vertices[faces[i * 3 + j]].y, vertices[faces[i * 3 + j]].z);
	}
	else {
		for (i = 0; i < nVert; i++)
			fprintf(obj_file, "v %f %f %f\n", vertices[i].x,
				vertices[i].y, vertices[i].z);
	}

	/* export the face normals */
	for (i = 0; i < nFaces; i++) {
		/* calculate cross product between v1-v0 and v2-v0 */
		v1 = vertices[faces[i * 3 + 1]];
		v2 = vertices[faces[i * 3 + 2]];
		v1.x -= vertices[faces[i * 3]].x;
		v1.y -= vertices[faces[i * 3]].y;
		v1.z -= vertices[faces[i * 3]].z;
		v2.x -= vertices[faces[i * 3]].x;
		v2.y -= vertices[faces[i * 3]].y;
		v2.z -= vertices[faces[i * 3]].z;
		normal = cross(&v1, &v2);

		/* normalise to unit length */
		scale = 1.0 / (sqrt(pow(normal.x, 2.0) + pow(normal.y, 2.0) + pow(normal.z, 2.0)) + 2.23e-9);
		normal.x *= scale;
		normal.y *= scale;
		normal.z *= scale;
		fprintf(obj_file, "vn %f %f %f\n", normal.x, normal.y, normal.z);
	}

	/* export the face indices */
	if (keepOnlyUsedVerticesFLAG) {
		for (i = 0; i < nFaces; i++) {
			/* vertices are in same order as the faces, and normals are in order */
			fprintf(obj_file, "f %lu//%lu %lu//%lu %lu//%lu\n",
				i * 3 + 1, i + 1,
				i * 3 + 1 + 1, i + 1,
				i * 3 + 2 + 1, i + 1);
		}
	}
	else {
		/* just normals are in order  */
		for (i = 0; i < nFaces; i++) {
			fprintf(obj_file, "f %lu//%lu %lu//%lu %lu//%lu\n",
				faces[i * 3] + 1, i + 1,
				faces[i * 3 + 1] + 1, i + 1,
				faces[i * 3 + 2] + 1, i + 1);
		}
	}
	fclose(obj_file);
}

void convhull_3d_export_m
(
	ch_vertex* const vertices,
	const unsigned long nVert,
	unsigned long* const faces,
	const unsigned long nFaces,
	char* const m_filename
)
{
	unsigned long i;
	char path[256] = { "\0" };
	memcpy(path, m_filename, strlen(m_filename));
	FILE* m_file;
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
	CV_STRCAT(path, ".m");
	fopen_s(&m_file, path, "wt");
#else
	m_file = fopen(strcat(path, ".m"), "wt");
#endif

	/* save face indices and vertices for verification in matlab: */
	fprintf(m_file, "vertices = [\n");
	for (i = 0; i < nVert; i++)
		fprintf(m_file, "%f, %f, %f;\n", vertices[i].x, vertices[i].y, vertices[i].z);
	fprintf(m_file, "];\n\n\n");
	fprintf(m_file, "faces = [\n");
	for (unsigned long f = 0; f < nFaces; f++) {
		fprintf(m_file, " %lu, %lu, %lu;\n",
			faces[3 * f + 0] + 1,
			faces[3 * f + 1] + 1,
			faces[3 * f + 2] + 1);
	}
	fprintf(m_file, "];\n\n\n");
	fclose(m_file);
}

void extractVerticesFromObjFile(char* const obj_filename, ch_vertex** out_vertices, unsigned long* out_nVert)
{
	FILE* obj_file;
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
	CV_STRCAT(obj_filename, ".obj");
	fopen_s(&obj_file, obj_filename, "r");
#else
	obj_file = fopen(strcat(obj_filename, ".obj"), "r");
#endif 

	/* determine number of vertices */
	unsigned long nVert = 0;
	char line[256];
	while (fgets(line, sizeof(line), obj_file)) {
		char* vexists = strstr(line, "v ");
		if (vexists != NULL)
			nVert++;
	}
	(*out_nVert) = nVert;
	(*out_vertices) = (ch_vertex*)malloc(nVert * sizeof(ch_vertex));

	/* extract the vertices */
	rewind(obj_file);
	unsigned long i = 0;
	unsigned long vertID, prev_char_isDigit, current_char_isDigit;
	char vert_char[256] = { 0 };
	while (fgets(line, sizeof(line), obj_file)) {
		char* vexists = strstr(line, "v ");
		if (vexists != NULL) {
			prev_char_isDigit = 0;
			vertID = -1;
			for (unsigned j = 0; j < strlen(line) - 1; j++) {
				if (isdigit(line[j]) || line[j] == '.' || line[j] == '-' || line[j] == '+' || line[j] == 'E' || line[j] == 'e') {
					vert_char[strlen(vert_char)] = line[j];
					current_char_isDigit = 1;
				}
				else
					current_char_isDigit = 0;
				if ((prev_char_isDigit && !current_char_isDigit) || j == strlen(line) - 2) {
					vertID++;
					if (vertID > 4) {
						/* not a valid file */
						free((*out_vertices));
						(*out_vertices) = NULL;
						(*out_nVert) = 0;
						return;
					}
					(*out_vertices)[i].v[vertID] = atof(vert_char);
					memset(vert_char, 0, 256 * sizeof(char));
				}
				prev_char_isDigit = current_char_isDigit;
			}
			i++;
		}
	}
}