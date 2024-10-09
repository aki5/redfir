
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "redfir.h"


// 0th order modified bessel function of the first kind.
// this isn't optimized for numerical performance, but
// it seems to work well enough.
static double
bessel0(double x)
{
	double kfac = 1.0;
	double halfx = 0.5*x;
	double sum = 0.0;
	for(int k = 0; k < 20; k++){
		kfac *= k;
		if(kfac == 0)
			kfac = 1;
		double term = pow(halfx, k) / kfac;
		sum += term*term;
	}
	return sum;
}

// double att = -6.0*log2(ripple);
// there's only one set of parameters for attenuations >50dB..
// attenuation beyond 21dB doesn't impact filter order.
// what's the tradeoff between less and more attenuation then?
// more attenuation should widen the transition band.. which means
// for the same transition we should need a wider filter.. but not?

void
redfir_kaiser(RedFirKaiser *kp, double transbw, double sfreq, double att_db)
{
	double tw = 2.0 * M_PI * transbw / sfreq;

	double beta;
	if(att_db <= 21.0){
		beta = 0.0;
	} else if(att_db <= 50.0){
		beta = 0.5842 * pow(att_db-21.0, 0.4) + 0.07886 * (att_db - 21.0);
	} else {
		beta = 0.1102 * (att_db - 8.7);
	}

	kp->beta = beta;
	kp->bessel_beta = bessel0(beta);
	// make sure it's odd, we use spectral inversion.
	kp->order = 1 | (int)ceil((att_db-7.95)/(2.285*tw));
}

// needs to be tested
double
redfir_kaiser_window(void *arg, int i, int n)
{
	RedFirKaiser *kp = (RedFirKaiser*)arg;
	double a = (n-1)/2.0;
	double naa = (i-a)/a;
	double naa2 = naa*naa;
	return bessel0(kp->beta * sqrt(1.0 - naa2)) / kp->bessel_beta;
}

double
redfir_blackman_window(void *arg, int i, int n)
{
	double rcp_n = 1.0f / (n-1);
	return 0.42f - 0.5f*cos(2.0f*M_PI*i*rcp_n) + 0.08f*cos(4.0f*M_PI*i*rcp_n);
}

double
redfir_blackman_harris_window(void *arg, int i, int n)
{
	double rcp_n = 1.0f / (n-1);
	return 0.35875f - 0.48829f*cos(2.0f*M_PI*i*rcp_n) + 0.14128f*cos(4.0f*M_PI*i*rcp_n) - 0.01168f*cos(6.0f*M_PI*i*rcp_n);
}

static double
sinc(double freq_ratio, int i, int n)
{
	double shift = (n & 1) == 0 ? 0.5f : 0.0f;
	double x = (i-n/2+shift)*2.0f*M_PI*freq_ratio;
	if(x == 0.0f)
		return 1.0f;
	return sin(x)/x;
}

int
redfir_init(RedFir *fp, int len)
{
	memset(fp, 0, sizeof fp[0]);

	fp->imp = malloc(len * sizeof fp->imp[0]);
	if(fp->imp == NULL)
		return -1;

	fp->samples = malloc(2 * len * sizeof fp->samples[0]);
	if(fp->samples == NULL){
		free(fp->imp);
		return -1;
	}

	memset(fp->imp, 0, len * sizeof fp->imp[0]);
	memset(fp->samples, 0, 2 * len * sizeof fp->samples[0]);
	fp->len = len;
	return 0;
}

void
redfir_free(RedFir *fp)
{
	free(fp->imp);
	free(fp->samples);
}

void
redfir_window(RedFir *ap, double (*winfunc)(void *arg, int i, int n), void *arg)
{
	for(int i = 0; i < ap->len; i++)
		ap->imp[i] = (*winfunc)(arg, i, ap->len);
}

void
redfir_sinc(RedFir *ap, double freq_ratio)
{
	for(int i = 0; i < ap->len; i++)
		ap->imp[i] = sinc(freq_ratio, i, ap->len);
}

int
redfir_add(RedFir *dst, RedFir *ap, RedFir *bp)
{
	if(dst->len == ap->len && ap->len == bp->len){
		for(int i = 0; i < ap->len; i++)
			dst->imp[i] = ap->imp[i] + bp->imp[i];
		return 0;
	}
	return -1;
}

int
redfir_mul(RedFir *dst, RedFir *ap, RedFir *bp)
{
	if(dst->len == ap->len && ap->len == bp->len){
		for(int i = 0; i < ap->len; i++)
			dst->imp[i] = ap->imp[i] * bp->imp[i];
		return 0;
	}
	return -1;
}

int
redfir_inverse(RedFir *dst, RedFir *ap)
{
	if(dst->len == ap->len){
		int i;
		for(i = 0; i < ap->len; i++)
			dst->imp[i] = -ap->imp[i];
		dst->imp[dst->len/2] += 1.0f;
		return 0;
	}
	return -1;
}

int
redfir_reverse(RedFir *dst, RedFir *ap)
{
	if(dst->len == ap->len){
		int i;
		for(i = 0; i < ap->len; i++){
			double tmp;
			if((i&1) != 0)
				tmp = -ap->imp[i];
			else
				tmp = ap->imp[i];
			dst->imp[i] = tmp;
		}
		return 0;
	}
	return -1;
}

int
redfir_normalize(RedFir *dst, RedFir *ap)
{
	if(dst->len == ap->len){
		double sum = ap->imp[0];
		for(int i = 1; i < ap->len; i++)
			sum += ap->imp[i];
		double rcp = 1.0f / sum;
		for(int i = 0; i < ap->len; i++)
			dst->imp[i] = rcp*ap->imp[i];
		return 0;
	}
	return -1;
}

// the ring buffer is implemented by inserting each sample at two locations,
// the current offset and also the filter length past the current offset.
// the buffer needs must be twice the filter length, but there's no modulo
// arithmetic in the loop (should be easy simd)
double
redfir_step(RedFir *fp, double sample)
{
	// insert new sample in the ring buffer
	int off = fp->off;
	int len = fp->len;
	fp->samples[off] = sample;
	fp->samples[off + len] = sample;

	// advance off to oldest sample.
	// off+len-1 is the last sample we'll use, the one just added.
	// 0: .0 1 2 0 1 2   [0], [3]
	// 1:  3.1 2 3 1 2   [1], [4]
	// 2:  3 4.2 3 4 2   [2], [5]
	// 3: .3 4 5 3 4 5   [0], [3]
	// 4:  6.4 5 6 4 5   [1], [4]
	off++;
	if(off == len)
		off = 0;

	double sum = fp->samples[off] * fp->imp[0];
	for(int i = 1; i < len; i++)
		sum += fp->samples[off + i] * fp->imp[i];

	fp->off = off;

	return sum;
}

void
redfir_reset(RedFir *ap)
{
	for(int i = 0; i < 2*ap->len; i++)
		ap->samples[i] = 0.0f;
}

int
redfir_scale(RedFir *dst, RedFir *ap, double scale)
{
	if(dst->len == ap->len){
		for(int i = 0; i < ap->len; i++)
			dst->imp[i] = scale * ap->imp[i];
		return 0;
	}
	return -1;
}

void
redfir_swap(RedFir *a, RedFir *b)
{
	RedFir c;
	memcpy(&c, a, sizeof c);
	memcpy(a, b, sizeof a[0]);
	memcpy(b, &c, sizeof c);
}


