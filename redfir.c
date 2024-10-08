
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

typedef struct RedFir RedFir;
struct RedFir {
	double *imp; // filter impulse
	double *samples; // ring buffer, 2*len
	int len;
	int off; // next sample index
};

static double
blackman_window(int i, int n)
{
	double rcp_n = 1.0f / (n-1);
	return 0.42f - 0.5f*cos(2.0f*M_PI*i*rcp_n) + 0.08f*cos(4.0f*M_PI*i*rcp_n);
}

static double
blackman_nuttall_window(int i, int n)
{
	double a0 = 0.3635819f;
	double a1 = 0.4891775f;
	double a2 = 0.1365995f;
	double a3 = 0.0106411f;
	double rcp_n = 1.0f / (n-1);
	double wn = a0 - a1 * cos(2.0*M_PI*i*rcp_n) + a2 * cos(4.0f*M_PI*i*rcp_n) - a3 * cos(6.0*M_PI*i*rcp_n);
	return wn;
}

static double
blackman_harris_window(int i, int n)
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
redfir_window(RedFir *ap, double (*winfunc)(int i, int n))
{
	for(int i = 0; i < ap->len; i++)
		ap->imp[i] = (*winfunc)(i, ap->len);
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

// [ 0 1 2 3 4 0 1 2 3 4 ]
//   x y z w v x y z w v
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

void
redfir_genfreq(double *dst, int len, double freqratio)
{
	for(int i = 0; i < len; i++){
		double x = i*2.0f*M_PI*freqratio;
		dst[i] = sin(x);
	}
}

void
redfir_scale(RedFir *dst, RedFir *ap, double scale)
{
	for(int i = 0; i < ap->len; i++)
		dst->imp[i] = scale * ap->imp[i];
}

void
redfir_swap(RedFir *a, RedFir *b)
{
	RedFir c;
	memcpy(&c, a, sizeof c);
	memcpy(a, b, sizeof a[0]);
	memcpy(b, &c, sizeof c);
}

void
redfir_freqplot(FILE *fp, RedFir *ap)
{
	double ratio = powf(2.0f, 1.0f/24.0f);
	enum { InputLen = 8192 };
	double *input = malloc(InputLen * sizeof input[0]);
//	for(double freq = 13.75f; freq < 10000.0f; freq *= ratio){
	for(double freq = 20.0f; freq < 10000.0f; freq += 25.0f){
		redfir_genfreq(input, InputLen, freq / 44100.0f);
		redfir_reset(ap);
		double maxout = 0.0f;
		for(int i = 0; i < InputLen; i++){
			double out = redfir_step(ap, input[i]);
			if(i > ap->len && maxout < out)
				maxout = out;
		}
		fprintf(fp, "%f %f\n", freq, 6.0*log2f(maxout));
	}
}

void
redfir_printf(FILE *fp, RedFir *ap)
{
	for(int i = 0; i < ap->len; i++){
		fprintf(fp, "%d %.6f\n", i, ap->imp[i]);
	}
}

int
main(int argc, char *argv[])
{
	RedFir win;
	RedFir lowpass;
	RedFir highpass;

	int freq = strtol(argv[1], NULL, 10);
	int flen = strtol(argv[2], NULL, 10);

	redfir_init(&lowpass, flen);
	redfir_init(&highpass, flen);
	redfir_init(&win, flen);

	//redfir_window(&win, &blackman_nuttall_window);
	redfir_window(&win, &blackman_harris_window);
	redfir_sinc(&lowpass, freq / 44100.0f);
	redfir_mul(&lowpass, &lowpass, &win);
	redfir_normalize(&lowpass, &lowpass);
	redfir_inverse(&highpass, &lowpass);

	double att_db = 96.0;
	fprintf(stderr, "# estimated order %.1f\n", 44100.0f / freq);

	FILE *fp = fopen("lp-freq.txt", "wb");
	redfir_freqplot(fp, &lowpass);
	fclose(fp);

	fp = fopen("hp-freq.txt", "wb");
	redfir_freqplot(fp, &highpass);
	fclose(fp);

	fp = fopen("lowpass.txt", "wb");
	redfir_printf(fp, &lowpass);
	fclose(fp);

	fp = fopen("window.txt", "wb");
	redfir_printf(fp, &win);
	fclose(fp);

	// generate step response and spit it out.

	fp = fopen("step.txt", "wb");
	double input = 1.0f;
	int j = 0;
	for(int i = 0; i < 1000; i++){
		if(j++ == 200){
			input = -input;
			j = 0;
		}
		double output = redfir_step(&lowpass, input);
		fprintf(fp, "%d %.3f\n", i, output);
	}
	fclose(fp);

	redfir_free(&lowpass);
	redfir_free(&highpass);
	redfir_free(&win);

	return 0;
}
