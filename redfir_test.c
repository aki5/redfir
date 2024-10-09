
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "redfir.h"

void
redfir_printf(FILE *fp, RedFir *ap)
{
	for(int i = 0; i < ap->len; i++){
		fprintf(fp, "%d %.20f\n", i, ap->imp[i]);
	}
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
redfir_freqplot(FILE *fp, RedFir *ap)
{
	double ratio = powf(2.0f, 1.0f/24.0f);
	enum { InputLen = 8192 };
	double *input = malloc(InputLen * sizeof input[0]);
	for(double freq = 50.0f; freq < 22050.0f; freq += 50.0f){
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

int
main(int argc, char *argv[])
{
	FILE *fp = NULL;
	RedFir window;
	RedFir lowpass;
	RedFir highpass;
	RedFir bafflestep;

	if(argc < 4){
		fprintf(stderr, "usage: %s freq att_db bw\n", argv[0]);
		fprintf(stderr, "	freq: cutoff frequency -6dB point\n");
		fprintf(stderr, "	att_db: stop-band attenuation in dB (ie. 96 for 16bit\n");
		fprintf(stderr, "	bw: transition bandwidth from full pass to specified stop band\n");
		exit(1);
	}

	int freq = strtol(argv[1], NULL, 10);
	int att_db = strtol(argv[2], NULL, 10);
	int transbw = strtol(argv[3], NULL, 10);

	RedFirKaiser kaiser;
	redfir_kaiser(&kaiser, transbw, 44100.0, att_db);
	fprintf(stderr, "# kaiser filter for %dhz band with -%ddB stop banda\n", transbw, att_db);
	fprintf(stderr, "# filter order %d\n", kaiser.order);

	redfir_init(&window, kaiser.order);
	redfir_window(&window, &redfir_kaiser_window, (void*)&kaiser);

	redfir_init(&lowpass, window.len);
	redfir_init(&highpass, window.len);
	redfir_init(&bafflestep, window.len);

	redfir_sinc(&lowpass, freq / 44100.0f);
	redfir_mul(&lowpass, &lowpass, &window);
	redfir_normalize(&lowpass, &lowpass);
	redfir_inverse(&highpass, &lowpass);

	redfir_scale(&bafflestep, &highpass, 0.5); // -6dB
	redfir_add(&bafflestep, &bafflestep, &lowpass);

	fp = fopen("lp-freq.txt", "wb");
	redfir_freqplot(fp, &lowpass);
	fclose(fp);

	fp = fopen("hp-freq.txt", "wb");
	redfir_freqplot(fp, &highpass);
	fclose(fp);

	fp = fopen("bs-freq.txt", "wb");
	redfir_freqplot(fp, &bafflestep);
	fclose(fp);

	fp = fopen("lowpass.txt", "wb");
	redfir_printf(fp, &lowpass);
	fclose(fp);

	fp = fopen("highpass.txt", "wb");
	redfir_printf(fp, &highpass);
	fclose(fp);

	fp = fopen("bafflestep.txt", "wb");
	redfir_printf(fp, &bafflestep);
	fclose(fp);

	fp = fopen("window.txt", "wb");
	redfir_printf(fp, &window);
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
	redfir_free(&window);

	return 0;
}
