typedef struct RedFir RedFir;
typedef struct RedFirKaiser RedFirKaiser;

struct RedFir {
	double *imp; // filter impulse
	double *samples; // ring buffer, 2*len
	int len;
	int off; // next sample index
};

// kaiser window parameters
struct RedFirKaiser {
	double beta;
	double bessel_beta;
	int order;
};

// constructor, destructor
int redfir_init(RedFir *fp, int len);
void redfir_free(RedFir *fp);

// initializers for making a window, sinc function
void redfir_kaiser(RedFirKaiser *kp, double transbw, double sfreq, double att_db);
void redfir_window(RedFir *ap, double (*winfunc)(void *arg, int i, int n), void *arg);
void redfir_sinc(RedFir *ap, double freq_ratio);

// window functions (for use with redfir_window)
double redfir_blackman_window(void *arg, int i, int n);
double redfir_blackman_harris_window(void *arg, int i, int n);
double redfir_kaiser_window(void *arg, int i, int n);

// operators, binary and unary.
int redfir_add(RedFir *dst, RedFir *ap, RedFir *bp);
int redfir_mul(RedFir *dst, RedFir *ap, RedFir *bp);
int redfir_scale(RedFir *dst, RedFir *ap, double scale);
int redfir_inverse(RedFir *dst, RedFir *ap);
int redfir_reverse(RedFir *dst, RedFir *ap);
int redfir_normalize(RedFir *dst, RedFir *ap);

// redfir_step adds a sample, executes the filter and outputs a sample.
// output sample is delayed by fp->len/2. 
double redfir_step(RedFir *fp, double sample);

// reset sample buffer to zeros.
void redfir_reset(RedFir *ap);
