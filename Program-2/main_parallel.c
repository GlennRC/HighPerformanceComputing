#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include "image_template.h"

// Globals
int p_num;
float* m_image;
float* g_kernel;
float* g_dv_kernel;
float* h_gradient;
float* v_gradient;
int gbl_img_w;
int gbl_img_h;
int gbl_ker_w;
float** h_temps;
float** v_temps;

//Thread data
typedef struct
{
    pthread_t thread;
    int tid;
    int chk_h;
    int offset;
} thread_data;

void print_matrix(float* matrix, int length)
{
    printf("Matrix: ");
    for(int i=0; i<length; i++)
    {
        printf("%f ", matrix[i]);
    }
    printf("\n");
}

int in_bounds(int x, int y, int w, int h)
{
    if(x>=0 && x<w && y>=0 && y<h)
        return 1;
    else
        return 0;
}

void gen_gaussians(float** gaus, float** gaus_dv, int* width, float sigma)
{
    float a = round(2.5*sigma-0.5);
    *width = 2*a+1;
    float sum=0;

    (*gaus) = (float*)malloc((*width)*sizeof(float));
    (*gaus_dv) = (float*)malloc((*width)*sizeof(float));

    for(int i=0; i<(*width); i++)
    {
        (*gaus)[i] = exp((-1*(i-a)*(i-a)) / (2*sigma*sigma));
        sum += (*gaus)[i];
    }
    for(int i=0; i<(*width); i++)
        (*gaus)[i] /= sum;

    sum=0;
    for(int i=0; i<(*width); i++)
    {
        (*gaus_dv)[i] = -1*(i-a)*exp((-1*(i-a)*(i-a))/(2*sigma*sigma));
        sum += -i*(*gaus_dv)[i];
    }
    for(int i=0; i<(*width); i++)
        (*gaus_dv)[i] /= sum;
}

void convolve(float* gradient, float* image, int img_w ,int img_h, float* kernel, int ker_w, int ker_h, int offset)
{
    for(int i=0; i<img_h; i++)
    {
        for(int j=0; j<img_w; j++)
        {
            float sum = 0;
            for(int k=0; k<ker_h; k++)
            {
                for(int m=0; m<ker_w; m++)
                {
                    int offseti = -1*floor(ker_h/2)+k;
                    int offsetj = -1*floor(ker_w/2)+m;
                    int index = (i+offseti)*img_w+(j+offsetj)+offset;
                    if(index >= offset && index < img_w*img_h+offset)
                        sum += image[index]*kernel[k*ker_w+m];
                }
            }
            gradient[i*img_w+j+offset] = sum;
        }
    }
}
void smooth_image(int tid, int chk_h, int offset)
{
    h_temps[tid] = (float*)malloc(gbl_img_w*gbl_img_h*sizeof(float));
    v_temps[tid] = (float*)malloc(gbl_img_w*gbl_img_h*sizeof(float));

    // Convole	
	convolve(v_temps[tid], m_image, gbl_img_w, chk_h, g_kernel, 1, gbl_ker_w, offset);
    convolve(h_temps[tid], m_image, gbl_img_w, chk_h, g_kernel, gbl_ker_w, 1, offset);
}

void differentiate(int tid, int chk_h, int offset)
{
	convolve(h_gradient, h_temps[tid], gbl_img_w, chk_h, g_dv_kernel, 1, gbl_ker_w, offset);
	convolve(v_gradient, v_temps[tid], gbl_img_w, chk_h, g_dv_kernel, gbl_ker_w, 1, offset);

    free(h_temps[tid]);
    free(v_temps[tid]);
}

void* process_image(void* arg)
{   	
	thread_data* t_data = (thread_data*) arg;
    int tid = t_data->tid;
    int chk_h = t_data->chk_h;
    int offset = t_data->offset;

    smooth_image(tid, chk_h, offset);
	differentiate(tid, chk_h, offset);

    return NULL;
}

int main(int argc, char** argv) {
    
    if(argc != 4) {
        printf("./exec  <full image path> <sigma:float> <number of pthreads>");
        exit(1);
    }

    char* img_name = argv[1];
    float sigma = atof(argv[2]);
    p_num = atoi(argv[3]);
    int chk_h;
    
    // Holds the thread id and its arguments 
    thread_data *t_handles = (thread_data*) malloc(p_num*sizeof(thread_data));

    // Generate the both the gaussian and gaussian derivative kernels
    gen_gaussians(&g_kernel, &g_dv_kernel, &gbl_ker_w, sigma);
    
    // Read the image in
    read_image_template(img_name, &m_image, &gbl_img_w, &gbl_img_h); 
    
    // the horizontal and vertical gradients
    h_gradient = (float*)malloc(gbl_img_w*gbl_img_h*sizeof(float));
    v_gradient = (float*)malloc(gbl_img_w*gbl_img_h*sizeof(float));
    // array of temp images with first convolution applied
    h_temps = (float**)malloc(p_num*sizeof(float*));
    v_temps = (float**)malloc(p_num*sizeof(float*));

    // Start the timer
	struct timeval start, end;
  	gettimeofday(&start, NULL);

    // chunk height
    chk_h = gbl_img_h/p_num;

    for(int i=0; i<p_num; i++) {
        // Initialize thread and arguments
        int offset = gbl_img_w*chk_h*i;
        (t_handles[i]) = (thread_data) { .tid=i, .chk_h=chk_h, .offset=offset}; 
        pthread_create(&t_handles[i].thread, NULL, process_image, &(t_handles[i]));
    }

    // Wait for threads to finish
    for(int i=0; i<p_num; i++) 
        pthread_join(t_handles[i].thread, NULL);

    // Stop timer
	gettimeofday(&end, NULL);
	printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
		  - (start.tv_sec * 1000000 + start.tv_usec)));

    char name[20]; 
    sprintf(name, "hg_%d.pgm", p_num);
    write_image_template(name, h_gradient, gbl_img_w, gbl_img_h);

    sprintf(name, "vg_%d.pgm", p_num);
    write_image_template(name, v_gradient, gbl_img_w, gbl_img_h);	

    free(m_image);
    free(g_kernel);
    free(g_dv_kernel);
    free(t_handles);
    free(h_gradient);
    free(v_gradient);

    return 0;
}
