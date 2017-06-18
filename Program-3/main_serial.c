#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "image_template.h"

// helper function for printing out the matrix
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

void magnitude_phase(float** mag, float** phase, float* h_grad, float* v_grad, int width, int height)
{
    (*mag) = (float*)malloc(width*height*sizeof(float));
    (*phase) = (float*)malloc(width*height*sizeof(float));

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
			int index = i*width+j;

         	(*mag)[index]=sqrt( pow(v_grad[index],2) + pow(h_grad[index],2));
      		(*phase)[index]=atan2(v_grad[index], h_grad[index]);
        }
    }
}

// generate the gaussian kernal masks
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
        sum += (*gaus)[i-1];
    }
    for(int i=0; i<(*width); i++) 
        (*gaus)[i] /= sum;

    for(int i=0; i<(*width); i++) 
    {
        (*gaus_dv)[i] = -1*(i-a)*exp((-1*(i-a)*(i-a))/(2*sigma*sigma));
        sum += -i*(*gaus_dv)[i-1];
    }
    for(int i=0; i<(*width); i++)
        (*gaus_dv)[i] /= sum;
}

void convolve(float** gradient, float* image, int img_w ,int img_h, float* kernel, int ker_w, int ker_h) 
{
    (*gradient) = (float*)malloc(img_h*img_w*sizeof(float));

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

                    if(in_bounds(i+offseti, j+offsetj, img_w, img_h)) 
                        sum += image[(i+offseti)*img_w+(j+offsetj)]*kernel[k*ker_w+m];
                }
            }
            (*gradient)[i*img_w+j] = sum;
        }
    }
}

int main(int argc, char** argv) 
{
    if(argc != 3) 
    {
        printf("./exec  <full image path> <sigma:float>");
        exit(1);
    }

    char* img_name = argv[1];
    float sigma = atof(argv[2]);
    float* image;
    float* gaussian_kernel;
    float* gaussian_dv_kernel;
    char g_img_name[20];
    float* horizontal_gradient;
    float* vertical_gradient;
    float* mag_image;
    float* phase_image;
    float* temp;
    int width;
    int size;

    printf("Generating Gaussian kernels...\n");
    gen_gaussians(&gaussian_kernel, &gaussian_dv_kernel, &width, sigma);

    print_matrix(gaussian_kernel, width);
    print_matrix(gaussian_dv_kernel, width);
 
    printf("Reading image: %s size: %d\n", img_name, size);
    read_image_template(img_name, &image, &size, &size);  
   
  	// Start the timer
    struct timeval start, end;
    gettimeofday(&start, NULL);
  
    printf("Convoluting image to horizontal gradient...\n");
    convolve(&temp, image, size, size, gaussian_kernel, width, 1);
    convolve(&horizontal_gradient, temp, size, size, gaussian_dv_kernel, 1, width);
    free(temp);

    printf("Convoluting image to veritical gradient...\n");
    convolve(&temp, image, size, size, gaussian_kernel, 1, width);
    convolve(&vertical_gradient, image, size, size, gaussian_dv_kernel, 1, width);
    free(temp);

    // Stop timer
    gettimeofday(&end, NULL);
    printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
          - (start.tv_sec * 1000000 + start.tv_usec)));

    sprintf(g_img_name, "hg.pgm");
    printf("Writing horizontal gradient: %s, size: %d\n", g_img_name, size);
    write_image_template(g_img_name, horizontal_gradient, size, size);

    sprintf(g_img_name, "vg.pgm");
    printf("Writing vertical gradient: %s, size: %d\n", g_img_name, size);
    write_image_template(g_img_name, vertical_gradient, size, size);

    magnitude_phase(&mag_image, &phase_image, horizontal_gradient, vertical_gradient, size, size);

    char mag_name[20];
    sprintf(mag_name, "mag_name.pgm");
    printf("Writing magnitude image\n");
    write_image_template(mag_name, mag_image, size, size);

    char phase_name[20];
    sprintf(phase_name, "phase_name.pgm");
    printf("Writing phase image\n");
    write_image_template(phase_name, phase_image, size, size);


    free(image);
    free(gaussian_kernel);
    free(gaussian_dv_kernel);
    free(horizontal_gradient);
    free(vertical_gradient);

    return 0;
}
