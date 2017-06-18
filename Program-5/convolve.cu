#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
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

__global__
void convolve(float** gradient, float* image, int img_w ,int img_h, float* kernel, int ker_w, int ker_h) 
{
    (*gradient) = (float*)malloc(img_h*img_w*sizeof(float));

    float sum=0;
    int i,j,k,l,m;

    i = blockIdx.x*blockDim.x + threadIdx.x;
    j = blockIdx.x*blockDim.y + threadIdx.y;

    if(i < img_h && j < width)
    { 
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


int main(int argc, char** argv) 
{
    char* img_name = argv[1];
    float sigma = atof(argv[2]);
    char name[20];

    //CPU device buffer for orginal image
    float* org_img;

    //GPU device buffer for orig image
    float* d_org_img;

    //CPU buffers for final output
    float* horizontal_gradient, *vertical_gradient;

    //GPU buffers for final output
    float* d_horizontal_gradient, *d_vertical_gradient;

    //GPU buffers to hold intermediate convolution results
    float* d_temp_horizontal, *d_temp_vertical;

    //CPU buffers to hold convolution masks    
    float* gaussian_kernel, *gaussian_dv_kernel;
    
    //GPU device buffers to store the convolution masks
    float* d_gaussian_kernel, *d_gaussian_dv_kernel;

    if(argc != 3) 
    {
        printf("./exec  <full image path> <sigma:float>");
        exit(1);
    }

    char* img_name = argv[1];
    float sigma = atof(argv[2]);
    char name[20];

    int ker_width;
    int width, height;

    gen_gaussians(&gaussian_kernel, &gaussian_dv_kernel, &ker_width, sigma);

    read_image_template<float>(img_name, &org_img, &width, &height);  
   
    cudaMalloc((void **)&d_org_img, sizeof(float)*width*height);    
    cudaMalloc((void **)&d_temp_horizontal, sizeof(float)*width*height);
    cudaMalloc((void **)&d_temp_vertical, sizeof(float)*width*height); 
    cudaMalloc((void **)&d_horizontal_gradient, sizeof(float)*width*height); 
    cudaMalloc((void **)&d_vertical_gradient, sizeof(float)*width*height); 

    cudaMalloc((void **)&d_gaussian_kernel, sizeof(float)*ker_width); 
    cudaMalloc((void **)&d_gaussian_kernel, sizeof(float)*ker_width); 

    //offload all of the data to GPU device for convolution
    cudaMemcpy(d_org_img, org_img, sizeof(float)*width*height, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gaussian_kernel, sizeof(float)*ker_width, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gaussian_dv_kernel, sizeof(flaot)*ker_width, cudaMemcpyHostToDevice);

    int block_dim = 16;
    dim3 dimGrid(ceil(height/block_dim), ceil(width/block_dim), 1);
    dim3 dimBlock(block_dim, block_dim, 1);

    convolve<<<dimGrid, dimBlock>>>(d_temp_horizontal, d_org_img, width, height, gaussian_kernel, ker_width, 1);
    convolve<<<dimGrid, dimBlock>>>(d_horizontal_gradient, d_temp_horizontal, width, height, gaussian_dv_kernel, 1, ker_width);

    convolve<<<dimGrid, dimBlock>>>(d_temp_vertical, d_org_img, width, height, gaussian_kernel, 1, ker_width);
    convolve<<<dimGrid, dimBlock>>>(d_vertical_gradient, d_temp_vertical, width, height, gaussian_dv_kernel, 1, ker_width);

    cudaMemcpy(horizontal_gradient, d_horizontal_gradient, sizeof(float)*width*height, cudaMemcpyDeviceToHost);
    cudaMemcpy(vertical_gradient, d_vertical_gradient, sizeof(float)*width*height, cudaMemcpyDeviceToHost);

    cudaThreadSynchronize();

    write_image_template<float>((char *)("horizontal_gradient.pgm"), horizontal_gradient, width*height);
    write_image_template<float>((char *)("vertical_gradient.pgm"), vertical_gradient, width*height);
   
    return 0;
}
