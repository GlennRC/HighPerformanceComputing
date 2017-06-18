#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "image_template.h"
#include <math.h>
#include <sys/time.h>
#define TILE_DIM 16

//A consolidated function that creates both the Gaussian kernel and the derivative kernel
void create_gaussians(float **gaussian_kernel,float **gaussian_deriv,int k_width,float sigma)
{
    int i;
    float sum=0;
    int a=k_width/2;

    *gaussian_kernel=(float *)malloc(sizeof(float)*k_width);   
    *gaussian_deriv=(float *)malloc(sizeof(float)*k_width);  

    //Create kernel
    sum=0;
    for(i=0;i<k_width;i++)
    {
        (*gaussian_kernel)[i]=exp((-1*(i-a)*(i-a))/(2*sigma*sigma));  
        sum+=(*gaussian_kernel)[i];
    }

    for(i=0;i<k_width;i++)
        (*gaussian_kernel)[i]/=sum;


    //Create derivative
    sum=0;
    for(i=0;i<k_width;i++)
    {
        (*gaussian_deriv)[i]=-1*(i-a)*exp((-1*(i-a)*(i-a))/(2*sigma*sigma));  
        sum-=i*((*gaussian_deriv)[i]);
    }

    for(i=0;i<k_width;i++)
        (*gaussian_deriv)[i]/=sum;
}
__global__
void convolve_v(float *in_image, int width, int height, float *mask, int mask_width, float *out_image)
{
    extern __shared__ float Ashared[];

    float sum;
    int offsetj;

    int i=blockIdx.x*blockDim.x+threadIdx.x;
    int j=blockIdx.y*blockDim.y+threadIdx.y;
    
    //Ashared[threadIdx.x*blockDim.x+threadIdx.y] = in_image[i*width+j];
    Ashared[threadIdx.x*blockDim.x+threadIdx.y] = (float)in_image[i*width+j];
    
    __syncthreads();
    sum=0;
    for(int k=0;k<mask_width;k++)
    {
        offsetj = (-1 * floor((double)(mask_width/2)) + k);
        if(threadIdx.y+offsetj>=0 && threadIdx.y+offsetj<blockDim.y)
        {
            sum+=(float)Ashared[threadIdx.x*blockDim.x+threadIdx.y+offsetj]*mask[k];
        }
        else //have faith that global accesses core cache
        {
            if(j+offsetj>=0 && j+offsetj<width)
            {
                sum+= in_image[i * width + j + offsetj] * mask[k];
            }
        }
    }
    __syncthreads();
    out_image[i*width+j]=(float)sum; 
}

__global__
void convolve_h(float *in_image, int width, int height, float *mask, int mask_width, float *out_image)
{
    extern __shared__ float Ashared[];

    float sum;
    int offseti;

    int i=blockIdx.x*blockDim.x+threadIdx.x;
    int j=blockIdx.y*blockDim.y+threadIdx.y;
    
    //Ashared[threadIdx.x*blockDim.x+threadIdx.y] = in_image[i*width+j];
    Ashared[threadIdx.x*blockDim.x+threadIdx.y] = (float)in_image[i*width+j];
    
    __syncthreads();
    sum=0;
    for(int k=0;k<mask_width;k++)
    {
        offseti = (-1 * floor((double)(mask_width/2)) + k);
        if(threadIdx.x+offseti>=0 && threadIdx.x+offseti<blockDim.x)
        {
            sum+=(float)Ashared[(threadIdx.x + offseti) * blockDim.x +threadIdx.y]*mask[k];
        }
        else //have faith that global accesses core cache
        {
            if(i+offseti>=0 && i+offseti<height) //Have faith in L2 but make sure data is there!
                sum+= in_image[(i + offseti) * width + j] * mask[k];
        }

    }
    __syncthreads();
    out_image[i*width+j]=(float)sum; 
}

int main(int argc, char **argv)
{

    //Declare all of the variable here
    float  *org_img;

    //GPU device buffer for original image
    float *d_org_img;

    //CPU host buffers for the final output 
    float  *vertical_gradient, *horizontal_gradient;

    //GPU buffers for the final result
    float *d_vertical_gradient, *d_horizontal_gradient;

    //GPU buffers to hold intermediate convolution results
    float *d_temp_horizontal, *d_temp_vertical;

    //CPU host buffers to store the convolution masks
    float *gaussian_kernel, *gaussian_deriv;

    //GPU device buffers to store the convolution masks
    float *d_gaussian_kernel, *d_gaussian_deriv;

    int width,height,k_width;
    float sigma,a;
    struct timeval start,end;
    
    if(argc!=3)
    {
        printf("\n The correct argument list is: exec <image file> <Sigma> \n");
        exit(0);
    }

    //obtain the parameters
    sigma=atof(argv[2]);
    a=ceil((float)(2.5*sigma-0.5));
    k_width=2*a+1;

    //CPU portion of the code that reads/prepares the input data
    read_image_template<float>(argv[1],&org_img,&width,&height);    

    //Computation starts here
    gettimeofday(&start,NULL);

    create_gaussians(&gaussian_kernel,&gaussian_deriv,k_width,sigma);

    horizontal_gradient=(float *)malloc(sizeof(float)*width*height);
    vertical_gradient=(float *)malloc(sizeof(float)*width*height);

    // CPU host mallocs for GPU buffers
    cudaMalloc((void **)&d_org_img,sizeof(float)*width*height);
    cudaMalloc((void **)&d_temp_horizontal,sizeof(float)*width*height);
    cudaMalloc((void **)&d_temp_vertical,sizeof(float)*width*height);
    cudaMalloc((void **)&d_horizontal_gradient,sizeof(float)*width*height);
    cudaMalloc((void **)&d_vertical_gradient,sizeof(float)*width*height);

    cudaMalloc((void **)&d_gaussian_kernel,sizeof(float)*k_width);
    cudaMalloc((void **)&d_gaussian_deriv,sizeof(float)*k_width);

    //Offload all of the data to GPU device for convolution
    cudaMemcpy(d_org_img,org_img,sizeof(float)*width*height,cudaMemcpyHostToDevice);

    cudaMemcpy(d_gaussian_kernel,gaussian_kernel,sizeof(float)*k_width,cudaMemcpyHostToDevice);
    cudaMemcpy(d_gaussian_deriv,gaussian_deriv,sizeof(float)*k_width,cudaMemcpyHostToDevice);


    //Horizontal gradient. vertical kernel then horizontal derivative
    int block_dim=16;
    dim3 dimGrid(ceil(height/block_dim),ceil(width/block_dim),1);
    dim3 dimBlock(block_dim,block_dim,1);

    convolve_v<<<dimGrid,dimBlock, sizeof(float)*TILE_DIM*TILE_DIM>>>(d_org_img,width,height,d_gaussian_kernel,k_width,d_temp_horizontal); 
	convolve_h<<<dimGrid,dimBlock, sizeof(float)*TILE_DIM*TILE_DIM>>>(d_temp_horizontal,width,height,d_gaussian_deriv,k_width,d_horizontal_gradient);
    
    //Vertical gradient. horizontal kernel then vertical derivative
    convolve_h<<<dimGrid,dimBlock, sizeof(float)*TILE_DIM*TILE_DIM>>>(d_org_img,width,height,d_gaussian_kernel,k_width,d_temp_vertical); 
	convolve_v<<<dimGrid,dimBlock, sizeof(float)*TILE_DIM*TILE_DIM>>>(d_temp_vertical,width,height,d_gaussian_deriv,k_width,d_vertical_gradient);

    //GPU to Host transfer of the final result
    cudaMemcpy(horizontal_gradient,d_horizontal_gradient,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
    cudaMemcpy(vertical_gradient,d_vertical_gradient,sizeof(float)*width*height,cudaMemcpyDeviceToHost);

    write_image_template<float>((char *)("horizontal_gradient.pgm"),horizontal_gradient,width,height);
    write_image_template<float>((char *)("horizontal_gradient.pgm"),vertical_gradient,width,height);

    cudaThreadSynchronize();

    gettimeofday(&end,NULL);

    printf("%d, %ld\n", width, ((end.tv_sec * 1000 + end.tv_usec/1000)
          - (start.tv_sec * 1000 + start.tv_usec/1000)));

    //free variables
    free(org_img);
    free(horizontal_gradient);
    free(vertical_gradient);
    free(gaussian_kernel);
    free(gaussian_deriv);

    cudaFree(d_org_img);
    cudaFree(d_gaussian_kernel);
    cudaFree(d_gaussian_deriv);
    cudaFree(d_temp_horizontal);
    cudaFree(d_temp_vertical);
    cudaFree(d_vertical_gradient);
    cudaFree(d_horizontal_gradient);
    return 0;
}
