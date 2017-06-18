#include <stdio.h>
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

void suppression(float** sup, float* mag, float* phase, int img_w, int img_h)
{
    (*sup) = (float*)malloc(sizeof(float)*img_h*img_w);
    memcpy(*sup, mag, sizeof(float)*img_h*img_w);
    // for each Gmag
    for(int i=0; i<img_h; i++)
    {
        for(int j=0; j<img_w; j++)
        {
            int index = i*img_w+j;
            float theta = phase[index];

            if(theta < 0)
                theta += M_PI;
            theta *= (180/M_PI);
            
            //compare Gphase with left and right of Gmag
            if(theta <= 22.5 || theta > 157.5)
            {
                //left and right
                if(j > 0 && j < img_w-1) 
                {
                    int left = i*img_w+(j-1);
                    int right = i*img_w+(j+1);
                    if(mag[index] < mag[left] && mag[index] < mag[right])
                        (*sup)[index] = 0;          
                }
            }
            else if(theta > 22.5 && theta <= 67.5)
            {
                //top-left bottom-right
                if(i > 0 && i < img_h-1 && j > 0 && j < img_w-1)
                {
                    int topleft = (i-1)*img_w+(j-1);
                    int bottomright = (i+1)*img_w+(j+1);
                    if(mag[index] < mag[topleft] && mag[index] < mag[bottomright])
                        (*sup)[index] = 0; 
                }
            }
            else if(theta > 67.5 && theta <= 112.5)
            {
                //top bottom
                if(i > 0 && i < img_h-1)
                {
                    int top = (i-1)*img_w+j;
                    int bottom = (i+1)*img_w+j;
                    if(mag[index] < mag[top] && mag[index] < mag[bottom])
                        (*sup)[index] = 0; 
                }
            }
            else if(theta > 112.5 && theta <= 157.5)
            {
                //top-right bottom-left
                if(i > 0 && i < img_h-1 && j > 0 && j < img_w-1)
                {
                    int topright = (i-1)*img_w+(j+1);
                    int bottomleft = (i+1)*img_w+(j-1);
                    if(mag[index] < mag[topright] && mag[index] < mag[bottomleft])
                        (*sup)[index] = 0; 
                }
            }

        }
    }
}

int cmpfunc(const void * a, const void * b)
{
    float x = *(float*)a;
    float y = *(float*)b;
    return x < y ? -1 : x == y ? 0 : 1;
}

void edge(float** hyst, float* sup, int img_h, int img_w)
 {
    size_t size = img_h*img_w*sizeof(float);
    (*hyst) = (float*)malloc(size);

    float* temp = (float*)malloc(size);
    memcpy(temp, sup, size);

	int len = img_w*img_h;
    qsort(temp, len, sizeof(float), cmpfunc);

	float h = len*.9;
	float l = h/5;	

    float t_high = temp[(int)h];
    float t_low = temp[(int)l];


    for(int i=0; i<len; i++) 
    {
        if(sup[i] >= t_high)
            (*hyst)[i] = 255;
        else if(sup[i] <= t_low)
            (*hyst)[i] = 0;
        else if(sup[i] > t_low && sup[i] < t_high)
            (*hyst)[i] = 125;
    }
}

void edge_link(float** edges, float* hyst, int img_h, int img_w)
{
    int len = img_h*img_w;
    int size = len*sizeof(float);

    *edges = (float*)malloc(size);

    for(int i=0; i<img_h; i++)
    {
        for(int j=0; j<img_w; j++)
        {
            int index = i*img_w+j;
            if(hyst[index] == 125)
            {
                int connected = 0;
                for(int k=-2; k<3; k++) 
                {
                    for(int l=-2; l<3; l++)
                    {
                        if(k==0 && l==0) continue;
                        if(i+k >= 0 && i+k < img_h && j+l >= 0 && j+l < img_w)
                        {
                            if(hyst[(i+k)*img_w+(j+l)] == 255)
                            {
                                (*edges)[index] = 255;
                                connected = 1;
                                break;   
                            }    
                        }
                    }
                    if(connected == 1) break;
                }
                if(connected == 0)
                    (*edges)[index] = 0;
            }
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
    char name[20];
    float* image;
    float* gaussian_kernel;
    float* gaussian_dv_kernel;
    float* horizontal_gradient;
    float* vertical_gradient;
    float* mag_image;
    float* phase_image;
    float* sup_image;
    float* hyst_image;
    float* edge_image;
    float* temp;
    int width;
    int size;

    gen_gaussians(&gaussian_kernel, &gaussian_dv_kernel, &width, sigma);

    read_image_template(img_name, &image, &size, &size);  
   
  	// Start the timer
    struct timeval start, end;
    gettimeofday(&start, NULL);
  
    convolve(&temp, image, size, size, gaussian_kernel, width, 1);
    convolve(&horizontal_gradient, temp, size, size, gaussian_dv_kernel, 1, width);
    free(temp);

    convolve(&temp, image, size, size, gaussian_kernel, 1, width);
    convolve(&vertical_gradient, image, size, size, gaussian_dv_kernel, 1, width);
    
    // Stop timer
    gettimeofday(&end, NULL);
    printf("1, %d, %ld\n", size, ((end.tv_sec * 1000000 + end.tv_usec)
          - (start.tv_sec * 1000000 + start.tv_usec)));
  
    free(temp);
    free(image);
    free(gaussian_kernel);
    free(gaussian_dv_kernel);
    free(horizontal_gradient);
    free(vertical_gradient);

    return 0;
}
