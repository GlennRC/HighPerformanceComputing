#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


void test_grad(float** t_grad, int width, int height, int num_sec) {
    int chk_h = height/num_sec;
    int color = 0;
    *t_grad = (float*)malloc(sizeof(float)*width*height);
    for(int i=0; i<height; i++) {
        if(i%chk_h == 0)
            color += (256/num_sec)-1;
        for(int j=0; j<width; j++) {
            (*t_grad)[i*width+j] = color;        
        }
    }
}

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

void convolve(float** gradient, float* image, int img_w, int img_h, float* kernel, int ker_w, int ker_h, int start, int end, int bounds, int debug) 
{
    (*gradient) = (float*)malloc(img_w*img_h*sizeof(float));
    
    for(int i=start; i<end; i++) 
    {
        for(int j=0; j<img_w; j++) 
        {
            float sum = 0;
            for(int k=0; k<ker_h; k++) 
            {
                for(int m=0; m<ker_w; m++) 
                {
                    int i_off = i+(-1*floor(ker_h/2)+k);
                    int j_off = j+(-1*floor(ker_w/2)+m);

                    if(i_off>=0 && i_off<bounds && j_off>=0 && j_off<img_w) { 
                        sum += image[i_off*img_w+j_off]*kernel[k*ker_w+m];
                    }
                }
            }
            (*gradient)[(i-start)*img_w+j] = sum;
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
                if(j >= 0 && j < img_w) 
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
                if(i >= 0 && i < img_h && j >= 0 && j < img_w)
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
                if(i >= 0 && i < img_h)
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
                if(i >= 0 && i < img_h && j >= 0 && j < img_w)
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

void edge_link(float** edges, float* hyst, int img_w, int img_h)
{
    int ker_w=3;
    *edges = (float*)malloc(sizeof(float)*img_w*img_h);

    for(int i=0; i<img_h; i++)
    {
        for(int j=0; j<img_w; j++)
        {
            int index = i*img_w+j;
            if(hyst[index] == 125)
            {
                int connected = 0;
                for(int k=0; k<ker_w; k++) 
                {
                    for(int m=0; m<ker_w; m++) 
                    {
                        int i_off = i+(-1*floor(ker_w/2)+k);
                        int j_off = j+(-1*floor(ker_w/2)+m);

                        if(i_off>=0 && i_off<img_h && j_off>=0 && j_off<img_w)
                        {
                            if(hyst[i_off*img_w+j_off] == 255)
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

float* mysendrcv(float** work_chunk, float* chunk, int num_ghost_rows, int chunk_width, int chunk_height, int* wrk_h, int rank, int comm_size)
{
    MPI_Status status;

    if(rank == 0 || rank==comm_size-1) {
        *wrk_h = chunk_height+num_ghost_rows;
        (*work_chunk) = (float*)malloc(sizeof(float)*chunk_width*(*wrk_h));
    }
    else {
        *wrk_h = chunk_height+2*num_ghost_rows;
        (*work_chunk) = (float*)malloc(sizeof(float)*chunk_width*(*wrk_h));
    }

    if(rank == 0)
    {
        //sendtag = senders's rank
        //rectag = source rank
        MPI_Sendrecv(chunk+(chunk_height-num_ghost_rows)*chunk_width, //sendby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank+1, //dest
                    rank, //send tag
                    (*work_chunk)+(chunk_height*chunk_width), //recvby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank+1, //source
                    rank+1, //recv tag
                    MPI_COMM_WORLD, //comm
                    &status);
    }
    else if(rank == comm_size-1)
    {
        MPI_Sendrecv(chunk, //sendby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank-1, //dest
                    rank, //send tag
                    (*work_chunk), //recvby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank-1, //source
                    rank-1, //rev tag
                    MPI_COMM_WORLD, //comm
                    &status);
    }
    else
    {
        MPI_Sendrecv(chunk, //sendby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank-1, //dest
                    rank, //send tag
                    *work_chunk, //recvby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank-1, //source
                    rank-1, //rev tag
                    MPI_COMM_WORLD, //comm
                    &status);

        MPI_Sendrecv(chunk+(chunk_height-num_ghost_rows)*chunk_width, //sendby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank+1, //dest
                    rank, //send tag
                    (*work_chunk)+((chunk_height+num_ghost_rows)*chunk_width), //recvby
                    num_ghost_rows*chunk_width, //send count
                    MPI_FLOAT, //send type
                    rank+1, //source
                    rank+1, //recv tag
                    MPI_COMM_WORLD, //comm
                    &status);
    }

    if(rank == 0)
        memcpy(*work_chunk, chunk, sizeof(float)*chunk_width*chunk_height);
    else
        memcpy((*work_chunk)+(num_ghost_rows*chunk_width), chunk, sizeof(float)*chunk_width*chunk_height);
}
