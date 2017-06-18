#include <sys/time.h>
#include "image_template.h"
#include "convolve.h"

int main(int argc, char** argv) 
{
    if(argc != 3) 
    {
        printf("./exec  <full image path> <sigma:float>");
        exit(1);
    }

    //Setup MPI
    int comm_size,comm_rank,rc;
    rc=MPI_Init(&argc, &argv);

    if(rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);

    //get communicator size and rank
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    
    char* img_name = argv[1];
    float sigma = atof(argv[2]);
    char name[20];
    float* image;
    float* chunk;
    float* work_chk;
    float* temp;
    float* g_ker;
    float* g_dv_ker;
    float* h_grad;
    float* v_grad;
    float* mag_img;
    float* phase_img;
    float* supp_img;
    float* hyst_img;
    float* edge_img;
    int num_ghst_rows;
    int ker_w, img_w, img_h, chk_h, wrk_h;
     
    gen_gaussians(&g_ker, &g_dv_ker, &ker_w, sigma);

    num_ghst_rows = ker_w/2;
  
	struct timeval start, end;
    if(comm_rank==0) {
    // Start the timer
    	gettimeofday(&start, NULL);
        read_image_template(img_name, &image, &img_w, &img_h);
    }
    
    MPI_Bcast(&img_w, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&img_h, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Set the hieght of the image chunk
    chk_h = img_h/comm_size;
    
    // The receiving chunk of image for each proc
    chunk = (float*)malloc(img_w*chk_h*sizeof(float));	

    // Send the chunks of image to each proc
    MPI_Scatter(    image,      //sendby
                    img_w*chk_h, //send count
                    MPI_FLOAT,  //send type
                    chunk,  //recvby
                    img_w*chk_h,  //recv count
                    MPI_FLOAT,  //send type
                    0,          //root
                    MPI_COMM_WORLD); //comm   

	int offset = num_ghst_rows;
    if(comm_rank==0) offset = 0;

    mysendrcv(&work_chk, chunk, num_ghst_rows, img_w, chk_h, &wrk_h, comm_rank, comm_size); 

	//Horizontal gradient
    convolve(&temp, work_chk, img_w, wrk_h, g_ker, ker_w, 1, 0, wrk_h, wrk_h, 0);    
    convolve(&h_grad, temp, img_w, chk_h, g_dv_ker, 1, ker_w, offset, offset+chk_h, wrk_h, 0);
    free(temp);

	//Veritical gradient 
    convolve(&temp, work_chk, img_w, wrk_h, g_ker, 1, ker_w, offset, wrk_h, wrk_h, 0);
    convolve(&v_grad, temp, img_w, chk_h, g_dv_ker, 1, ker_w, offset, offset+chk_h, wrk_h, 0);

    magnitude_phase(&mag_img, &phase_img, h_grad, v_grad, img_w, chk_h);

    suppression(&supp_img, mag_img, phase_img, img_w, chk_h);

    edge(&hyst_img, supp_img, img_w, chk_h);

    edge_link(&edge_img, hyst_img, img_w, chk_h);

    float* chunks = NULL;

    if(comm_rank == 0)
        chunks = (float*)malloc(sizeof(float)*img_w*img_h);

    MPI_Gather( edge_img,
                img_w*chk_h,
                MPI_FLOAT,
                chunks,
                img_w*chk_h,
                MPI_FLOAT,
                0,
                MPI_COMM_WORLD);

    if(comm_rank == 0) {
    	gettimeofday(&end, NULL);
    	printf("%d, %d, %ld\n", comm_size, img_w, ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
        sprintf(name, "edge.pgm");
        write_image_template(name, chunks, img_w, img_h);
    }
    

    MPI_Finalize();
 
    return 0;
}
