#include <stdio.h>
#include <string.h>
#include <time.h>
#include "image.h"

int main(int argc, char** argv) {
    
    time_t start = time(NULL);
    int p_num;
    char* name;
    int size;
    if(argc == 4) {
        p_num = atoi(argv[1]);
        size = atoi(argv[2]);
        name = argv[3];
    }
    else {
        printf("Please enter number of processors, size, and name of file\n");
        exit(1);
    }
    
    int width = size;
    int height = size;
    int image_size = width*height;
    int offset = width*height/p_num;
    int *m_image;

    printf("Reading in image\n");
    read_image_template(name, &m_image, &width, &height);
    
    int p_width = width+2;
    int p_height = (height/p_num)+2;
    int chunk_size = p_width*p_height;
    
    for(int p=0; p<p_num; p++) {
        
        printf("Creating Chunk: %d\n", p);
        int *chunk = (int*)malloc(sizeof(int)*chunk_size);
       
        printf("Copying image into chunk: %d\n", p); 
        for(int row=0; row<p_height; row++) {
            for(int col=0; col<p_width; col++) {

                // transposed chunk index
                int c_index = row*p_width+col;
                // transposed image index
                int m_index = (row-1)*width+(col-1);
                
                if((row == 0 && col == 0) ||
                    (row == 0 && col ==  p_width-1) ||
                    (row == p_height-1 && col == 0) ||
                    (row == p_height-1 && col == p_width-1)) {
                    //This is a corner of pixel; fill it with gray
                    chunk[c_index] = 127; 
                    continue; 
                }
                if (row == 0) {
                    if (p == 0) {
                        // duplicate bottom row
                        chunk[c_index] = m_image[row*width+(col-1)+(offset*p)];
                    }
                    else {
                        // grab the last row of the previous chunk
                        chunk[c_index] = m_image[(height/p_num-1)*width+(col-1)+(offset*(p-1))];
                    }
                }
                else if(row == p_height-1) {
                    if (p == p_num-1) {
                        // duplicate top row
                        chunk[c_index] = m_image[(row-2)*width+(col-1)+(offset*p)];
                    } 
                    else { 
                        // grab the top row of the subsequent chunk
                        chunk[c_index] = m_image[(col-1)+(offset*(p+1))];
                    }   
                }
                //Mirror edges
                else if (col == 0) {
                    chunk[row*p_width+col] = m_image[(row-1)*width+(offset*p)];
                }
                else if(col == p_width-1) {
                    chunk[row*p_width+col] = m_image[(row-1)*width+(width-1)+(offset*p)];
                }
                else {
                    int num = 0;
                    chunk[row*p_width+col] = m_image[m_index+(offset*p)];
                }
            }
        }

        char img_name[20];
        sprintf(img_name, "op_%d_%d.pbm", p, p_num);
         
        printf("Writing chunk to file: %s\n", img_name);
        write_image_template(img_name, chunk, p_width, p_height);
        free(chunk);
    }
    
    time_t end = time(NULL); 
    printf("Total time: %ld sec", end-start); 
    
    return 0;
}
