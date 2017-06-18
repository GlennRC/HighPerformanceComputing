/* This program was originally written by
Sumedh Naik (now at Intel) at Clemson University
as a part of his thesis titled, "Connecting Architectures,
fitness, Optimizations and Performance using an Anisotropic
Diffusion Filter. This header was also used
in Dr. Pallipuram's dissertation work. */

#include "stdio.h" 
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define BUFFER 512

// Function Declaration


void read_image(char *name, unsigned char **image, int *im_width, int *im_height);
void read_image_template(char *name, int **image, int *im_width, int *im_height);
void write_image(char *name, unsigned char *image, int im_width, int im_height);
void write_image_template(char *name, int *image, int im_width, int im_height);

//Function Definition

/*Call this function alone to read images*/

void read_image_template(char *name, int **image, int *im_width, int *im_height)
{
        unsigned char *temp_img;

        int i;
        
        // Reads in the image that can be accessed by temp_img
        read_image(name, &temp_img, im_width, im_height);

        // Allocates space for a 2D array
        *image=(int *)malloc(sizeof(int)*(*im_width)*(*im_height));

        // Copies temp_img to image
        for(i=0;i<(*im_width)*(*im_height);i++)
        {
                // typecaststhe unsigned chars to ints
                (*image)[i]=(int)temp_img[i];
        }

        //Frees the temp image
        free(temp_img);
}

void read_image(char *name, unsigned char **image, int *im_width, int *im_height)
{
	FILE *fip;
	char buf[BUFFER];
    // pointer to tokenized image
	char *parse;
	int im_size;
	
    // open the file
	fip=fopen(name,"rb");
    // check if file pointer is null
	if(fip==NULL)
	{
		fprintf(stderr,"ERROR:Cannot open %s\n",name);
		exit(0);
	}
    // read in line from image file
	fgets(buf,BUFFER,fip);
	do
	{
		fgets(buf,BUFFER,fip);
	}
	while(buf[0]=='#');
	parse=strtok(buf," ");
	(*im_width)=atoi(parse);

	parse=strtok(NULL,"\n");
	(*im_height)=atoi(parse);

	fgets(buf,BUFFER,fip);
	parse=strtok(buf," ");
	
	im_size=(*im_width)*(*im_height);
	(*image)=(unsigned char *)malloc(sizeof(unsigned char)*im_size);
	fread(*image,1,im_size,fip);
	
	fclose(fip);
}

/* Call this function alone to write an image*/

void write_image_template(char *name, int *image, int im_width, int im_height)
{
	    int i;
        // create a temp image that will be written
        unsigned char *temp_img=(unsigned char*)malloc(sizeof(unsigned char)*im_width*im_height);
        for(i=0;i<(im_width*im_height);i++)
        {
                // typecast int to char
                temp_img[i]=(unsigned char)image[i];
        }
        // write image
        write_image(name,temp_img,im_width,im_height);

        free(temp_img);
}


void write_image(char *name, unsigned char *image, int im_width, int im_height)
{
	FILE *fop; 
	int im_size=im_width*im_height;
	
	fop=fopen(name,"w+");
	fprintf(fop,"P5\n%d %d\n255\n",im_width,im_height);
	fwrite(image,sizeof(unsigned char),im_size,fop);
	
	//fclose(fop);
}

