<F6>When image processing applications are parallelized, often times a single processor inputs the image, splits it into chunks, and preprocesses the chunks prior to distribution across the cluster. This can be viewed as the inherently sequential portion (σ) of your parallel program. In this assignment, you’ll write this sequential program and also find the time, σ.

Write a C program that:

Reads a PGM image using functions given in the header file image.h
Splits the image into p chunks (for distribution to p =2, 4, 8, 16 processors) in row-wise fashion. Note that the last chunk may have an additional row if the height_of_image%p is not zero.
Simulates Send/recv the “ghost” rows to/from the adjacent chunks.
Pads the chunks with “ghost” columns.
Writes each of the chunks to a PGM image.
NO REPORT for this assignment.
Coding Specifications:

1. The header file: image.h is given to you. Go through it and concisely comment what the different functions are doing.

2. Image given to you is 2-D, however you may use 1-D arrays to store the relevant images.

3. Save the padded chunks to files. Use the following nomenclature:

<op>_<chunk number>_<p>.pgm

For example: op_2_16.pgm is 3rd chunk (chunks start from 0) when 16 chunks are created.

4. Time your code:

oYou will time (in microseconds) the splitting and adding ghosts to your chunks. This is your inherently sequential time, σ. Insert the timer functions accordingly. Use gettimeofday function. For example, look at:

http://www.cs.loyola.edu/~ (Links to an external site.)jglenn/702/S2008/Projects/P3/time.html (Links to an external site.)

What to Turn In?

Alpha – On 01/30, commit your codes developed so far to BitBucket.

Beta -- On 02/02, commit your final code to BitBucket.
