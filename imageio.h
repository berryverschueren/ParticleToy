#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <stdlib.h>
#include <png.h>

//
// COMPILATION
//
// Should be compiled with -ltiff and -lpng

// Sets tbe width and height to the appropriate values and mallocs
// a char *buffer loading up the values in row-major, RGBA format.
// The memory associated with the buffer can be deallocated with free().
// If there was an error reading file, then 0 is returned, and
// width = height = -1. 
unsigned char *loadImageRGBA(unsigned char *fileName, int *width, int *height);

// Saves image given by buffer with specicified with and height
// to the given file name, returns true on success, false otherwise.
// The image format is RGBA.
bool saveImageRGBA(char *fileName, unsigned char *buffer, int width, int height);

// returns index into image buffer for given coordinate
#define indxRGBA(X,Y,W) (((Y) * (W) + (X)) * 4)

#endif
