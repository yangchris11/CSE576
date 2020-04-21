#include "image.h"
#include <iostream>

using namespace std;


// SELF
const int& clamp(const int& val, const int& low, const int& high){
  return (val < low) ? low : (high < val) ? high : val;
}

// HW0 #1
// const Image& im: input image
// int x,y: pixel coordinates
// int ch: channel of interest
// returns the 0-based location of the pixel value in the data array
int pixel_address(const Image& im, int x, int y, int ch)
  {
    // TODO: calculate and return the index
    return ch * im.w * im.h + y * im.w + x; 
  }

// HW0 #1
// const Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
// returns the value of the clamped pixel at channel ch
float get_clamped_pixel(const Image& im, int x, int y, int ch)
  {
    // TODO: clamp the coordinates and return the correct pixel value
    x = clamp(x, 0, im.w-1);
    y = clamp(y, 0, im.h-1);
    ch = clamp(ch, 0, im.c-1);
    return im.data[pixel_address(im, x, y, ch)];
  }


// HW0 #1
// Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
void set_pixel(Image& im, int x, int y, int c, float value)
  {
    // TODO: Only set the pixel to the value if it's inside the image
    if(x>=0 && x<im.w && y>=0 && y<im.h && c>=0 && c<im.c){
      im.data[pixel_address(im, x, y, c)] = value;
    }
  }



// HW0 #2
// Copies an image
// Image& to: destination image
// const Image& from: source image
void copy_image(Image& to, const Image& from)
  {
    // TODO: populate the remaining fields in 'to' and copy the data
    // You might want to check how 'memcpy' function works
    // allocating data for the new image
    to.data=(float*)calloc(from.w*from.h*from.c,sizeof(float));
    to.w = from.w;
    to.h = from.h;
    to.c = from.c;
    memcpy(to.data, from.data,from.w*from.h*from.c*sizeof(float));
  }
