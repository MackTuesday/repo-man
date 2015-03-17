/*
 *
 *  Generalized Iterated Function System
 *  Header
 *
 *  Brent Lehman
 *  1 February 2015
 *
 *
 *
 */

//#ifndef GeneralizedIFSC_h
//#define GeneralizedIFSC_h

#include <stdlib.h>


#define MAX_PARAMETERS      24
#define MAX_FUNCTIONS        8
#define MAX_CACHE_SIZE 1048576
#define BITS_PER_CHANNEL     8

static unsigned int brentRandN = 427856372;
#define B_RAND_MAX 4294967296
#define B_RAND() (brentRandN = brentRandN*1103515425 + 12345)
#define SELECT(n) (((B_RAND()>>16)*(n)) >> 16)


enum
{
    IF_QUADRATIC,
    IF_NUMTYPES
};

enum
{
    MAP_IDENTITY,       //  0
    MAP_SINE,           //  1
    MAP_SPHERICAL,      //  2
    MAP_SWIRL,          //  3
    MAP_HORSESHOE,      //  4
    MAP_POLAR,          //  5
    MAP_HANDKERCHIEF,   //  6
    MAP_HEART,          //  7
    MAP_DISC,           //  8
    MAP_SPIRAL,         //  9
    MAP_HYPERBOLIC,     // 10
    MAP_DIAMOND,        // 11
    MAP_EX,             // 12
    MAP_JULIA,          // 13
    MAP_BENT,           // 14
    MAP_WAVES,          // 15
    MAP_FISHEYE,        // 16
    MAP_POPCORN,        // 17
    MAP_EXPONENTIAL,    // 18
    MAP_POWER,          // 19
    MAP_COSINE,         // 20
    MAP_RINGS,          // 21
    MAP_FAN,            // 22
    MAP_BLOB,           // 23
    MAP_EYEFISH,        // 24
    MAP_BUBBLE,         // 25
    MAP_CYLINDER,       // 26
    MAP_CURL,           // 27
    MAP_RECTANGLES,     // 28
    MAP_TANGENT,        // 29
    MAP_SECANT,         // 30
    MAP_NUMTYPES
};

struct gifs_buffer_t
{
    char*    data;
    unsigned size;
};

struct gifs_t
{
    int    num_functions;
    int*   function_types;
    int*   map_types;
    float* param_sets;
    int*   probs;
    int*   color_sets;

    int    cache_size;
    int*   sequence;
    float* x_cache;
    float* y_cache;

    float  x_last;
    float  y_last;

    struct gifs_buffer_t buffer;
};

void  gifs_safe_free(void* ptr);

void  gifs_init(struct gifs_t* ifsptr, size_t cache_size);
void  gifs_cleanup(struct gifs_t* ifsptr);
int   gifs_perturb(struct gifs_t* ifsptr);
int   gifs_compute(struct gifs_t* ifsptr, int count);
void  gifs_imitate(struct gifs_t* ifsptr, const struct gifs_t* srcptr);
float gifs_gen_new_param(int param_index);

struct gifs_buffer_t* gifs_export_buffer(struct gifs_t* ifsptr);
int gifs_import_buffer(struct gifs_t* ifsptr, void* buffer);

//#endif
