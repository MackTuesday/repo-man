#include           <float.h>
#include           <stdio.h>
#include          <string.h>
#include            <math.h>
#include            <time.h>
#include "GeneralizedIFSC.h"

#pragma warning(disable:4996)

//#define SELECT(n) (((B_RAND()>>16)*((unsigned)(n))) >> 16)
#define RANDFLOAT() ((float)(B_RAND() / (float)B_RAND_MAX))

void init(struct gifs_t* ifsptr);
void compute_buffer_size(struct gifs_t* ifsptr);

void gifs_safe_free(void* ptr)
{
    if (NULL != ptr)
    {
        free(ptr);
    }
}

void gifs_init(struct gifs_t* ifsptr, size_t cache_size)
{
    size_t malloc_size = 0;

    ifsptr->cache_size = (0 != cache_size) ? cache_size : MAX_CACHE_SIZE;

    malloc_size = sizeof(int) * MAX_FUNCTIONS;
    ifsptr->function_types = (int*)malloc(malloc_size);
    memset(ifsptr->function_types, 0, malloc_size);

    malloc_size = sizeof(int) * MAX_FUNCTIONS;
    ifsptr->map_types = (int*)malloc(malloc_size);
    memset(ifsptr->map_types, 0, malloc_size);

    malloc_size = sizeof(float) * MAX_FUNCTIONS * MAX_PARAMETERS;
    ifsptr->param_sets = (float*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    malloc_size = sizeof(int) * MAX_FUNCTIONS;
    ifsptr->probs = (int*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    malloc_size = sizeof(int) * MAX_FUNCTIONS * 3;
    ifsptr->color_sets = (int*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    malloc_size = sizeof(int) * ifsptr->cache_size;
    ifsptr->sequence = (int*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    malloc_size = sizeof(float) * ifsptr->cache_size;
    ifsptr->x_cache = (float*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    malloc_size = sizeof(float) * ifsptr->cache_size;
    ifsptr->y_cache = (float*)malloc(malloc_size);
    memset(ifsptr->param_sets, 0, sizeof(malloc_size));

    init(ifsptr);
}

void init(struct gifs_t* ifsptr)
{
    int i, j;

    ifsptr->num_functions = 3;

    for (i = 0; i < MAX_FUNCTIONS; i++)
    {
        ifsptr->function_types[i] = IF_QUADRATIC;
        ifsptr->map_types[i] = MAP_IDENTITY;

        for (j = 0; j < MAX_PARAMETERS; j++)
        {
            ifsptr->param_sets[i*MAX_PARAMETERS+j] = 0.0;
        }

        ifsptr->probs[i] = 256;

        ifsptr->color_sets[i*3  ] = 0;
        ifsptr->color_sets[i*3+1] = 0;
        ifsptr->color_sets[i*3+2] = 0;
    }

    for (i = 0; i < ifsptr->cache_size; i++)
    {
        ifsptr->sequence[i] = 2;
        ifsptr->x_cache[i] = 0.0;
        ifsptr->y_cache[i] = 0.0;
    }
/*
    ifsptr->param_sets[ 3] = 0.5;
    ifsptr->param_sets[ 5] = 0.0;
    ifsptr->param_sets[10] = 0.5;
    ifsptr->param_sets[11] = (float)sqrt(3.0) * 0.25f;
    ifsptr->color_sets[0] = (1<<BITS_PER_CHANNEL)-1;

    ifsptr->param_sets[15] = 0.5;
    ifsptr->param_sets[17] = -0.5;
    ifsptr->param_sets[22] = 0.5;
    ifsptr->param_sets[23] = -(float)sqrt(3.0) * 0.25f;
    ifsptr->color_sets[4] = (1<<BITS_PER_CHANNEL)-1;

    ifsptr->param_sets[27] = 0.5;
    ifsptr->param_sets[29] = 0.5;
    ifsptr->param_sets[34] = 0.5;
    ifsptr->param_sets[35] = -(float)sqrt(3.0) * 0.25f;
    ifsptr->color_sets[8] = (1<<BITS_PER_CHANNEL)-1;
*/
    ifsptr->param_sets[ 3] = 0.5;
    ifsptr->param_sets[ 5] = 0.0;
    ifsptr->param_sets[10] = 0.5;
    ifsptr->param_sets[11] = (float)sqrt(3.0) * 0.25f;
    ifsptr->param_sets[15] = 1.0;
    ifsptr->param_sets[22] = 1.0;
    ifsptr->color_sets[0] = (1<<BITS_PER_CHANNEL)-1;

    ifsptr->param_sets[27] = 0.5;
    ifsptr->param_sets[29] = -0.5;
    ifsptr->param_sets[34] = 0.5;
    ifsptr->param_sets[35] = -(float)sqrt(3.0) * 0.25f;
    ifsptr->param_sets[39] = 1.0;
    ifsptr->param_sets[46] = 1.0;
    ifsptr->color_sets[4] = (1<<BITS_PER_CHANNEL)-1;

    ifsptr->param_sets[51] = 0.5;
    ifsptr->param_sets[53] = 0.5;
    ifsptr->param_sets[58] = 0.5;
    ifsptr->param_sets[59] = -(float)sqrt(3.0) * 0.25f;
    ifsptr->param_sets[63] = 1.0;
    ifsptr->param_sets[70] = 1.0;
    ifsptr->color_sets[8] = (1<<BITS_PER_CHANNEL)-1;

    memset(ifsptr->sequence, 0, sizeof(int)   * ifsptr->cache_size);
    memset(ifsptr->x_cache,  0, sizeof(float) * ifsptr->cache_size);
    memset(ifsptr->y_cache,  0, sizeof(float) * ifsptr->cache_size);

    ifsptr->x_last = 0.0f;
    ifsptr->y_last = 0.0f;

    ifsptr->buffer.data = NULL;
    ifsptr->buffer.size = 0;
}

void gifs_cleanup(struct gifs_t* ifsptr)
{
    free(ifsptr->function_types);
    free(ifsptr->map_types);
    free(ifsptr->param_sets);
    free(ifsptr->probs);
    free(ifsptr->color_sets);

    gifs_safe_free(ifsptr->sequence);
    gifs_safe_free(ifsptr->x_cache);
    gifs_safe_free(ifsptr->y_cache);
    gifs_safe_free(ifsptr->buffer.data);
}

int gifs_perturb(struct gifs_t* ifsptr)
{
    int j, k;
    int change_type;
    int selection_a;
    int selection_b;
    int max_channel_value = 1 << BITS_PER_CHANNEL;

    struct gifs_t temp_ifs;
    gifs_init(&temp_ifs, 1);
    gifs_imitate(&temp_ifs, ifsptr);

pick_again:

    change_type = SELECT(8);

    switch (change_type)
    {
        case  0:
        // Alter each value by a small amount
        for (j = 0; j < ifsptr->num_functions; j++)
        {
            for (k = 0; k < MAX_PARAMETERS; k++)
            {
                if (0.0f == ifsptr->param_sets[j*MAX_PARAMETERS + k])
                {
                    ifsptr->param_sets[j*MAX_PARAMETERS + k] = gifs_gen_new_param(k);
                }
                else
                {
                    ifsptr->param_sets[j*MAX_PARAMETERS + k] *= 1.5f * RANDFLOAT() + 0.5f;
                }
            }
        }
        break;

        case  1:
        // Change the number of active functions
        if ((B_RAND() < (B_RAND_MAX>>1) && 0 != ifsptr->num_functions) || MAX_FUNCTIONS == ifsptr->num_functions)
        {
            selection_a = SELECT(ifsptr->num_functions);
            if (selection_a == ifsptr->num_functions)
            {
                ifsptr->num_functions--;
            }
            else
            {
                ifsptr->num_functions--;
                ifsptr->function_types[selection_a] = ifsptr->function_types[ifsptr->num_functions];
                memcpy(ifsptr->param_sets + selection_a*MAX_PARAMETERS,
                       ifsptr->param_sets + ifsptr->num_functions*MAX_PARAMETERS, MAX_PARAMETERS*sizeof(float));
                memcpy(ifsptr->color_sets + selection_a*3, ifsptr->color_sets + ifsptr->num_functions*3, 3*sizeof(int));
            }
        }
        else
        {
            for (k = 0; k < MAX_PARAMETERS; k++)
            {
                ifsptr->param_sets[MAX_PARAMETERS*ifsptr->num_functions + k] = gifs_gen_new_param(k);
            }

            for (k = 0; k < 3; k++)
            {
                ifsptr->color_sets[3*ifsptr->num_functions + k] = SELECT(max_channel_value);
            }
            ifsptr->num_functions++;
        }
        break;

        case  2:
        // Swap the sign of some parameter
        selection_a = SELECT(ifsptr->num_functions);
        selection_b = SELECT(MAX_PARAMETERS);
        ifsptr->param_sets[selection_a*MAX_PARAMETERS+selection_b] = 
            -ifsptr->param_sets[selection_a*MAX_PARAMETERS+selection_b];
        break;

        case  3:
        // Change the color of some attractor
        selection_a = SELECT(ifsptr->num_functions);
        ifsptr->color_sets[selection_a*3  ] = SELECT(max_channel_value);
        ifsptr->color_sets[selection_a*3+1] = SELECT(max_channel_value);
        ifsptr->color_sets[selection_a*3+2] = SELECT(max_channel_value);
        break;

        case  4:
        // Regenerate one of the functions
        selection_a = SELECT(ifsptr->num_functions);
        for (k = 0; k < MAX_PARAMETERS; k++)
        {
            ifsptr->param_sets[selection_a*MAX_PARAMETERS + k] = gifs_gen_new_param(k);
        }
        break;

        case  5:
        // Change the map type of one of the functions
        selection_a = SELECT(ifsptr->num_functions);
        selection_b = SELECT(MAP_NUMTYPES*2-1);
        selection_b = (selection_b >= MAP_NUMTYPES) ? MAP_IDENTITY : selection_b;
        if (selection_b == ifsptr->map_types[selection_a])
        {
            goto pick_again;
        }
        else
        {
            ifsptr->map_types[selection_a] = selection_b;
        }
        break;

        case  6:
        // Either zero out or regenerate some of the quadratic parameters of one of the functions
        selection_a = SELECT(ifsptr->num_functions);
        selection_b = selection_a*MAX_PARAMETERS + SELECT(4) * 6;
        if (0.0 == ifsptr->param_sets[selection_b] && 0.0 == ifsptr->param_sets[selection_b+1] &&
            0.0 == ifsptr->param_sets[selection_b+2])
        {
            ifsptr->param_sets[selection_b] = gifs_gen_new_param(selection_b);
            ifsptr->param_sets[selection_b+1] = gifs_gen_new_param(selection_b+1);
            ifsptr->param_sets[selection_b+2] = gifs_gen_new_param(selection_b+2);
        }
        else
        {
            ifsptr->param_sets[selection_b] = 0.0;
            ifsptr->param_sets[selection_b+1] = 0.0;
            ifsptr->param_sets[selection_b+2] = 0.0;
        }
        break;

        case  7:
        // Either reset one of the function stages to identity or regenerate it
        selection_a = SELECT(ifsptr->num_functions);
        selection_b = selection_a*MAX_PARAMETERS + SELECT(4) * 6;
        if (0.0 == ifsptr->param_sets[selection_b]   && 0.0 == ifsptr->param_sets[selection_b+1] &&
            0.0 == ifsptr->param_sets[selection_b+2] &&
            ((0 == selection_a%2 && 1.0 == ifsptr->param_sets[selection_b+3] && 0.0 == ifsptr->param_sets[selection_b+4]) ||
             (1 == selection_a%2 && 0.0 == ifsptr->param_sets[selection_b+3] && 1.0 == ifsptr->param_sets[selection_b+4])) &&
            0.0 == ifsptr->param_sets[selection_b+5])
        {
            ifsptr->param_sets[selection_b] = gifs_gen_new_param(selection_b);
            ifsptr->param_sets[selection_b+1] = gifs_gen_new_param(selection_b+1);
            ifsptr->param_sets[selection_b+2] = gifs_gen_new_param(selection_b+2);
            ifsptr->param_sets[selection_b+3] = gifs_gen_new_param(selection_b+3);
            ifsptr->param_sets[selection_b+4] = gifs_gen_new_param(selection_b+4);
            ifsptr->param_sets[selection_b+5] = gifs_gen_new_param(selection_b+5);
        }
        else
        {
            ifsptr->param_sets[selection_b] = 0.0;
            ifsptr->param_sets[selection_b+1] = 0.0;
            ifsptr->param_sets[selection_b+2] = 0.0;
            if (0 == selection_a%2)
            {
                ifsptr->param_sets[selection_b+3] = 0.0;
                ifsptr->param_sets[selection_b+4] = 1.0;
            }
            else
            {
                ifsptr->param_sets[selection_b+3] = 1.0;
                ifsptr->param_sets[selection_b+4] = 0.0;
            }
            ifsptr->param_sets[selection_b+5] = 0.0;
        }
        break;

        default:
        break;
    }

    if (10000 > gifs_compute(ifsptr, 10000))
    {
        gifs_imitate(ifsptr, &temp_ifs);
        goto pick_again;
    }

    gifs_cleanup(&temp_ifs);

    return 1;
}


#define COMPUTE_MAP(mapid, xlast, ylast, tempx, tempy, oneOverR, r, theta, cosine, sine, tempa, tempb, tempc) \
    switch (mapid) \
    { \
        case MAP_SINE: \
        tempx = (float)sin(xlast); \
        tempy = (float)sin(ylast); \
        break; \
 \
        case MAP_SPHERICAL: \
        r = xlast * xlast + ylast * ylast; \
        if (r > 0.0f) \
        { \
            oneOverR = 1.0f / r; \
            tempx = xlast * oneOverR; \
            tempy = ylast * oneOverR; \
        } \
        else \
        { \
            tempx = 1.0e30f; \
            tempy = 1.0e30f; \
        } \
        break; \
 \
        case MAP_SWIRL: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast) + r; \
        tempx = r * (float)cos(theta); \
        tempy = r * (float)sin(theta); \
        break; \
 \
        case MAP_HORSESHOE: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        if (r > 0.0f) \
        { \
            cosine = xlast / r; \
            sine   = ylast / r; \
            tempx  = r * (cosine * cosine - 0.5f); \
            tempx += tempx; \
            tempy  = r * sine * cosine; \
            tempy += tempy; \
        } \
        else \
        { \
            tempx = 0.0f; \
            tempy = 0.0f; \
        } \
        break; \
 \
        case MAP_POLAR: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast); \
        tempx = theta * oneOverPi; \
        tempy = r - 1.0f; \
        break; \
 \
        case MAP_HANDKERCHIEF: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast); \
        tempx = r * (float)cos(theta + r); \
        tempy = r * (float)sin(theta - r); \
        break; \
 \
        case MAP_HEART: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast); \
        tempx =  r * (float)sin(theta * r); \
        tempy = -r * (float)cos(theta * r); \
        break; \
 \
        case MAP_DISC: \
        r = (float)sqrt(xlast * xlast + ylast * ylast) * localPi;   \
        theta = (float)atan2(ylast, xlast) * oneOverPi; \
        tempx = theta * (float)sin(r); \
        tempy = theta * (float)cos(r); \
        break; \
 \
        case MAP_SPIRAL: \
        r = (float)sqrt(xlast * xlast + ylast * ylast) * localPi; \
        if (0.0f != r) \
        { \
            theta = (float)atan2(ylast, xlast) * oneOverPi; \
            oneOverR = 1.0f / r; \
            tempx = oneOverR * ((float)cos(theta) + (float)sin(r)); \
            tempy = oneOverR * ((float)sin(theta) - (float)cos(r)); \
        } \
        else \
        { \
            tempx = 1.0e+30f; \
            tempy = 1.0e+30f; \
        } \
        break; \
 \
        case MAP_HYPERBOLIC: \
        r = (float)sqrt(xlast * xlast + ylast * ylast) * localPi; \
        if (0.0f != r) \
        { \
            theta = (float)atan2(ylast, xlast) * oneOverPi; \
            oneOverR = 1.0f / r; \
            tempx = oneOverR * (float)sin(theta); \
            tempy = r * (float)cos(theta); \
        } \
        else \
        { \
            tempx = 1.0e+30f; \
            tempy = 0.0f; \
        } \
        break; \
 \
        case MAP_DIAMOND: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast); \
        tempx = (float)sin(theta) * (float)cos(r); \
        tempy = (float)cos(theta) * (float)sin(r); \
        break; \
 \
        case MAP_EX: \
        r = (float)sqrt(xlast * xlast + ylast * ylast); \
        theta = (float)atan2(ylast, xlast); \
        cosine = (float)cos(theta - r); \
        sine   = (float)sin(theta + r); \
        tempx  = r * sine * sine * sine; \
        tempy  = r * cosine * cosine * cosine; \
        break; \
 \
        case MAP_JULIA: \
        r = (float)sqrt(sqrt(xlast * xlast + ylast * ylast)); \
        theta = (float)atan2(ylast, xlast) * 0.5f + SELECT(2) * localPi; \
        tempx = r * (float)cos(theta); \
        tempy = r * (float)sin(theta); \
        break; \
 \
        case MAP_BENT: \
        tempx = xlast; \
        tempy = ylast; \
        if (xlast < 0.0f) \
            tempx *= 2.0f; \
        if (ylast < 0.0f) \
            tempy *= 0.5f; \
        break; \
 \
        case MAP_WAVES: \
        if (paramptr[5] != 0) \
            tempx = xlast + paramptr[4]  * (float)sin(ylast/(paramptr[5]*paramptr[5])); \
        else \
            tempx = xlast + paramptr[4]  * (float)sin(ylast); \
        if (paramptr[11] != 0) \
            tempy = ylast + paramptr[10] * (float)sin(xlast/(paramptr[11]*paramptr[11])); \
        else \
            tempy = ylast + paramptr[10] * (float)sin(xlast); \
        break; \
 \
        case MAP_FISHEYE: \
        oneOverR = 2.0f / (1.0f + (float)sqrt(xlast*xlast + ylast*ylast)); \
        tempx = ylast * oneOverR; \
        tempy = ylast * oneOverR; \
        break; \
 \
        case MAP_POPCORN: \
        tempx = xlast + paramptr[5]  * (float)sin(tan(3.0f * ylast)); \
        tempy = ylast + paramptr[11] * (float)sin(tan(3.0f * ylast)); \
        break; \
 \
        case MAP_EXPONENTIAL: \
        r = (float)exp(xlast - 1.0f); \
        tempx = r * (float)cos(localPi * ylast); \
        tempy = r * (float)sin(localPi * ylast); \
        break; \
\
        case MAP_POWER: \
        theta = (float)atan2(ylast, xlast); \
        tempy = (float)sin(theta); \
        r = (float)sqrt(xlast*xlast + ylast*ylast); \
        r = (float)pow(r, tempy); \
        tempx = (float)cos(theta) * r; \
        tempy *= r; \
        break; \
 \
        case MAP_COSINE: \
        tempx = (float)(cos(localPi * xlast) * cosh(ylast)); \
        tempy = -(float)(sin(localPi * xlast) * sinh(ylast)); \
        break; \
 \
        case MAP_RINGS: \
        tempa = paramptr[5] * paramptr[5]; \
        if (tempa != 0.0f) \
        { \
            r = (float)sqrt(xlast*xlast + ylast*ylast); \
            tempb = (float)fmod(r + tempa, tempa*2.0f) - tempa + r * (1.0f - tempa); \
            tempx = xlast * tempb; \
            tempy = ylast * tempb; \
        } \
        else \
        { \
            tempx = 1.0e30f; \
            tempy = 1.0e30f; \
        } \
        break; \
 \
        case MAP_FAN: \
        tempa = localPi * paramptr[5] * paramptr[5]; \
        tempb = tempa * 0.5f; \
        theta = (float)atan2(ylast, xlast); \
        if (tempa != 0.0f) \
        { \
            tempc = (float)fmod(theta + paramptr[11], tempa); \
            r = (float)sqrt(xlast*xlast + ylast*ylast); \
            if (tempc > tempb) \
                tempb = -tempb; \
            tempx = r * (float)cos(theta + tempb); \
            tempy = r * (float)sin(theta + tempb); \
        } \
        else \
        { \
            tempx = r * (float)cos(theta + tempb); \
            tempy = r * (float)sin(theta + tempb); \
        } \
        break; \
 \
        case MAP_BLOB: \
        theta = (float)atan2(ylast, xlast); \
        tempa = 1.4f + 0.7f * ((float)sin(5.0f * theta) + 1.0f); \
        tempx = tempa * xlast; \
        tempy = tempa * ylast; \
        break; \
 \
        case MAP_EYEFISH: \
        r = (float)sqrt(xlast*xlast + ylast*ylast); \
        tempa = 2.0f / (r + 1.0f); \
        tempx = tempa * xlast; \
        tempy = tempa * ylast; \
        break; \
 \
        case MAP_BUBBLE: \
        tempa = 4.0f / (xlast*xlast + ylast*ylast + 4.0f); \
        tempx = tempa * xlast; \
        tempy = tempa * ylast; \
        break; \
 \
        case MAP_CYLINDER: \
        tempx = (float)sin(xlast); \
        tempy = ylast; \
        break; \
 \
        case MAP_CURL: \
        tempa = 1 + 1.1f * (xlast + xlast*xlast - ylast*ylast); \
        tempb = 1.1f * ylast * (1.0f + 2.0f * xlast); \
        tempc = 1.0f / (tempa*tempa + tempb*tempb); \
        tempx = tempc * (tempa*xlast + tempb*ylast); \
        tempy = tempc * (-tempb*xlast + tempa*ylast); \
        break; \
 \
        case MAP_RECTANGLES: \
        tempx = (2.0f * (float)floor(xlast/2.7f) + 1.0f) * 2.7f - xlast; \
        tempy = (2.0f * (float)floor(ylast/2.7f) + 1.0f) * 2.7f - ylast; \
        break; \
 \
        case MAP_TANGENT: \
        tempx = (float)sin(xlast) / ((float)cos(ylast) + 1.0e-12f); \
        if (0.0f != (float)cos(xlast)) { \
            tempy = (float)tan(xlast); \
        } \
        else \
        { \
            tempy = 1.0e30f; \
        } \
        break; \
 \
        case MAP_SECANT: \
        tempx = xlast; \
        tempy = 1.0f / (1.4f * (float)cos(1.4f * (float)sqrt(xlast*xlast + ylast*ylast))); \
        break; \
 \
        default: \
        break; \
    }


int gifs_compute(struct gifs_t* ifsptr, int count)
{
    float xlast = ifsptr->x_last;
    float ylast = ifsptr->y_last;

    const float localPi   = 3.1415927f;
    const float oneOverPi = 1.0f / localPi;

    int     i = 0;
    int     j;
    float*  x = ifsptr->x_cache;
    float*  y = ifsptr->y_cache;
    float   r = 0.0f;
    float   oneOverR = 1.0e30f;
    float   theta = 0.0f;
    float   cosine, sine;
    int     fxn;
    int     fxnid;
    int     mapid;
    float*  paramptr;
    int*    seqptr = ifsptr->sequence;
    float   tempa;
    float   tempb;
    float   tempc;
    float   tempx;
    float   tempy;
    int     bounded = 1;
    int     thingy = ((int)time(NULL) & 0xf) * 17;

    for (j = 0; j < thingy; j++)
    {
        B_RAND();
    }

    if ((count <= 0) || (count > ifsptr->cache_size))
    {
        count = ifsptr->cache_size;
    }

try_again:

    for (j = 0; j < 100; j++)
    {
        fxn = SELECT(ifsptr->num_functions);
        fxnid = ifsptr->function_types[fxn];
        paramptr = ifsptr->param_sets + fxn*MAX_PARAMETERS;

        tempx = (paramptr[ 0] * xlast + paramptr[ 1] * ylast + paramptr[ 3]) * xlast +
                                       (paramptr[ 2] * ylast + paramptr[ 4]) * ylast + paramptr[ 5];

        tempy = (paramptr[ 6] * xlast + paramptr[ 7] * ylast + paramptr[ 9]) * xlast +
                                       (paramptr[ 8] * ylast + paramptr[10]) * ylast + paramptr[11];

        // xlast and ylast are used as temp variables here.
        xlast = tempx;
        ylast = tempy;

        mapid = ifsptr->map_types[fxn];
        COMPUTE_MAP(mapid, xlast, ylast, tempx, tempy, oneOverR, r, theta, cosine, sine, tempa, tempb, tempc)
        
        xlast = tempx;
        ylast = tempy;

        tempx = (paramptr[12] * xlast + paramptr[13] * ylast + paramptr[15]) * xlast +
                                       (paramptr[14] * ylast + paramptr[16]) * ylast + paramptr[17];
        tempy = (paramptr[18] * xlast + paramptr[19] * ylast + paramptr[21]) * xlast +
                                       (paramptr[20] * ylast + paramptr[22]) * ylast + paramptr[23];

        // xlast and ylast are used as temp variables here.
        xlast = tempx;
        ylast = tempy;

        if ((xlast > 1e3) || (xlast < -1e3) || (ylast > 1e3) || (ylast < -1e3))
        {
            xlast = RANDFLOAT();
            ylast = RANDFLOAT();
            bounded = 0;
            break;
        }
    }
    
    if (0 == bounded)
        return i;

    for (; i<count; i++, x++, y++, seqptr++) // initialization happens up top
    {
        fxn = SELECT(ifsptr->num_functions);
        fxnid = ifsptr->function_types[fxn];
        *seqptr = fxn;
        paramptr = ifsptr->param_sets + fxn*MAX_PARAMETERS;

        tempx = (paramptr[ 0] * xlast + paramptr[ 1] * ylast + paramptr[ 3]) * xlast +
                                       (paramptr[ 2] * ylast + paramptr[ 4]) * ylast + paramptr[ 5];
        tempy = (paramptr[ 6] * xlast + paramptr[ 7] * ylast + paramptr[ 9]) * xlast +
                                       (paramptr[ 8] * ylast + paramptr[10]) * ylast + paramptr[11];
        xlast = tempx;
        ylast = tempy;

        mapid = ifsptr->map_types[fxn];
        COMPUTE_MAP(mapid, xlast, ylast, tempx, tempy, oneOverR, r, theta, cosine, sine, tempa, tempb, tempc)

        xlast = tempx;
        ylast = tempy;

        *x = (paramptr[12] * xlast + paramptr[13] * ylast + paramptr[15]) * xlast +
                                    (paramptr[14] * ylast + paramptr[16]) * ylast + paramptr[17];

        *y = (paramptr[18] * xlast + paramptr[19] * ylast + paramptr[21]) * xlast +
                                    (paramptr[20] * ylast + paramptr[22]) * ylast + paramptr[23];

        // xlast and ylast are used as temp variables here.
        xlast = *x;
        ylast = *y;

        if ((xlast > 1e3) || (xlast < -1e3) || (ylast > 1e3) || (ylast < -1e3))
        {
            xlast = (float)B_RAND()/(float)B_RAND_MAX;
            ylast = (float)B_RAND()/(float)B_RAND_MAX;
            bounded = 0;
            break;
        }
    }

    if (0 == bounded)
    {
        bounded = 1;
    	goto try_again;
    }

    ifsptr->x_last = xlast;
    ifsptr->y_last = ylast;

    if (_isnan(ifsptr->x_last) || !_finite(ifsptr->x_last))
    {
        ifsptr->x_last = RANDFLOAT() - 0.5f;
    }

    if (_isnan(ifsptr->y_last) || !_finite(ifsptr->y_last))
    {
        ifsptr->y_last = RANDFLOAT() - 0.5f;
    }

    return i;
}

void gifs_imitate(struct gifs_t* ifsptr, const struct gifs_t* srcptr)
{
    ifsptr->num_functions = srcptr->num_functions;
    memcpy(ifsptr->function_types, srcptr->function_types, sizeof(unsigned) * MAX_FUNCTIONS);
    memcpy(ifsptr->map_types,      srcptr->map_types,      sizeof(unsigned) * MAX_FUNCTIONS);
    memcpy(ifsptr->param_sets, srcptr->param_sets, sizeof(float) * MAX_FUNCTIONS*MAX_PARAMETERS);
    memcpy(ifsptr->color_sets, srcptr->color_sets, sizeof(unsigned) * MAX_FUNCTIONS * 3);
    memset(ifsptr->sequence, 0, sizeof(unsigned) * ifsptr->cache_size);
    memset(ifsptr->x_cache, 0, sizeof(float) * ifsptr->cache_size);
    memset(ifsptr->y_cache, 0, sizeof(float) * ifsptr->cache_size);
}

float gifs_gen_new_param(int paramIndex)
{
    return RANDFLOAT() - 0.5f;
}

void compute_buffer_size(struct gifs_t* ifsptr)
{
    if (0 == ifsptr->buffer.size)
    {
        unsigned buffer_size = 0;

        buffer_size += sizeof(int);                                            // ifsptr->num_functions
        buffer_size += sizeof(int) *   ifsptr->num_functions;                  // ifsptr->function_types
        buffer_size += sizeof(int) *   ifsptr->num_functions;                  // ifsptr->map_types
        buffer_size += sizeof(float) * ifsptr->num_functions * MAX_PARAMETERS; // ifsptr->param_sets
        buffer_size += sizeof(int) *   ifsptr->num_functions * 3;              // ifsptr->color_sets

        ifsptr->buffer.size = buffer_size;
    }
}


struct gifs_buffer_t* gifs_export_buffer(struct gifs_t* ifsptr)
{
    char* buffer_pointer;
    unsigned copySize;

    if (NULL == ifsptr->buffer.data)
    {
        ifsptr->buffer.data = (char*)malloc(ifsptr->buffer.size);
    }

    buffer_pointer = ifsptr->buffer.data;
    copySize = sizeof(int);
    memcpy(buffer_pointer, &ifsptr->num_functions, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions;
    memcpy(buffer_pointer, ifsptr->function_types, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions;
    memcpy(buffer_pointer, ifsptr->map_types, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(float) * ifsptr->num_functions * MAX_PARAMETERS;
    memcpy(buffer_pointer, ifsptr->param_sets, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions * 3;
    memcpy(buffer_pointer, ifsptr->color_sets, copySize);

    return &ifsptr->buffer;
}


int gifs_import_buffer(struct gifs_t* ifsptr, void* buffer)
{
    char* buffer_pointer;
    unsigned copySize;

    if (NULL == buffer)
    {
        return 0;
    }

    buffer_pointer = (char*)buffer;
    copySize = sizeof(int);
    memcpy(&ifsptr->num_functions, buffer_pointer, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions;
    memcpy(ifsptr->function_types, buffer_pointer, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions;
    memcpy(ifsptr->map_types, buffer_pointer, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(float) * ifsptr->num_functions * MAX_PARAMETERS;
    memcpy(ifsptr->param_sets, buffer_pointer, copySize);

    buffer_pointer += copySize;
    copySize = sizeof(int) * ifsptr->num_functions * 3;
    memcpy(ifsptr->color_sets, buffer_pointer, copySize);

    return 1;
}
