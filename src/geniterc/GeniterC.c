#include  <stdio.h>
#include <stdlib.h>
#include   <math.h>
#include   <time.h>

#define _WIN32_WINNT 0x0500 

#include         <windows.h>
#include             "sdl.h"
#include "GeneralizedIFSC.h"

#pragma warning(disable: 4996)

#define ENABLE_VIDEO

#define  MAX_COMPUTE MAX_CACHE_SIZE
#define BASE_COMPUTE 32768
#define   IMAGE_SIZE   768
#define TOTAL_PIXELS IMAGE_SIZE*IMAGE_SIZE


enum
{
    TWEAKMODE_NONE,
    TWEAKMODE_SCALE,
    TWEAKMODE_XY,
    TWEAKMODE_GAMMA,
    TWEAKMODE_VIBRANCE,
    TWEAKMODE_CONTRAST,
    TWEAKMODE_NUMTYPES
};

enum
{
    RENDERMODE_OFF,
    RENDERMODE_ONESHOT,
    RENDERMODE_CONTINUOUS,
    RENDERMODE_NUMTYPES
};


struct engine_state_t
{
    struct gifs_t* ifss;
    struct gifs_t* curr_ifs;
    unsigned sequence_length;
    unsigned num_gifss;
    unsigned num_rows;
    unsigned num_columns;
    unsigned num_pixels;
    float scale;
    float translation_x;
    float translation_y;
    float gamma;
    float vibrance;
    float contrast;
    unsigned* r_plane;
    unsigned* g_plane;
    unsigned* b_plane;
    float*    n_plane;
    HANDLE console;
    unsigned tweak_mode;
    int show_calc_rate;
    BOOL calc_rate_avail;
    float num_kpts;
    float num_ksecs;
    int aa_on;
    unsigned render_on;
    int pts_to_next_render;
    int user_req_exit;
};


void collect_pts(struct engine_state_t* state_ptr);
void paint_pts(struct engine_state_t* state_ptr);
void handle_event(SDL_Event* sdl_event_ptr, struct engine_state_t* state_ptr);
void handle_kbd_event(SDLKey key, struct engine_state_t* state_ptr);
void handle_tweak(SDLKey key, struct engine_state_t* state_ptr);
void change_curr_ifs(struct engine_state_t* state_ptr, unsigned new_index);
void reset_gfx(struct engine_state_t* state_ptr);
void print_help(struct engine_state_t* state_ptr);
void print_params(struct engine_state_t* state_ptr);
void print_details(struct engine_state_t* state_ptr);
int save_params(struct engine_state_t* state_ptr, char* file_name);
int load_params(struct engine_state_t* state_ptr, char* file_name);
float zoom_level(struct engine_state_t* state_ptr);

unsigned g_workspc[IMAGE_SIZE * IMAGE_SIZE];


int main(int argc, char* argv[])
{
    struct engine_state_t engine_state;
    size_t gifs_size = 0;
    size_t malloc_size = 0;
    LARGE_INTEGER counter_freq, counter_a, counter_b;
    unsigned i = 0;
    int result = 0xdeadbeef;
    SDL_Event sdl_event;
    unsigned clock_start = 0;
    int event_occurred = 0;
    double totalSeconds = 0.0;
    double totalPoints  = 0.0;

    memset(&engine_state, 0, sizeof(struct engine_state_t));

    engine_state.num_gifss     = 10;
    engine_state.num_rows      = IMAGE_SIZE;
    engine_state.num_columns   = IMAGE_SIZE;
    engine_state.num_pixels    = engine_state.num_rows * engine_state.num_columns;
    engine_state.scale         = IMAGE_SIZE * 0.5f;
    engine_state.translation_x = 0.0f;
    engine_state.translation_y = 0.0f;
    engine_state.gamma         = 4.0f;
    engine_state.vibrance      = 1.0f;
    engine_state.contrast      = 2.0f;
    
    gifs_size = sizeof(struct gifs_t);
    engine_state.ifss = malloc(sizeof(struct gifs_t) * engine_state.num_gifss);

    malloc_size = sizeof(unsigned) * engine_state.num_pixels;
    engine_state.r_plane = (unsigned*)malloc(malloc_size);
    memset(engine_state.r_plane, 0, malloc_size);

    malloc_size = sizeof(unsigned) * engine_state.num_pixels;
    engine_state.g_plane = (unsigned*)malloc(malloc_size);
    memset(engine_state.g_plane, 0, malloc_size);

    malloc_size = sizeof(unsigned) * engine_state.num_pixels;
    engine_state.b_plane = (unsigned*)malloc(malloc_size);
    memset(engine_state.b_plane, 0, malloc_size);

    malloc_size = sizeof(float) * engine_state.num_pixels;
    engine_state.n_plane = (float*)malloc(malloc_size);
    memset(engine_state.n_plane, 0, malloc_size);

    engine_state.show_calc_rate = 0;
    engine_state.calc_rate_avail = QueryPerformanceFrequency(&counter_freq);
    engine_state.user_req_exit = 0;

    engine_state.render_on = RENDERMODE_CONTINUOUS;
    engine_state.pts_to_next_render = MAX_COMPUTE;

    for (i = 0; i < engine_state.num_gifss; i++)
    {
        gifs_init(&engine_state.ifss[i], 0);
    }

    engine_state.curr_ifs = &engine_state.ifss[1];

    if (load_params(&engine_state, "lastsave.dat"))
    {
        gifs_imitate(&engine_state.ifss[0], engine_state.curr_ifs);
        for (i = 2; i < engine_state.num_gifss; i++)
        {
            gifs_imitate(&engine_state.ifss[i], engine_state.curr_ifs);
        }
    }

    gifs_perturb(&engine_state.ifss[0]);
    for (i = 2; i < engine_state.num_gifss; i++)
    {
        gifs_perturb(&engine_state.ifss[i]);
    }

    AllocConsole();
    engine_state.console = GetStdHandle(STD_OUTPUT_HANDLE);
    freopen("CON", "w", stdout);
    printf("\nWelcome to Geniter. For help, click on the graphics window and hit 'h'. This\n");
    printf("text window is resizeable.\n");

    printf("\nTo quit, hit 'q'. The X button in the title bar doesn't work yet!\n");

    result = SDL_Init(SDL_INIT_VIDEO);
    SDL_WM_SetCaption("Geniter", NULL);
    SDL_SetVideoMode(IMAGE_SIZE, IMAGE_SIZE, 32, SDL_SWSURFACE);

    if (engine_state.calc_rate_avail)
    {
        QueryPerformanceCounter(&counter_a);
    }

    do {
        if (0 != SDL_PollEvent(&sdl_event) && SDL_KEYDOWN == sdl_event.type)
        {
            event_occurred = 1;
            handle_event(&sdl_event, &engine_state);
            clock_start = SDL_GetTicks();
        }

        if (0 != event_occurred)
        {
            if ((SDL_GetTicks() - clock_start) > 500)
            {
                event_occurred = 0;
            }
            else
            {
                Sleep(50);
            }
        }

        if (0 != event_occurred)
        {
            if (engine_state.calc_rate_avail)
            {
                QueryPerformanceCounter(&counter_b);
                counter_a = counter_b;
            }
        }
        else
        {
            engine_state.sequence_length = gifs_compute(engine_state.curr_ifs, BASE_COMPUTE);
            collect_pts(&engine_state);
            engine_state.pts_to_next_render -= BASE_COMPUTE;

            if (RENDERMODE_ONESHOT == engine_state.render_on ||
                    (RENDERMODE_CONTINUOUS == engine_state.render_on &&
                     0 >= engine_state.pts_to_next_render))
            {
                paint_pts(&engine_state);

                if (0 != engine_state.calc_rate_avail)
                {
                    float numSeconds = 0;
                    QueryPerformanceCounter(&counter_b);
                    numSeconds = (float)(counter_b.QuadPart - counter_a.QuadPart) / (float)counter_freq.QuadPart;
                    engine_state.num_kpts += BASE_COMPUTE * 0.001f;
                    engine_state.num_ksecs += numSeconds * 0.001f;
                    counter_a = counter_b;
                }

                if (RENDERMODE_ONESHOT == engine_state.render_on)
                {
                    engine_state.render_on = RENDERMODE_OFF;
                }
                else
                {
                    engine_state.pts_to_next_render += MAX_COMPUTE;
                }
            }
        }
    }
    while (0 == engine_state.user_req_exit);

    SDL_Quit();
    FreeConsole();

    return 0;
}


void collect_pts(struct engine_state_t* state_ptr)
{
    unsigned i;
    unsigned plane_index;
    //unsigned neighborIndex;
    float x, y;
    int x_truncated, y_truncated;
    unsigned curr_r = 0;
    unsigned curr_g = 0;
    unsigned curr_b = 0;
    float x_skew = 0.0f;
    float y_skew = 0.0f;
    float x_skew_absolute = 0.0f;
    float y_skew_absolute = 0.0f;
    float radius = 0.4f;
    float spill_area = 0.0f;
    float x_chord_measure = 0.0f;
    float y_chord_measure = 0.0f;
    float corner_measure = 0.0f;
    float translate_x =  state_ptr->translation_x + state_ptr->num_rows / state_ptr->scale * 0.5f;
    float translate_y = -state_ptr->translation_y + state_ptr->num_rows / state_ptr->scale * 0.5f;
/*
    if (state_ptr->aa_on)
    {
      for (i = 0; i < state_ptr->sequence_length; i++)
        {
            curr_r += state_ptr->curr_ifs->state_ptr[state_ptr->curr_ifs->sequence[i]*3];
            curr_r >>= 1;
            curr_g += state_ptr->curr_ifs->state_ptr[state_ptr->curr_ifs->sequence[i]*3+1];
            curr_g >>= 1;
            curr_b += state_ptr->curr_ifs->state_ptr[state_ptr->curr_ifs->sequence[i]*3+2];
            curr_b >>= 1;

            x = (translate_x + state_ptr->curr_ifs->x_cache[i]) * state_ptr->scale;
            x_truncated = (int)(x + 1) - 1;

            if ((x_truncated < 1) || (x_truncated >= IMAGE_SIZE-1))
            {
                continue;
            }

            y = (translate_y - state_ptr->curr_ifs->y_cache[i]) * state_ptr->scale;
            y_truncated = (int)(y + 1) - 1;

            if ((y_truncated < 1) || (y_truncated >= IMAGE_SIZE-1))
            {
                continue;
            }

            plane_index = x_truncated + y_truncated * IMAGE_SIZE;

            s_r_plane[plane_index] += curr_r;
            s_g_plane[plane_index] += curr_g;
            s_b_plane[plane_index] += curr_b;

            state_ptr->n_plane[plane_index]++;

            x_chord_measure = 0.0f;
            y_chord_measure = 0.0f;
            corner_measure = 0.0f;
    
            x_skew = x - x_truncated - 0.5f;
            x_skew_absolute = fabs(x_skew);
            x_chord_measure = x_skew_absolute + radius - 0.5f;
            //x_chord_measure *= x_chord_measure;
    
            y_skew = y - y_truncated - 0.5f;
            y_skew_absolute = fabs(y_skew);
            y_chord_measure = y_skew_absolute + radius - 0.5f;
            //y_chord_measure *= y_chord_measure;
    
            corner_measure = x_chord_measure * y_chord_measure;
            x_chord_measure -= corner_measure;
            y_chord_measure -= corner_measure;
    
            spill_area = x_chord_measure + y_chord_measure + corner_measure;
    
            s_r_plane[plane_index] -= (unsigned)(spill_area * curr_r);
            s_g_plane[plane_index] -= (unsigned)(spill_area * curr_g);
            s_b_plane[plane_index] -= (unsigned)(spill_area * curr_b);
    
            state_ptr->n_plane[plane_index] -= spill_area;
    
            if (0.0f < x_chord_measure)
            {
                if (x_skew > 0.0f)
                {
                    neighborIndex = plane_index + 1;
                    s_r_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_r);
                    s_g_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_g);
                    s_b_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_b);
                    state_ptr->n_plane[neighborIndex] += x_chord_measure;
                }
                else
                {
                    neighborIndex = plane_index - 1;
                    s_r_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_r);
                    s_g_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_g);
                    s_b_plane[neighborIndex] += (unsigned)(x_chord_measure * curr_b);
                    state_ptr->n_plane[neighborIndex] += x_chord_measure;
                }
            }
    
            if (0.0f < y_chord_measure)
            {
                if (y_skew > 0.0f)
                {
                    neighborIndex = plane_index + IMAGE_SIZE;
                    s_r_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_r);
                    s_g_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_g);
                    s_b_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_b);
                    state_ptr->n_plane[neighborIndex] += y_chord_measure;
                }
                else
                {
                    neighborIndex = plane_index - IMAGE_SIZE;
                    s_r_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_r);
                    s_g_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_g);
                    s_b_plane[neighborIndex] += (unsigned)(y_chord_measure * curr_b);
                    state_ptr->n_plane[neighborIndex] += y_chord_measure;
                }
            }
    
            if (0.0f < corner_measure)
            {
                if (x_skew > 0.0f)
                {
                    if (y_skew > 0.0f)
                    {
                        neighborIndex = plane_index + (IMAGE_SIZE+1);
                        s_r_plane[neighborIndex] += (unsigned)(corner_measure * curr_r);
                        s_g_plane[neighborIndex] += (unsigned)(corner_measure * curr_g);
                        s_b_plane[neighborIndex] += (unsigned)(corner_measure * curr_b);
                        state_ptr->n_plane[neighborIndex] += corner_measure;
                    }
                    else
                    {
                        neighborIndex = plane_index - (IMAGE_SIZE-1);
                        s_r_plane[neighborIndex] += (unsigned)(corner_measure * curr_r);
                        s_g_plane[neighborIndex] += (unsigned)(corner_measure * curr_g);
                        s_b_plane[neighborIndex] += (unsigned)(corner_measure * curr_b);
                        state_ptr->n_plane[neighborIndex] += corner_measure;
                    }
                }
                else
                {
                    if (y_skew > 0.0f)
                    {
                        neighborIndex = plane_index + (IMAGE_SIZE-1);
                        s_r_plane[neighborIndex] += (unsigned)(corner_measure * curr_r);
                        s_g_plane[neighborIndex] += (unsigned)(corner_measure * curr_g);
                        s_b_plane[neighborIndex] += (unsigned)(corner_measure * curr_b);
                        state_ptr->n_plane[neighborIndex] += corner_measure;
                    }
                    else
                    {
                        neighborIndex = plane_index - (IMAGE_SIZE+1);
                        s_r_plane[neighborIndex] += (unsigned)(corner_measure * curr_r);
                        s_g_plane[neighborIndex] += (unsigned)(corner_measure * curr_g);
                        s_b_plane[neighborIndex] += (unsigned)(corner_measure * curr_b);
                        state_ptr->n_plane[neighborIndex] += corner_measure;
                    }
                }
            }
        }
    }
    else // antialiasing is off
    { */
        for (i = 0; i < state_ptr->sequence_length; i++)
        {
            curr_r += state_ptr->curr_ifs->color_sets[state_ptr->curr_ifs->sequence[i]*3];
            curr_r >>= 1;
            curr_g += state_ptr->curr_ifs->color_sets[state_ptr->curr_ifs->sequence[i]*3+1];
            curr_g >>= 1;
            curr_b += state_ptr->curr_ifs->color_sets[state_ptr->curr_ifs->sequence[i]*3+2];
            curr_b >>= 1;
    
            x = (translate_x + state_ptr->curr_ifs->x_cache[i]) * state_ptr->scale;
            x_truncated = (int)(x + 1) - 1;
    
            if ((x_truncated < 1) || (x_truncated >= (IMAGE_SIZE-1)))
            {
                continue;
            }
    
            y = (translate_y - state_ptr->curr_ifs->y_cache[i]) * state_ptr->scale;
            y_truncated = (int)(y + 1) - 1;
    
            if ((y_truncated < 1) || (y_truncated >= (IMAGE_SIZE-1)))
            {
                continue;
            }
    
            plane_index = x_truncated + y_truncated * IMAGE_SIZE;
    
            state_ptr->r_plane[plane_index] += curr_r;
            state_ptr->g_plane[plane_index] += curr_g;
            state_ptr->b_plane[plane_index] += curr_b;
    
            state_ptr->n_plane[plane_index]++;
        }
   /* } */
}

static unsigned char s_r_plane[IMAGE_SIZE * IMAGE_SIZE];
static unsigned char s_g_plane[IMAGE_SIZE * IMAGE_SIZE];
static unsigned char s_b_plane[IMAGE_SIZE * IMAGE_SIZE];

void paint_pts(struct engine_state_t* state_ptr)
{
    float max = 0;
    unsigned max_pixel = 0;
    unsigned i;
    float count_correction = 0.0f;
    float corrected_n_plane = 0.0f;
    float level_correction = 0.0f;
    float gamma_recip_minus_one;
    float r_correction = 0.0f;
    float g_correction = 0.0f;
    float b_correction = 0.0f;
    float a_correction = 0.0f;
    float normed_r = 0.0f;
    float normed_g = 0.0f;
    float normed_b = 0.0f;
    float temp = 0.0f;
    float one_minus_vibrance = 1.0f - state_ptr->vibrance;
    SDL_Surface* surface = NULL;
    unsigned* pixels = NULL;
    unsigned int r, g, b;
    float adj_vibrance;
    float adj_contrast;

    for (i = 0; i < state_ptr->num_pixels; i++)
    {
        if (state_ptr->n_plane[i] > max)
        {
            max = state_ptr->n_plane[i];
            max_pixel = i;
        }
    }

    count_correction = 3.0f / (float)max; // 3 has nicer results than 1. Don't know why.
    level_correction = count_correction / 256.0f;
    gamma_recip_minus_one = 1.0f / state_ptr->gamma - 1.0f;
    surface = SDL_GetVideoSurface();
    pixels = (unsigned*)(surface->pixels);

    for (i = 0; i < state_ptr->num_pixels; i++)
    {
        g_workspc[i] = 0;

        if (0.0f == state_ptr->n_plane[i])
        {
            continue;
        }

        //temp = state_ptr->n_plane[i] * count_correction;
        //a_correction = exp(log(temp) * gammaReciprocal) / temp;
        //a_correction = exp(log(temp) * gammaReciprocal - log(temp));
        //a_correction = exp(log(temp) * (gammaReciprocal - 1.0f));
        corrected_n_plane = state_ptr->n_plane[i] * count_correction;
        a_correction = (float)exp(log(corrected_n_plane) * gamma_recip_minus_one);
        adj_vibrance = state_ptr->vibrance * a_correction;
        adj_contrast = state_ptr->contrast * 256.0f;

        if (0 != state_ptr->r_plane[i])
        {
            normed_r = state_ptr->r_plane[i] * level_correction;
            r_correction = (float)exp(log(normed_r) * gamma_recip_minus_one);
            r = (int)((adj_vibrance + one_minus_vibrance * r_correction) *
                       normed_r * adj_contrast);
            r = (r >= 0x80000000) ? 0 : r;
            r = (r >= 255) ? 255 : r;
        }
        else
        {
            r = 0;
        }

        s_r_plane[i] = (char)r;

        if (0 != state_ptr->g_plane[i])
        {
            normed_g = state_ptr->g_plane[i] * level_correction;
            g_correction = (float)exp(log(normed_g) * gamma_recip_minus_one);
            g = (int)((adj_vibrance + one_minus_vibrance * g_correction) *
                       normed_g * adj_contrast);
            g = (g >= 0x80000000) ? 0 : g;
            g = (g >= 255) ? 255 : g;
        }
        else
        {
            g = 0;
        }

        s_g_plane[i] = (char)g;

        if (0 != state_ptr->b_plane[i])
        {
            normed_b = state_ptr->b_plane[i] * level_correction;
            b_correction = (float)exp(log(normed_b) * gamma_recip_minus_one);
            b = (int)((adj_vibrance + one_minus_vibrance * b_correction) *
                       normed_b * adj_contrast);
            b = (b >= 0x80000000) ? 0 : b;
            b = (b >= 255) ? 255 : b;
        }
        else
        {
            b = 0;
        }

        s_b_plane[i] = (char)b;
    }

    if (0 != state_ptr->aa_on)
    {
        const float max_kernel_radius = 4.0f;
        int max_kernel_rad_int = (int)max_kernel_radius;
        int rows, cols;
        for (rows = max_kernel_rad_int; rows < IMAGE_SIZE-max_kernel_rad_int; rows++)
        {
            for (cols = max_kernel_rad_int; cols < IMAGE_SIZE-max_kernel_rad_int; cols++)
            {
                int index = rows*IMAGE_SIZE + cols;
                float kernel_size = max_kernel_radius - (float)log(state_ptr->n_plane[index] + 1.0f);
                int kernel_size_int = 0;
                float norming_scalar = 0.0f;
                float r_sum = 0.0f;
                float g_sum = 0.0f;
                float b_sum = 0.0f;
                float n_sum = 0.0f;
                int kernel_index = 0;
                int dx, dy;
                kernel_size = (kernel_size > 0.0f) ? kernel_size : 0.0f;
                kernel_size_int = (int)kernel_size;
                norming_scalar = 1.0f / ((kernel_size+1.0f) * (kernel_size+1.0f));
                kernel_index = index - kernel_size_int * (IMAGE_SIZE + 1);
                for (dy = -kernel_size_int; dy <= kernel_size_int; dy++)
                {
                    for (dx = -kernel_size_int; dx <= kernel_size_int; dx++)
                    {
                        float kernel_value = (1.0f - (dx*dx + dy*dy) * norming_scalar);
                        kernel_value = (kernel_value > 0.0f) ? kernel_value : 0.0f;
                        r_sum += kernel_value * s_r_plane[kernel_index];
                        g_sum += kernel_value * s_g_plane[kernel_index];
                        b_sum += kernel_value * s_b_plane[kernel_index];
                        n_sum += kernel_value;
                        kernel_index++;
                    }
                    kernel_index += IMAGE_SIZE - (kernel_size_int<<1) - 1;
                }

                n_sum = (n_sum >= 0.0f) ? n_sum : 1.0f;
                n_sum = 1.0f / n_sum;

                r = (int)(r_sum * n_sum);
                r = (r >= 0) ? r : 0;
                r = (r <= 255) ? r : 255;

                g = (int)(g_sum * n_sum);
                g = (g >= 0) ? g : 0;
                g = (g <= 255) ? g : 255;

                b = (int)(b_sum * n_sum);
                b = (b >= 0) ? b : 0;
                b = (b <= 255) ? b : 255;

                pixels[index] = (r << surface->format->Rshift) |
                                (g << surface->format->Gshift) |
                                (b << surface->format->Bshift);
            }
        }
    }
    else
    {
        unsigned i;
        for (i = 0; i < state_ptr->num_pixels; i++)
        {
            pixels[i] = (s_r_plane[i] << surface->format->Rshift) |
                        (s_g_plane[i] << surface->format->Gshift) |
                        (s_b_plane[i] << surface->format->Bshift);
        }
    }

    SDL_UpdateRect(surface, 0, 0, 0, 0);
}


void handle_event(SDL_Event* sdl_event_ptr, struct engine_state_t* state_ptr)
{
    switch (sdl_event_ptr->type)
    {
        case SDL_KEYDOWN:
        handle_kbd_event(sdl_event_ptr->key.keysym.sym, state_ptr);
        break;
    }
}


void handle_kbd_event(SDLKey key, struct engine_state_t* state_ptr)
{
    unsigned i;

    switch (key)
    {
        case SDLK_0:
        change_curr_ifs(state_ptr, 0);
        printf("Displaying image number 10.\n");
        break;

        case SDLK_1:
        change_curr_ifs(state_ptr, 1);
        printf("Displaying image number 1.\n");
        break;

        case SDLK_2:
        change_curr_ifs(state_ptr, 2);
        printf("Displaying image number 2.\n");
        break;

        case SDLK_3:
        change_curr_ifs(state_ptr, 3);
        printf("Displaying image number 3.\n");
        break;

        case SDLK_4:
        change_curr_ifs(state_ptr, 4);
        printf("Displaying image number 4.\n");
        break;

        case SDLK_5:
        change_curr_ifs(state_ptr, 5);
        printf("Displaying image number 5.\n");
        break;

        case SDLK_6:
        change_curr_ifs(state_ptr, 6);
        printf("Displaying image number 6.\n");
        break;

        case SDLK_7:
        change_curr_ifs(state_ptr, 7);
        printf("Displaying image number 7.\n");
        break;

        case SDLK_8:
        change_curr_ifs(state_ptr, 8);
        printf("Displaying image number 8.\n");
        break;

        case SDLK_9:
        change_curr_ifs(state_ptr, 9);
        printf("Displaying image number 9.\n");
        break;

        case SDLK_a:
        state_ptr->aa_on = (0 == state_ptr->aa_on) ? 1 : 0;
        if (0 != state_ptr->aa_on)
        {
            printf("\nAnti-aliasing on\n");
        }
        else
        {
            printf("\nAnti-aliasing off\n");
        }
        break;

        case SDLK_c:
        state_ptr->tweak_mode = TWEAKMODE_CONTRAST;
        printf("\nEntering contrast tweak mode.\nConstrast is %.3f\n", state_ptr->contrast);
        break;

        case SDLK_d:
        print_details(state_ptr);
        break;

        case SDLK_g:
        state_ptr->tweak_mode = TWEAKMODE_GAMMA;
        printf("\nEntering gamma tweak mode.\nGamma is %.3f\n", state_ptr->gamma);
        break;

        case SDLK_h:
        print_help(state_ptr);
        break;

        case SDLK_l:
        print_params(state_ptr);
        break;

        case SDLK_m: // 'm' for 'move'
        state_ptr->tweak_mode = TWEAKMODE_XY;
        printf("\nEntering move mode.\nPosition is (%.3f, %.3f)\n", state_ptr->translation_x, state_ptr->translation_y);
        break;

        case SDLK_p:
        for (i = 0; i < state_ptr->num_gifss; i++)
        {
            if (&state_ptr->ifss[i] != state_ptr->curr_ifs)
            {
                gifs_imitate(&state_ptr->ifss[i], state_ptr->curr_ifs);
                gifs_perturb(&state_ptr->ifss[i]);
            }
            else
            {
                printf("Filling other slots with perturbed versions of %u.\n",
                       (i+state_ptr->num_gifss-1) % state_ptr->num_gifss + 1);
            }
        }
        break;

        case SDLK_q:
        state_ptr->user_req_exit = 1;
        break;

        case SDLK_r:
        if (state_ptr->calc_rate_avail && state_ptr->num_ksecs > 0.0f)
        {
            printf("\nAverage rate is %d points per second.\n", (int)(state_ptr->num_kpts / state_ptr->num_ksecs + 0.5f));
        }
        else
        {
            printf("Sorry, the calculation rate isn't available.\n");
        }
        break;

        case SDLK_s:
        {
            char file_name[64];
            sprintf(file_name, "%u.dat", (unsigned)time(NULL));
            printf("\nSave to %s ", file_name);
            if (save_params(state_ptr, file_name))
            {
                printf("succeeded.\n");
            }
            else
            {
                printf("failed.\n");
            }
            printf("Save to lastsave.dat ", file_name);
            if (save_params(state_ptr, "lastsave.dat"))
            {
                printf("succeeded.\n");
            }
            else
            {
                printf("failed.\n");
            }
        }
        break;

        case SDLK_u:
        if (RENDERMODE_OFF == state_ptr->render_on)
        {
            state_ptr->render_on = RENDERMODE_CONTINUOUS;
            state_ptr->pts_to_next_render = BASE_COMPUTE;
            printf("\nContinuous update is now on.\n");
        }
        else if (RENDERMODE_CONTINUOUS == state_ptr->render_on)
        {
            state_ptr->render_on = RENDERMODE_OFF;
            printf("\nContinuous update is now off.\n");
        }
        break;

        case SDLK_v:
        state_ptr->tweak_mode = TWEAKMODE_VIBRANCE;
        printf("\nEntering vibrance tweak mode.\nVibrance is %.3f\n", state_ptr->vibrance);
        break;

        case SDLK_w:
        {
            char file_name[64];
            sprintf(file_name, "%u.bmp", (unsigned)time(NULL));
            printf("\nWindow capture to %s ", file_name);
            if (0 == SDL_SaveBMP(SDL_GetVideoSurface(), file_name))
            {
                printf("succeeded.\n");
            }
            else
            {
                printf("failed.\n");
            }
        }
        break;

        case SDLK_z: // 'z' for 'zoom'
        state_ptr->tweak_mode = TWEAKMODE_SCALE;
        printf("\nEntering zoom mode.\nZoom is %.3f\n", zoom_level(state_ptr));
        break;

        case SDLK_RETURN:
        if (RENDERMODE_OFF == state_ptr->render_on)
        {
            state_ptr->render_on = RENDERMODE_ONESHOT;
        }
        break;

        case SDLK_UP:
        case SDLK_DOWN:
        case SDLK_LEFT:
        case SDLK_RIGHT:
        handle_tweak(key, state_ptr);
        break;

        case SDLK_ESCAPE:
        state_ptr->tweak_mode = TWEAKMODE_NONE;
        printf("\nExiting tweak mode.\n");
        break;
    }
}


void handle_tweak(SDLKey key, struct engine_state_t* state_ptr)
{
    float sign = 0.0f;

    if (TWEAKMODE_NONE == state_ptr->tweak_mode)
    {
        return;
    }

    switch (key)
    {
        case SDLK_UP:
        case SDLK_RIGHT:
        sign = 1.0f;
        break;

        case SDLK_DOWN:
        case SDLK_LEFT:
        sign = -1.0f;
        break;

        default:
        return;
    }

    switch (state_ptr->tweak_mode)
    {
        case TWEAKMODE_SCALE:
        state_ptr->scale *= (float)exp(log(2.0f) * 0.1f * sign);
        printf("Zoom level is now %.3f.\n", state_ptr->scale / state_ptr->num_rows * 2.0f);
        reset_gfx(state_ptr);
        SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
        break;

        case TWEAKMODE_XY:
        if ((SDLK_LEFT == key) || (SDLK_RIGHT == key))
        {
            state_ptr->translation_x += state_ptr->num_rows / state_ptr->scale * 0.05f * sign;
            printf("X translation is now %.3f.\n", state_ptr->translation_x);
        }
        else
        {
            state_ptr->translation_y += state_ptr->num_rows / state_ptr->scale * 0.05f * sign;
            printf("Y translation is now %.3f.\n", state_ptr->translation_y);
        }
        reset_gfx(state_ptr);
        SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
        break;

        case TWEAKMODE_GAMMA:
        state_ptr->gamma += 0.1f * sign;
        printf("Gamma is now %.3f.\n", state_ptr->gamma);
        SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
        break;

        case TWEAKMODE_VIBRANCE:
        state_ptr->vibrance += 0.1f * sign;
        printf("Vibrance is now %.3f.\n", state_ptr->vibrance);
        SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
        break;

        case TWEAKMODE_CONTRAST:
        state_ptr->contrast *= (float)exp(log(2.0f) * 0.1f * sign);
        printf("Contrast is now %.3f.\n", state_ptr->contrast);
        SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
        break;

        default:
        break;
    }
}


void change_curr_ifs(struct engine_state_t* state_ptr, unsigned new_index)
{
    if ((new_index >= state_ptr->num_gifss) || (&state_ptr->ifss[new_index] == state_ptr->curr_ifs))
    {
        return;
    }

    if (RENDERMODE_OFF == state_ptr->render_on)
    {
        state_ptr->render_on = RENDERMODE_ONESHOT;
    }

    state_ptr->curr_ifs = &state_ptr->ifss[new_index];

    reset_gfx(state_ptr);
    SDL_FillRect(SDL_GetVideoSurface(), NULL, 0);
}


void reset_gfx(struct engine_state_t* state_ptr)
{
    memset(state_ptr->r_plane, 0, state_ptr->num_pixels * sizeof(unsigned));
    memset(state_ptr->g_plane, 0, state_ptr->num_pixels * sizeof(unsigned));
    memset(state_ptr->b_plane, 0, state_ptr->num_pixels * sizeof(unsigned));
    memset(state_ptr->n_plane, 0, state_ptr->num_pixels * sizeof(float));
    memset(s_r_plane, 0, IMAGE_SIZE*IMAGE_SIZE * sizeof(char));
    memset(s_g_plane, 0, IMAGE_SIZE*IMAGE_SIZE * sizeof(char));
    memset(s_b_plane, 0, IMAGE_SIZE*IMAGE_SIZE * sizeof(char));
}


void print_help(struct engine_state_t* state_ptr)
{
    printf("\n");

    printf("\nThe keyboard controls work only when the graphics window is on top.\n");

    printf("\nThere are 10 fractals in memory. Hit digits 1-9 or 0 to view them.\n");

    printf("\nWhen you find one you like, hit 'p' to perturb all of the fractals except for\n");
    printf("the one showing. You'll find that the fractals in the other slots will have\n");
    printf("changed. Maybe you like one of the new ones even more. So hit 'p' to keep that\n");
    printf("one and perturb all of the others. If you don't see any improvements, hit 'p'\n'");
    printf("with the old one showing to get some new possibilities. Keep doing this to see\n");
    printf("what wonderful forms emerge over many iterations.\n");

    printf("\nYou can also adjust gamma, contrast, and 'vibrance' settings. Hit 'g', 'c', or\n");
    printf("'v' to enter the adjustment (i.e. 'tweak') mode of your choice and use the\n");
    printf("arrow keys to change the value.\n");

    printf("\nTo move or zoom the image, hit 'm' or 'z' to enter the desired tweak mode.\n");
    printf("Again, use the arrow keys to make the adjustment.\n");

    printf("\nHit 'l' to get a listing of some of the current parameter values.\n");

    printf("\nHit 'd' to get some information about the displayed fractal.\n");

    printf("\nHit 's' to save the fractal's parameters. The next time Geniter starts, it will\n");
    printf("look for a file called 'lastsave.dat' and load those parameters automatically.\n");

    printf("\nHit 'a' to toggle anti-aliasing.\n");

    printf("\nHit 'u' to toggle continuous graphics updates. When they're on they slow\n");
    printf("fractal generation somewhat.  When they're off, the Enter key will cause the\n");
    printf("graphics to update once.\n");

    printf("\nTo get a window capture, hit 'w'.\n");

    printf("\nHit 'r' to get a measure of the rate of computation.\n");
}


void print_params(struct engine_state_t* state_ptr)
{
    printf("\n");
    printf("Currently showing slot %u.\n",
           (state_ptr->curr_ifs - state_ptr->ifss + state_ptr->num_gifss - 1) % state_ptr->num_gifss + 1);
    printf("Continuous updates are %s.\n", RENDERMODE_CONTINUOUS == state_ptr->render_on ? "on" : "off");
    printf("Anti-aliasing is %s.\n", (0 != state_ptr->aa_on) ? "on" : "off");
    printf("Contrast level is %.3f.\n", state_ptr->contrast);
    printf("Gamma level is %.3f.\n", state_ptr->gamma);
    printf("Zoom level is %.3f.\n", zoom_level(state_ptr));
    printf("Position is (%.3f, %.3f).\n", state_ptr->translation_x, state_ptr->translation_y);
    printf("Vibrance level is %.3f.\n", state_ptr->vibrance);
}


void print_details(struct engine_state_t* state_ptr)
{
    int i;
    char* map_names[MAP_NUMTYPES] =
    { "Identity", "Sine", "Spherical", "Horseshoe", "Polar", "Handkerchief", "Heart", "Disc", "Spiral",
      "Hyperbolic", "Diamond", "Ex", "Julia", "Bent", "Waves", "Fisheye", "Popcorn", "Exponential",
      "Power", "Cosine", "Rings", "Fan", "Blob", "Eyefish", "Bubble", "Cylinder", "Rectangles", "Tangent",
      "Secant" };

    printf("\n");
    for (i = 0; i < state_ptr->curr_ifs->num_functions; i++)
    {
        float* param_set = &state_ptr->curr_ifs->param_sets[i*MAX_PARAMETERS];
        if (0.0f == param_set[0] && 0.0f == param_set[1] && 0.0f == param_set[2])
        {
            if (1.0f == param_set[3] && 0.0f == param_set[4] && 0.0f == param_set[5])
            {
                printf("Identity/");
            }
            else
            {
                printf("Affine/");
            }
        }
        else
        {
            printf("Quadratic/");
        }
        if (0.0f == param_set[6] && 0.0f == param_set[7] && 0.0f == param_set[8])
        {
            if (0.0f == param_set[9] && 1.0f == param_set[10] && 0.0f == param_set[11])
            {
                printf("Identity");
            }
            else
            {
                printf("Affine");
            }
        }
        else
        {
            printf("Quadratic");
        }

        printf(" - %s map - ", map_names[state_ptr->curr_ifs->map_types[i]]);

        if (0.0f == param_set[12] && 0.0f == param_set[13] && 0.0f == param_set[14])
        {
            if (1.0f == param_set[15] && 0.0f == param_set[16] && 0.0f == param_set[17])
            {
                printf("Identity/");
            }
            else
            {
                printf("Affine/");
            }
        }
        else
        {
            printf("Quadratic/");
        }
        if (0.0f == param_set[18] && 0.0f == param_set[19] && 0.0f == param_set[20])
        {
            if (0.0f == param_set[21] && 1.0f == param_set[22] && 0.0f == param_set[23])
            {
                printf("Identity\n");
            }
            else
            {
                printf("Affine\n");
            }
        }
        else
        {
            printf("Quadratic\n");
        }
    }
}


int save_params(struct engine_state_t* state_ptr, char* file_name)
{
    FILE* save_file;
    struct gifs_buffer_t* ifs_buffer;

    if (NULL == file_name)
    {
        return 0;
    }

    save_file = fopen(file_name, "wb");
    if (NULL == save_file)
    {
        return 0;
    }

    ifs_buffer = gifs_export_buffer(state_ptr->curr_ifs);
    if (NULL == ifs_buffer)
    {
        fclose(save_file);
        return 0;
    }

    fwrite(&ifs_buffer->size,         sizeof(ifs_buffer->size),         1, save_file);
    fwrite(&ifs_buffer->data,         ifs_buffer->size,                 1, save_file);
    fwrite(&state_ptr->aa_on,         sizeof(state_ptr->aa_on),         1, save_file);
    fwrite(&state_ptr->contrast,      sizeof(state_ptr->contrast),      1, save_file);
    fwrite(&state_ptr->gamma,         sizeof(state_ptr->gamma),         1, save_file);
    fwrite(&state_ptr->num_columns,   sizeof(state_ptr->num_columns),   1, save_file);
    fwrite(&state_ptr->num_rows,      sizeof(state_ptr->num_rows),      1, save_file);
    fwrite(&state_ptr->render_on,     sizeof(state_ptr->render_on),     1, save_file);
    fwrite(&state_ptr->scale,         sizeof(state_ptr->scale),         1, save_file);
    fwrite(&state_ptr->translation_x, sizeof(state_ptr->translation_x), 1, save_file);
    fwrite(&state_ptr->translation_y, sizeof(state_ptr->translation_y), 1, save_file);
    fwrite(&state_ptr->vibrance,      sizeof(state_ptr->vibrance),      1, save_file);

    fclose(save_file);
    return 1;
}


int load_params(struct engine_state_t* state_ptr, char* file_name)
{
    FILE* load_file;
    unsigned buffer_size = 0;
    char* import_buffer;

    if (NULL == file_name)
    {
        return 0;
    }

    load_file = fopen(file_name, "rb");
    if (NULL == load_file)
    {
        return 0;
    }

    fread(&buffer_size, sizeof(unsigned), 1, load_file);
    import_buffer = (char*)malloc(buffer_size);
    fread(import_buffer, buffer_size, 1, load_file);
    if (0 == gifs_import_buffer(state_ptr->curr_ifs, (void*)import_buffer))
    {
        fclose(load_file);
        return 0;
    }

    if (1 != fread(&state_ptr->aa_on, sizeof(state_ptr->aa_on), 1, load_file))
    {
        fclose(load_file);
        return 1;
    }

    fread(&state_ptr->contrast,      sizeof(state_ptr->contrast),      1, load_file);
    fread(&state_ptr->gamma,         sizeof(state_ptr->gamma),         1, load_file);
    fread(&state_ptr->num_columns,   sizeof(state_ptr->num_columns),   1, load_file);
    fread(&state_ptr->num_rows,      sizeof(state_ptr->num_rows),      1, load_file);
    fread(&state_ptr->render_on,     sizeof(state_ptr->render_on),     1, load_file);
    fread(&state_ptr->scale,         sizeof(state_ptr->scale),         1, load_file);
    fread(&state_ptr->translation_x, sizeof(state_ptr->translation_x), 1, load_file);
    fread(&state_ptr->translation_y, sizeof(state_ptr->translation_y), 1, load_file);
    fread(&state_ptr->vibrance,      sizeof(state_ptr->vibrance),      1, load_file);

    fclose(load_file);
    return 1;
}


float zoom_level(struct engine_state_t* state_ptr)
{
    return state_ptr->scale / 
           ((state_ptr->num_columns > state_ptr->num_rows) ? state_ptr->num_rows : state_ptr->num_columns) * 2.0f;
}
