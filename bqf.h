
#ifndef BQF_INCLUDE_H
#define BQF_INCLUDE_H

#ifdef BQF_STATIC
#define BQF_API static
#else
#define BQF_API extern
#endif // BQF_STATIC

struct bqf_f32_filter;
struct bqf_f32_state_df1;
struct bqf_f32_state_tdf2;

BQF_API void bqf_f32_df1(float *y, float const *x,
                         struct bqf_f32_filter const *filter,
                         struct bqf_f32_state_df1 *state, unsigned int n);

BQF_API void bqf_f32_df1_i(float *yx, struct bqf_f32_filter const *filter,
                           struct bqf_f32_state_df1 *state, unsigned int n);

BQF_API void bqf_f32_tdf2(float *y, float const *x,
                          struct bqf_f32_filter const *filter,
                          struct bqf_f32_state_tdf2 *state, unsigned int n);

BQF_API void bqf_f32_tdf2_i(float *yx, struct bqf_f32_filter const *filter,
                            struct bqf_f32_state_tdf2 *state, unsigned int n);

#undef BQF_API

#endif // BQF_INCLUDE_H

#ifdef BQF_IMPLEMENTATION

#ifdef BQF_STATIC
#define BQF_IMPL static
#else
#define BQF_IMPL
#endif // BQF_STATIC

struct bqf_f32_filter
{
    float b0;
    float b1;
    float b2;
    float a1;
    float a2;
};

struct bqf_f32_state_df1
{
    float x1;
    float x2;
    float y1;
    float y2;
};

struct bqf_f32_state_tdf2
{
    float s1;
    float s2;
};

BQF_IMPL void bqf_f32_df1(float *const y, float const *const x,
                          struct bqf_f32_filter const *const filter,
                          struct bqf_f32_state_df1 *const state,
                          unsigned int const n)
{
    float const b0 = filter->b0;
    float const b1 = filter->b1;
    float const b2 = filter->b2;
    float const a1 = filter->a1;
    float const a2 = filter->a2;
    float x1 = state->x1;
    float x2 = state->x2;
    float y1 = state->y1;
    float y2 = state->y2;

    for (unsigned int i = 0; i < n; i++)
    {
        float const x0 = x[i];
        float const y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = y0;
        x2 = x1;
        x1 = x0;
        y[i] = y0;
    }

    state->x1 = x1;
    state->x2 = x2;
    state->y1 = y1;
    state->y2 = y2;
}

BQF_IMPL void bqf_f32_df1_i(float *const yx,
                            struct bqf_f32_filter const *const filter,
                            struct bqf_f32_state_df1 *const state,
                            unsigned int const n)
{
    float const b0 = filter->b0;
    float const b1 = filter->b1;
    float const b2 = filter->b2;
    float const a1 = filter->a1;
    float const a2 = filter->a2;
    float x1 = state->x1;
    float x2 = state->x2;
    float y1 = state->y1;
    float y2 = state->y2;

    for (unsigned int i = 0; i < n; i++)
    {
        float const x0 = yx[i];
        float const y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = y0;
        x2 = x1;
        x1 = x0;
        yx[i] = y0;
    }

    state->x1 = x1;
    state->x2 = x2;
    state->y1 = y1;
    state->y2 = y2;
}

BQF_IMPL void bqf_f32_tdf2(float *const y, float const *const x,
                           struct bqf_f32_filter const *const filter,
                           struct bqf_f32_state_tdf2 *const state,
                           unsigned int const n)
{
    float const b0 = filter->b0;
    float const b1 = filter->b1;
    float const b2 = filter->b2;
    float const a1 = filter->a1;
    float const a2 = filter->a2;
    float s1 = state->s1;
    float s2 = state->s2;

    for (unsigned int i = 0; i < n; i++)
    {
        float const x0 = x[i];
        float const y0 = b0 * x0 + s1;
        s1 = b1 * x0 - a1 * y0 + s2;
        s2 = b2 * x0 - a2 * y0;
        y[i] = y0;
    }

    state->s1 = s1;
    state->s2 = s2;
}

BQF_IMPL void bqf_f32_tdf2_i(float *const yx,
                             struct bqf_f32_filter const *const filter,
                             struct bqf_f32_state_tdf2 *const state,
                             unsigned int const n)
{
    float const b0 = filter->b0;
    float const b1 = filter->b1;
    float const b2 = filter->b2;
    float const a1 = filter->a1;
    float const a2 = filter->a2;
    float s1 = state->s1;
    float s2 = state->s2;

    for (unsigned int i = 0; i < n; i++)
    {
        float const x0 = yx[i];
        float const y0 = b0 * x0 + s1;
        s1 = b1 * x0 - a1 * y0 + s2;
        s2 = b2 * x0 - a2 * y0;
        yx[i] = y0;
    }

    state->s1 = s1;
    state->s2 = s2;
}

#undef BQF_IMPL

#endif // BQF_IMPLEMENTATION
