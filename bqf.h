
#ifndef BQF_INCLUDE_H
#define BQF_INCLUDE_H

#ifndef BQF_CDEC
#ifdef BQF_STATIC
#define BQF_CDEC static
#else
#define BQF_CDEC extern
#endif // BQF_STATIC
#endif // BQF_CDEC

struct bqf_filter;
struct bqf_state_df1;
struct bqf_state_tdf2;

#ifndef bq_sample
#define bq_sample float
#endif

BQF_CDEC void bqf_df1(bq_sample *samples_out,
                      bq_sample const *samples_in,
                      unsigned int samples_count,
                      struct bqf_filter const *filters,
                      struct bqf_state_df1 *states,
                      unsigned int filters_count);

BQF_CDEC void bqf_df1_i(bq_sample *samples,
                        unsigned int samples_count,
                        struct bqf_filter const *filters,
                        struct bqf_state_df1 *states,
                        unsigned int filters_count);

BQF_CDEC void bqf_tdf2(bq_sample *samples_out,
                       bq_sample const *samples_in,
                       unsigned int samples_count,
                       struct bqf_filter const *filters,
                       struct bqf_state_tdf2 *states,
                       unsigned int filters_count);

BQF_CDEC void bqf_tdf2_i(bq_sample *samples,
                         unsigned int samples_count,
                         struct bqf_filter const *filters,
                         struct bqf_state_tdf2 *states,
                         unsigned int filters_count);

#endif // BQF_INCLUDE_H

#ifdef BQF_IMPLEMENTATION

#ifndef BQF_CDEF
#ifdef BQF_STATIC
#define BQF_CDEF static
#else
#define BQF_CDEF
#endif // BQF_STATIC
#endif // BQF_CDEF

struct bqf_filter
{
    bq_sample b0;
    bq_sample b1;
    bq_sample b2;
    bq_sample a1;
    bq_sample a2;
};

struct bqf_state_df1
{
    bq_sample x1;
    bq_sample x2;
    bq_sample y1;
    bq_sample y2;
};

struct bqf_state_tdf2
{
    bq_sample s1;
    bq_sample s2;
};

BQF_CDEF void bqf_df1(bq_sample *samples_out,
                      bq_sample const *samples_in,
                      unsigned int samples_count,
                      struct bqf_filter const *filters,
                      struct bqf_state_df1 *states,
                      unsigned int filters_count)
{
    for (unsigned int f = 0; f < filters_count; f++)
    {
        bq_sample const b0 = filters[f].b0;
        bq_sample const b1 = filters[f].b1;
        bq_sample const b2 = filters[f].b2;
        bq_sample const a1 = filters[f].a1;
        bq_sample const a2 = filters[f].a2;

        bq_sample x1 = states[f].x1;
        bq_sample x2 = states[f].x2;
        bq_sample y1 = states[f].y1;
        bq_sample y2 = states[f].y2;

        for (unsigned int i = 0; i < samples_count; i++)
        {
            bq_sample const x0 = samples_in[i];
            bq_sample const y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            y2 = y1;
            y1 = y0;
            x2 = x1;
            x1 = x0;
            samples_out[i] = y0;
        }

        states[f].x1 = x1;
        states[f].x2 = x2;
        states[f].y1 = y1;
        states[f].y2 = y2;
    }
}

BQF_CDEF void bqf_df1_i(bq_sample *samples,
                        unsigned int samples_count,
                        struct bqf_filter const *filters,
                        struct bqf_state_df1 *states,
                        unsigned int filters_count)
{
    for (unsigned int f = 0; f < filters_count; f++)
    {
        bq_sample const b0 = filters[f].b0;
        bq_sample const b1 = filters[f].b1;
        bq_sample const b2 = filters[f].b2;
        bq_sample const a1 = filters[f].a1;
        bq_sample const a2 = filters[f].a2;

        bq_sample x1 = states[f].x1;
        bq_sample x2 = states[f].x2;
        bq_sample y1 = states[f].y1;
        bq_sample y2 = states[f].y2;

        for (unsigned int i = 0; i < samples_count; i++)
        {
            bq_sample const x0 = samples[i];
            bq_sample const y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            y2 = y1;
            y1 = y0;
            x2 = x1;
            x1 = x0;
            samples[i] = y0;
        }

        states[f].x1 = x1;
        states[f].x2 = x2;
        states[f].y1 = y1;
        states[f].y2 = y2;
    }
}

BQF_CDEF void bqf_tdf2(bq_sample *samples_out,
                       bq_sample const *samples_in,
                       unsigned int samples_count,
                       struct bqf_filter const *filters,
                       struct bqf_state_tdf2 *states,
                       unsigned int filters_count)
{
    for (unsigned int f = 0; f < filters_count; f++)
    {
        bq_sample const b0 = filters[f].b0;
        bq_sample const b1 = filters[f].b1;
        bq_sample const b2 = filters[f].b2;
        bq_sample const a1 = filters[f].a1;
        bq_sample const a2 = filters[f].a2;

        bq_sample s1 = states[f].s1;
        bq_sample s2 = states[f].s2;

        for (unsigned int i = 0; i < samples_count; i++)
        {
            bq_sample const x0 = samples_in[i];
            bq_sample const y0 = b0 * x0 + s1;
            s1 = b1 * x0 - a1 * y0 + s2;
            s2 = b2 * x0 - a2 * y0;
            samples_out[i] = y0;
        }

        states[f].s1 = s1;
        states[f].s2 = s2;
    }
}

BQF_CDEF void bqf_tdf2_i(bq_sample *samples,
                         unsigned int samples_count,
                         struct bqf_filter const *filters,
                         struct bqf_state_tdf2 *states,
                         unsigned int filters_count)
{
    for (unsigned int f = 0; f < filters_count; f++)
    {
        bq_sample const b0 = filters[f].b0;
        bq_sample const b1 = filters[f].b1;
        bq_sample const b2 = filters[f].b2;
        bq_sample const a1 = filters[f].a1;
        bq_sample const a2 = filters[f].a2;

        bq_sample s1 = states[f].s1;
        bq_sample s2 = states[f].s2;

        for (unsigned int i = 0; i < samples_count; i++)
        {
            bq_sample const x0 = samples[i];
            bq_sample const y0 = b0 * x0 + s1;
            s1 = b1 * x0 - a1 * y0 + s2;
            s2 = b2 * x0 - a2 * y0;
            samples[i] = y0;
        }

        states[f].s1 = s1;
        states[f].s2 = s2;
    }
}

#endif // BQF_IMPLEMENTATION
