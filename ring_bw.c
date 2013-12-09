/*******************************************************************************
 *  Bandwidth-Latency-Benchmark
 *
 *  Authors: Rolf Rabenseifner
 *           Gerrit Schulz
 *           Michael Speck
 *
 *  Copyright (c) 2003 HLRS, University of Stuttgart
 *
 *******************************************************************************
 *  Additional modifications by Glenn K. Lockwood
 *
 *  San Diego Supercomputer Center, December 2013
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

/* global vars */
double wtick;

#define WTICK_FACTOR 10

/* Message Tags */
#define PING 100
#define PONG 101
#define NEXT_CLIENT 102
#define TO_RIGHT 200
#define TO_LEFT  201

typedef struct {
    int    msglen;
    double ring_lat;
    double ring_bwidth;
    double rand_lat;
    double rand_bwidth;
} BenchmarkResult;

/* measurement results, used only on rank 0 */
static void SumLongLong(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) 
{
    int i, n = *len;
    long long *invecll = (long long *) invec, 
            *inoutvecll = (long long *) inoutvec;
    for (i = n; i; i--, invecll++, inoutvecll++) 
        *inoutvecll += *invecll;
    return;
}

static void cross_ping_pong_set(
    int client_rank_low,
    int client_rank_high,
    int client_rank_stride,
    int server_rank_low,
    int server_rank_high,
    int server_rank_stride,
    int msg_length,
    int loop_length,
    int number_of_measurements,
    int flag,
    double *latency_min,
    double *latency_avg,
    double *latency_max,
    double *bandwidth_min,
    double *bandwidth_avg,
    double *bandwidth_max,
    long long *total_number_of_pairs )
{
    MPI_Status status;
    int    client_rank, server_rank;
    int    i_meas;
    int    i_loop, i;
    unsigned char *sndbuf, *rcvbuf;
    double end_time, start_time, lat_one_meas;
    double *local_results;
    double lat, bw;
    int    result_index;
    long long number_of_results;
    int    size, myrank;
    double loc_latency_min;
    double loc_latency_avg;
    double loc_latency_max;
    double loc_bandwidth_min;
    double loc_bandwidth_avg;
    double loc_bandwidth_max;
    MPI_Op sumll;
    int meas_ok;

    register int base;

    /* get number of processors and own rank */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* check the benchmark parameter */
    if (client_rank_low < 0) client_rank_low = 0;
    if (client_rank_high >= size) client_rank_high = size-1;
    client_rank_high = (client_rank_high-client_rank_low) /
                        client_rank_stride*client_rank_stride + client_rank_low;
    if (server_rank_low < 0) server_rank_low = 0;
    if (server_rank_high >= size) server_rank_high = size-1;
    server_rank_low = server_rank_high -
        (server_rank_high-server_rank_low)/server_rank_stride*server_rank_stride;

    local_results = (double *) malloc( ((server_rank_high -
        server_rank_low)/server_rank_stride+1) * number_of_measurements *
        sizeof(double) );

    /* set the initial result index*/
    result_index = 0;

    /* get memory for the send/recv buffer */
    sndbuf = (unsigned char *) malloc (msg_length);
    rcvbuf = (unsigned char *) malloc (msg_length);

    number_of_results = 0;

    /* do the measurements */
    for (i_meas=0; i_meas < number_of_measurements; i_meas++)
    {
        result_index = 0;
        for (client_rank=client_rank_low; client_rank <= client_rank_high;
            client_rank += client_rank_stride)
        {
            /* the following message receives a token indicating the right to send
            * messages to server processes
            */
            if ((myrank == client_rank) && (client_rank > client_rank_low))
                MPI_Recv (rcvbuf, 0, MPI_BYTE, client_rank - client_rank_stride,
                    NEXT_CLIENT, MPI_COMM_WORLD, &status);

        /* measurement loop */
        for (server_rank = server_rank_low; 
            server_rank <= server_rank_high;
            server_rank += server_rank_stride)
        {
            if (((flag <= 0) && (server_rank > client_rank)) 
            ||  ((flag >= 0) && (server_rank < client_rank)))
            {
                if (myrank == client_rank)
                {
                    do
                    {
                        meas_ok = 0;
                        /* communicate loop_length to server_rank */
                        MPI_Send (&loop_length, 1, MPI_INT,
                            server_rank, PING, MPI_COMM_WORLD);

                        for (i_loop = -1; i_loop < loop_length; i_loop++)
                        {
                            if (i_loop == 0) 
                                start_time = MPI_Wtime ();

                            /* send ping from client_rank to server_rank */
                            base = (i_loop + myrank + 1)&0x7f; /* = mod 128 */
                            sndbuf[0] = base; sndbuf[msg_length-1] = base+1;
                            MPI_Send (sndbuf, msg_length, MPI_BYTE,
                              server_rank, PING, MPI_COMM_WORLD);

                            /* recv pong from server_rank */
                            MPI_Recv (rcvbuf, msg_length, MPI_BYTE,
                                server_rank, PONG, MPI_COMM_WORLD, &status);

                            /* check returned values must be +13 of origin */
                            if (rcvbuf[0] != base+13 || rcvbuf[msg_length-1] != base + 14 )
                            {
                                printf( "[%d]: ERROR: expected %u and %u as first and last byte, but got %u and %u instead\n",
                                    myrank, base+13, base+14,
                                rcvbuf[0], rcvbuf[msg_length-1] ); 
                                fflush(stdout );
                            }
                        }
                        end_time = MPI_Wtime();
                        lat_one_meas = end_time - start_time;

                        if (lat_one_meas < WTICK_FACTOR * wtick)
                        {
                            if (loop_length == 1) loop_length = 2;
                            else loop_length = loop_length * 1.5;
                        }
                        else meas_ok = 1;
                            MPI_Send (&meas_ok, 1, MPI_INT, server_rank, PING, MPI_COMM_WORLD);
                    } while (!meas_ok);

                    fflush(stdout);

                    /* workaround to fix problems with MPI_Wtime granularity */
                    if (!lat_one_meas)
                    {
                        static int complain = 0;
                        lat_one_meas = wtick;
                        if (complain != loop_length)
                        {
#define MSG "In " __FILE__ ", routine bench_lat_bw, the 3rd parameter to cross_ping_pong_controlled was %d; increase it.\n"
                            fprintf(stderr, MSG, loop_length);
                            printf( MSG, loop_length );
#undef MSG
                        }
                        complain = loop_length;
                    }

                    /* store measurement results in the list */
                    local_results [i_meas*number_of_results + result_index] = lat_one_meas / (loop_length*2);
                    result_index++;
                }

                if (myrank == server_rank)
                {
                    do
                    {
                        meas_ok = 0;
                        /* recv the loop_length from client_rank */
                        MPI_Recv (&loop_length, 1, MPI_INT,
                            client_rank, PING, MPI_COMM_WORLD, &status);

                        for (i_loop = -1; i_loop < loop_length; i_loop++)
                        {
                            /* recv ping from client_rank */
                            MPI_Recv (rcvbuf, msg_length, MPI_BYTE,
                                client_rank, PING,
                                MPI_COMM_WORLD, &status);

                            /* server returns received value + const */
                            sndbuf[0] =             rcvbuf[0] + 13;
                            sndbuf[msg_length-1] =  rcvbuf[msg_length-1] + 13;

                            /* send pong from server_rank to client_rank */
                            MPI_Send (sndbuf, msg_length, MPI_BYTE, client_rank, PONG,
                                MPI_COMM_WORLD);
                        }
                        MPI_Recv(&meas_ok, 1, MPI_INT, client_rank, PING,
                            MPI_COMM_WORLD, &status);
                    } while(!meas_ok);
                }
            }
        }

            /* the following message sends a token indicating the right to send
            * messages to server processes
            */
            if (myrank == client_rank
            &&  client_rank < client_rank_high )
                MPI_Send(sndbuf, 0, MPI_BYTE, client_rank + client_rank_stride,
                    NEXT_CLIENT, MPI_COMM_WORLD);
            MPI_Bcast (sndbuf, 0, MPI_BYTE, client_rank_high, MPI_COMM_WORLD);
        }
        number_of_results = result_index;
    }

    /* free the send/recv buffer */
    free(sndbuf);
    free(rcvbuf);

    /* compute local min, max and avg on all client processes */
    /* gather minimal latency for all indexes in first measurement of all measurements */
    for ( i = 0; i < number_of_results; i++ )
        for (i_meas = 1; i_meas < number_of_measurements; i_meas++)
            if ( local_results[i_meas*number_of_results+i] < local_results[i] )
                local_results[i] = local_results[i_meas*number_of_results+i];

    loc_latency_min = 1e99;
    loc_latency_avg = 0;
    loc_latency_max = 0;
    loc_bandwidth_min = 1e99;
    loc_bandwidth_avg = 0;
    loc_bandwidth_max = 0;
    for (i=0; i < number_of_results; i++)
    {
        lat = local_results[i];  bw = msg_length / lat;

        if ( myrank == 0 ) 
        {
            printf ( "[%d] i=%d, lat=%10.6fms, bw=%10.6fMB/s\n", myrank, i, lat*1e3, bw/1e6); 
            fflush( stdout );
        }

        if (lat < (loc_latency_min))  loc_latency_min = lat;
        loc_latency_avg = loc_latency_avg + lat;
        if (lat > (loc_latency_max))  loc_latency_max = lat;
        if (bw < (loc_bandwidth_min))  loc_bandwidth_min = bw;
        loc_bandwidth_avg = loc_bandwidth_avg + bw;
        if (bw > (loc_bandwidth_max))  loc_bandwidth_max = bw;
    }

    if ( myrank == 0 ) 
    {
        printf ( "[%d] Latency   min / avg / max: %10.6f / %10.6f / %10.6f msecs\n",
              myrank, loc_latency_min * 1e3, loc_latency_avg / number_of_results * 1e3, loc_latency_max * 1e3);  
              fflush( stdout );
        printf ( "[%d] Bandwidth min / avg / max: %10.3f / %10.3f / %10.3f MByte/s\n\n",
              myrank, loc_bandwidth_min / 1e6, loc_bandwidth_avg / number_of_results / 1e6, loc_bandwidth_max / 1e6);  
              fflush( stdout );
    }

    /* free the local result list */
    free (local_results);

    /* send all local results to process 0 */
    MPI_Op_create( SumLongLong, 1, &sumll );
    MPI_Reduce (&number_of_results, total_number_of_pairs, 1, MPI_LONG_LONG_INT, sumll, 0,
                MPI_COMM_WORLD);
    MPI_Op_free( &sumll );
    MPI_Reduce (&loc_latency_min, latency_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce (&loc_latency_avg, latency_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce (&loc_latency_max, latency_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce (&loc_bandwidth_min, bandwidth_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce (&loc_bandwidth_avg, bandwidth_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce (&loc_bandwidth_max, bandwidth_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    /* compute global average on process 0 */
    if ((myrank == 0) && (*total_number_of_pairs > 0))
    {
        *latency_avg= *latency_avg / (*total_number_of_pairs);
        *bandwidth_avg= *bandwidth_avg / (*total_number_of_pairs);
    }

    /* print the results */
    if (myrank == 0)
    {
        printf ( "Message Length: %d\n", msg_length);
        printf ( "Latency   min / avg / max: %10.6f / %10.6f / %10.6f msecs\n",
            *latency_min * 1e3, *latency_avg * 1e3, *latency_max * 1e3);
        printf ( "Bandwidth min / avg / max: %10.3f / %10.3f / %10.3f MByte/s\n\n",
            *bandwidth_min / 1e6, *bandwidth_avg / 1e6, *bandwidth_max / 1e6);
        fflush( stdout );
    }
}

static void cross_ping_pong_controlled(
  double max_time,
  int    msg_length,
  int    loop_length,
  int    number_of_measurements,
  double *latency_min,
  double *latency_avg,
  double *latency_max,
  double *bandwidth_min,
  double *bandwidth_avg,
  double *bandwidth_max,
  long long *number_of_pairs
)
{
    int    size, myrank, i;
    double l_dum_min, l_dum_max; /* dummies */
    double b_dum_min, b_dum_avg, b_dum_max; /* dummies */
    long long dum_num_results; /* dummies */
    int    stride;
    double lat_msg;
    int    max_pings, not_prime;
    long long max_pairs;

    /* basic MPI initialization */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    cross_ping_pong_set( 0,0,1,  size-1,size-1,1,
                    msg_length, loop_length, number_of_measurements, 0,
                    &l_dum_min,  &lat_msg,  &l_dum_max,
                    &b_dum_min,  &b_dum_avg,  &b_dum_max, &dum_num_results);

    if ( myrank == 0 ) 
    {
        if (lat_msg*2*(loop_length+1) >= WTICK_FACTOR*wtick)
        {
            max_pairs = max_time / (lat_msg*2*(loop_length+1)*number_of_measurements);
            printf( "MPI_Wtime granularity is ok.\n");
        }
        else
        {
            max_pairs = max_time / (WTICK_FACTOR*wtick*number_of_measurements);
            printf( "Use MPI_Wtick for estimation of max pairs\n");
            fflush( stdout );
        }
        max_pings = (int)sqrt( (double)max_pairs );
        if ( max_pings < 5 ) max_pings = 5;
        stride = 1.0 * size / max_pings + 0.9;
        if ( stride < 1 ) stride = 1;
        if ( stride == 2) stride = 3;
        if ( stride > 3 ) 
        {
            while ( 1 ) 
            {
                not_prime = 0;
                for ( i = 2;  i < stride; i++ )
                {
                    if ( (stride % i) == 0 ) 
                    {
                        not_prime = 1;
                        break;
                    }
                    if ( not_prime )
                    {
                        if ( stride > (size/3) ) 
                            break;
                        else 
                            stride++;
                    }
                    else
                        break;
                }
            }
        }

        printf( "message size:                         %10d\n", msg_length );
        printf( "max time :                            %10.6f secs\n", max_time );
        printf( "latency for msg:                      %10.6f msecs\n", lat_msg*1e3 );
        printf( "estimation for ping pong:             %10.6f msecs\n", lat_msg*2*(loop_length+1)*number_of_measurements*1e3);
        printf( "max number of ping pong pairs       = %10.0f\n", 1.0*max_pairs );
        printf( "max client pings = max server pongs = %10d\n", max_pings );
        printf( "stride for latency                  = %10d\n", stride );
        fflush( stdout );
    }
    MPI_Bcast ( &stride, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cross_ping_pong_set( 0, size-1, stride, 0, size-1, stride,
                    msg_length, loop_length, number_of_measurements, 0,
                    latency_min, latency_avg, latency_max,
                    bandwidth_min, bandwidth_avg, bandwidth_max, number_of_pairs);
    return;
}

static void ring_lat_bw_loop( int msglen, int measurements, int loop_length_proposal,
  int rand_pattern_count, BenchmarkResult *result )
{
    int i_meas, i_pat, i_loop, i, j;
    double start_time, end_time, lat_sendrecv, lat_nonblocking;
    double *latencies; /* measurements * (rand_pattern_count+1) */
    double *max_latencies; /* reduced from all processors with MPI_MAX on rank 0 */
    double avg_latency; /* of random pattern rings */
    int *ranks; /* communication pattern, order of processors */
    int size, myrank, left_rank, right_rank;
    MPI_Request requests[4];
    MPI_Status statuses[4];
    unsigned char *sndbuf_left, *sndbuf_right, *rcvbuf_left, *rcvbuf_right;
    long seedval;
    double rcp = 1.0 / RAND_MAX;
    int loop_length;
    int meas_ok, meas_ok_recv;

    register int base;

    /* get number of processors and own rank */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* alloc memory and init with 0 */
    latencies     = (double *)malloc( measurements * (rand_pattern_count+1) * sizeof( *latencies ) );
    max_latencies = (double *)malloc( measurements * (rand_pattern_count+1) * sizeof( *max_latencies ) );
    ranks = (int *)malloc( size * sizeof( *ranks ) );
    sndbuf_left  = (unsigned char *)malloc( msglen );
    sndbuf_right = (unsigned char *)malloc( msglen );
    rcvbuf_left  = (unsigned char *)malloc( msglen );
    rcvbuf_right = (unsigned char *)malloc( msglen );

    /* init pseudo-random with time seed */
    /*seedval = (long)(time((time_t *) 0));*/
    seedval = (long)time(NULL) ^ getpid();

    if (myrank==0) 
    { 
        printf( "seedval = %ld\n", seedval ); 
        fflush( stdout ); 
    }


    /* benchmark */
    for ( i_meas = 0; i_meas < measurements; i_meas++ ) 
    {
        srand(seedval);
        for ( i_pat = 0; i_pat < rand_pattern_count+1; i_pat++ ) 
        {
            /* build pattern at rank 0 and broadcast to all */
            if ( myrank == 0 ) 
            {
                if ( i_pat > 0 )    /* random pattern */
                { 
                    for (i=0; i<size; i++) 
                        ranks[i] = -1;

                    for (i=0; i<size; i++) 
                    {
                        j = (int)(rand() * rcp * size);
                        while (ranks[j] != -1) j = (j+1) % size;
                        ranks[j] = i;
                    }
                }
                else    /* naturally ordered ring */
                {
                    for (i=0; i<size; i++) 
                        ranks[i] = i;
                }

                if ( i_meas == 0 ) 
                {
                    printf( "i_pat=%3d: ",i_pat);
                    for (i=0; i<size; i++) 
                        printf( " %2d",ranks[i]);
                    printf(  "\n" );  
                    fflush( stdout );
                }
            }
            MPI_Bcast(ranks, size, MPI_INT, 0, MPI_COMM_WORLD);

            /* get rank of left and right partner. therefore find myself (myrank)
             * in pattern first. */
            for ( i = 0; i < size; i++ )
            {
                if ( ranks[i] == myrank ) /* will definitely be found */
                { 
                    left_rank = ranks[(i-1+size)%size];
                    right_rank = ranks[(i+1)%size];
                }
            }

            do
            {
                meas_ok = 0;
                MPI_Allreduce (&loop_length_proposal, &loop_length, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                loop_length_proposal = loop_length;

                /* loop communication */
                for ( i_loop = -1; i_loop < loop_length; i_loop++ ) 
                {
                    if ( i_loop == 0 ) start_time = MPI_Wtime();

                    /* communicate to left and right partner */
                    base = (i_loop + myrank + 1)&0x7f; /* = mod 128 */
                    sndbuf_right[0] = base; sndbuf_right[msglen-1] = base+1;
                    sndbuf_left[0]  = base+2; sndbuf_left[msglen-1]  = base+3;

                    MPI_Sendrecv(
                        sndbuf_right, msglen, MPI_BYTE,
                        right_rank, TO_RIGHT,
                        rcvbuf_left, msglen, MPI_BYTE,
                        left_rank, TO_RIGHT,
                        MPI_COMM_WORLD, &(statuses[0]) );
                    MPI_Sendrecv(
                        sndbuf_left, msglen, MPI_BYTE,
                        left_rank, TO_LEFT,
                        rcvbuf_right, msglen, MPI_BYTE,
                        right_rank, TO_LEFT,
                        MPI_COMM_WORLD, &(statuses[1]) );

                    /* check whether bytes are received correctly */
                    base = (i_loop + left_rank + 1)&0x7f; /* = mod 128 */
                    if ( rcvbuf_left[0] != base || rcvbuf_left[msglen-1] != base+1 )
                    {
                        printf( "[%d]: ERROR: from right: expected %u and %u as first and last byte, but got %u and %u instead\n",
                            myrank, base, base+1,
                            rcvbuf_left[0], rcvbuf_left[msglen-1] ); 
                        fflush( stdout );
                    }
                    base = (i_loop + right_rank + 1)&0x7f; /* = mod 128 */
                    if ( rcvbuf_right[0] != base+2 || rcvbuf_right[msglen-1] != base + 3 )
                    {
                        printf( "[%d]: ERROR: from right: expected %u and %u as first and last byte, but got %u and %u instead\n",
                        myrank, base+2, base+3,
                        rcvbuf_right[0], rcvbuf_right[msglen-1] ); 
                        fflush( stdout );
                    }
                }
                end_time = MPI_Wtime();
                if ((end_time-start_time) < WTICK_FACTOR * wtick)
                {
                    if (loop_length_proposal == 1) loop_length_proposal = 2;
                    else loop_length_proposal = loop_length_proposal * 1.5;
                }
                else 
                    meas_ok=1;
                MPI_Allreduce (&meas_ok, &meas_ok_recv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                meas_ok = meas_ok_recv;
            } while (!meas_ok);

            lat_sendrecv = (end_time-start_time) / (2 * loop_length);

            /* communication loop with non-blocking routines, and previous loop_length */
            for ( i_loop = -1; i_loop < loop_length; i_loop++ ) 
            {
                if ( i_loop == 0 ) start_time = MPI_Wtime();

                /* communicate to left and right partner */
                base = (i_loop + myrank + 1)&0x7f; /* = mod 128 */
                sndbuf_right[0] = base; sndbuf_right[msglen-1] = base+1;
                sndbuf_left[0]  = base+2; sndbuf_left[msglen-1]  = base+3;

                /* irecv left */
                MPI_Irecv(
                    rcvbuf_left, msglen, MPI_BYTE,
                    left_rank, TO_RIGHT,
                    MPI_COMM_WORLD, &requests[0] );
                /* irecv right */
                MPI_Irecv(
                    rcvbuf_right, msglen, MPI_BYTE,
                    right_rank, TO_LEFT,
                    MPI_COMM_WORLD, &requests[1] );
                /* isend right */
                MPI_Isend(
                    sndbuf_right, msglen, MPI_BYTE,
                    right_rank, TO_RIGHT,
                    MPI_COMM_WORLD, &requests[2] );
                /* isend left */
                MPI_Isend(
                    sndbuf_left, msglen, MPI_BYTE,
                    left_rank, TO_LEFT,
                    MPI_COMM_WORLD, &requests[3] );
                /* waitall */
                MPI_Waitall( 4, requests, statuses );

                /* check whether both transfers were done right */
                base = (i_loop + left_rank + 1)&0x7f; /* = mod 128 */
                if ( rcvbuf_left[0] != base || rcvbuf_left[msglen-1] != base+1 )
                {
                    printf( "[%d]: ERROR: from right: expected %u and %u as first and last byte, but got %u and %u instead\n",
                    myrank, base, base+1,
                    rcvbuf_left[0], rcvbuf_left[msglen-1] ); 
                    fflush( stdout );
                }
                base = (i_loop + right_rank + 1)&0x7f; /* = mod 128 */
                if ( rcvbuf_right[0] != base+2 || rcvbuf_right[msglen-1] != base + 3 )
                {
                    printf( "[%d]: ERROR: from right: expected %u and %u as first and last byte, but got %u and %u instead\n",
                    myrank, base+2, base+3,
                    rcvbuf_right[0], rcvbuf_right[msglen-1] ); 
                    fflush( stdout );
                }
            }
            end_time = MPI_Wtime();
            lat_nonblocking = (end_time-start_time) / ( 2 * loop_length );

            /* workaround to fix problems with MPI_Wtime granularity */
            if (!lat_nonblocking)
            {
                static int complain = 0;
                lat_nonblocking = wtick;
                if (complain != loop_length)
                {
                    fprintf( stderr, "In " __FILE__ ", routine bench_lat_bw, the 3rd parameter to ring_lat_bw_loop was %d; increase it.\n", loop_length);
                    printf( "In " __FILE__ ", routine bench_lat_bw, the 3rd parameter to ring_lat_bw_loop was %d; increase it.\n", loop_length);
                }
                complain = loop_length;
            }
            latencies[i_meas*(rand_pattern_count+1)+i_pat] =
                (lat_sendrecv < lat_nonblocking ? lat_sendrecv : lat_nonblocking);
        }
    }

    /* reduce all vectors to get maximum vector at rank 0 */
    MPI_Reduce(
        latencies, max_latencies,
        measurements * (rand_pattern_count+1), MPI_DOUBLE,
        MPI_MAX, 0, MPI_COMM_WORLD );

    /* get minimal measurement from vector as final measurement and compute 
     * latency and bandwidth */
    if ( myrank == 0 ) 
    {
        /* reduce measurements to first minimal measurement */
        for ( i_pat = 0; i_pat < rand_pattern_count+1; i_pat++ )
        {
            /* minimal latencies over all measurements */
            for (i_meas = 1; i_meas < measurements; i_meas++)
            {
                if (max_latencies[i_meas*(rand_pattern_count+1)+i_pat] < max_latencies[i_pat])
                max_latencies[i_pat] = max_latencies[i_meas*(rand_pattern_count+1)+i_pat];
            }
        }

        /* get average latency of random rings by geometric means */
        avg_latency = 0;
        for ( i_pat = 1; i_pat < rand_pattern_count+1; i_pat++ )
            avg_latency += log( max_latencies[i_pat] );
        avg_latency = avg_latency / rand_pattern_count;
        avg_latency = exp( avg_latency );

        /* compute final benchmark results */
        result->msglen = msglen;
        result->ring_lat = max_latencies[0];
        result->ring_bwidth = msglen / max_latencies[0];
        result->rand_lat = avg_latency;
        result->rand_bwidth = msglen / avg_latency;
    }

    /* free memory */
    free( ranks );
    free( latencies );
    free( max_latencies );
    free( sndbuf_left );
    free( sndbuf_right );
    free( rcvbuf_left );
    free( rcvbuf_right );

    if (myrank == 0)
    {
        printf( "Message Size:               %13d Byte\n",   result->msglen );
        printf( "Natural Order Latency:      %13.6f msec\n", result->ring_lat*1e3 );
        printf( "Natural Order Bandwidth:    %13.6f MB/s\n", result->ring_bwidth/1e6 );
        printf( "Avg Random Order Latency:   %13.6f msec\n", result->rand_lat*1e3 );
        printf( "Avg Random Order Bandwidth: %13.6f MB/s\n", result->rand_bwidth/1e6 );
        printf( "\n" );
        fflush( stdout );
    }
    return;
}

static void bench_lat_bw(
    double max_time_for_latency,   /* for ping pong */
    double max_time_for_bandwidth, /* for ping pong */
    int    *msg_length_for_lat,
    int    *msg_length_for_bw,
    double *latency_min, /* */
    double *latency_avg, /* ping pong measurement latency */
    double *latency_max, /* */
    long long *number_of_pairs_for_lat, /* ping pong */
    double *bandwidth_min, /* */
    double *bandwidth_avg, /* ping pong measurement bandwidth */
    double *bandwidth_max, /* */
    long long *number_of_pairs_for_bw, /* ping pong */
    double *ring_lat, /* naturally ordered ring latency */
    double *rand_lat, /* randomly  ordered ring latency */
    double *ring_bw,  /* randomly  ordered ring bandwidth */
    double *rand_bw  /* naturally ordered ring bandwidth */
)
{
    double l_dum_min, l_dum_avg, l_dum_max; /* dummies */
    double b_dum_min, b_dum_avg, b_dum_max; /* dummies */
    BenchmarkResult result_lat, result_bw;
    double wtick_recv;

    int size, myrank;
    double wtime_total, wtime_cross_lat, wtime_cross_bw, wtime_ring_lat, wtime_ring_bw;

    *msg_length_for_lat = 8;
    *msg_length_for_bw  = 2000000;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* get the granularity of MPI_Wtime, but don't trust MPI_Wtick!! */
    wtick = MPI_Wtick();
    if (wtick < 0) wtick = -wtick;

    MPI_Allreduce (&wtick, &wtick_recv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    wtick = wtick_recv;

    if (myrank == 0)
    {
        printf( "MPI_Wtime granularity.\n");
        printf( "Max. MPI_Wtick is %f sec\n", wtick);
    }

    if (wtick < 1e-6) wtick = 1e-6;
    if (wtick > 0.01) wtick = 0.01;

    if (myrank == 0)
    {
        printf( "wtick is set to   %f sec  \n\n", wtick);
        fflush( stdout );
    }

    /* 
     * ping pong 
     */
    wtime_total     = - MPI_Wtime();
    wtime_cross_lat = - MPI_Wtime();


    cross_ping_pong_controlled( 
        max_time_for_latency, *msg_length_for_lat, 8, 5,
        latency_min, latency_avg, latency_max,
        &b_dum_min, &b_dum_avg, &b_dum_max,
        number_of_pairs_for_lat );

    wtime_cross_lat +=  MPI_Wtime();
    wtime_cross_bw  = - MPI_Wtime();

    cross_ping_pong_controlled( 
        max_time_for_bandwidth, *msg_length_for_bw, 1, 2,
        &l_dum_min, &l_dum_avg, &l_dum_max,
        bandwidth_min, bandwidth_avg, bandwidth_max,
        number_of_pairs_for_bw );

    wtime_cross_bw  +=  MPI_Wtime();


    /* 
     * ring
     */
    wtime_ring_lat = - MPI_Wtime();

    ring_lat_bw_loop( *msg_length_for_lat, 8, 5, 30, &result_lat );
    *ring_lat = result_lat.ring_lat;
    *rand_lat = result_lat.rand_lat;

    wtime_ring_lat +=  MPI_Wtime();
    wtime_ring_bw  = - MPI_Wtime();

    ring_lat_bw_loop( *msg_length_for_bw,  3, 2, 10, &result_bw );
    *ring_bw = result_bw.ring_bwidth;
    *rand_bw = result_bw.rand_bwidth;

    wtime_ring_bw  +=  MPI_Wtime();
    wtime_total    +=  MPI_Wtime();

    if (myrank==0)
    { 
        printf( "Execution time (wall clock)      = %9.3f sec on %d processes\n", wtime_total, size);
        printf( " - for cross ping_pong latency   = %9.3f sec\n", wtime_cross_lat);
        printf( " - for cross ping_pong bandwidth = %9.3f sec\n", wtime_cross_bw );
        printf( " - for ring latency              = %9.3f sec\n", wtime_ring_lat);
        printf( " - for ring bandwidth            = %9.3f sec\n", wtime_ring_bw );
        fflush( stdout );
    }
    return;
}

int main(int argc, char **argv) 
{
    int msg_length_for_lat;
    int msg_length_for_bw;
    double ring_lat, rand_lat;
    double ring_bw,  rand_bw;
    int size, myrank;
    double max_time_for_latency;
    double max_time_for_bandwidth;
    double latency_min;
    double latency_avg;
    double latency_max;
    double bandwidth_min;
    double bandwidth_avg;
    double bandwidth_max;
    long long number_of_pairs_for_lat, number_of_pairs_for_bw;
    int myRank, commSize;

    MPI_Init( &argc, &argv );

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0 )
    {
        printf ( "\n------------------------------------------------------------------\n");
        printf ( "Latency-Bandwidth-Benchmark R1.5.1 (c) HLRS, University of Stuttgart\n");
        printf ( "Written by Rolf Rabenseifner, Gerrit Schulz, and Michael Speck, Germany\n\n");
        printf ( "-----------------\n\n");
        fflush( stdout );
    }

  /* The following timings are used for the cross ping pong.
     Additionally, about 300 seconds (on a 100 MB/s) are necessary
     for benchmarking the ring patterns. */
    max_time_for_latency   = 10.0 /*sec*/;
    max_time_for_bandwidth = 30.0 /*sec*/;
    bench_lat_bw( 
        max_time_for_latency, 
        max_time_for_bandwidth,
        &msg_length_for_lat, 
        &msg_length_for_bw,
        &latency_min, 
        &latency_avg, 
        &latency_max,
        &number_of_pairs_for_lat,
        &bandwidth_min, 
        &bandwidth_avg, 
        &bandwidth_max,
        &number_of_pairs_for_bw,
        &ring_lat, 
        &rand_lat, 
        &ring_bw, 
        &rand_bw );

    if (myrank == 0 )
    {
        printf (  "\n------------------------------------------------------------------\n");
        printf (  "Latency-Bandwidth-Benchmark R1.5.1 (c) HLRS, University of Stuttgart\n");
        printf (  "Written by Rolf Rabenseifner, Gerrit Schulz, and Michael Speck, Germany\n\n");

        printf(   "Major Benchmark results:\n" );
        printf(   "------------------------\n\n" );
        printf(   "Max Ping Pong Latency:            %13.6f msecs\n", latency_max*1e3 );
        printf(   "Randomly Ordered Ring Latency:    %13.6f msecs\n", rand_lat*1e3 );
        printf(   "Min Ping Pong Bandwidth:          %13.6f MB/s\n", bandwidth_min/1e6 );
        printf(   "Naturally Ordered Ring Bandwidth: %13.6f MB/s\n", ring_bw/1e6 );
        printf(   "Randomly  Ordered Ring Bandwidth: %13.6f MB/s\n", rand_bw/1e6 );

        printf (  "\n------------------------------------------------------------------\n");

        printf(  "\nDetailed benchmark results:\n" );
        printf(  "Ping Pong:\n" );
        printf ( "Latency   min / avg / max: %10.6f / %10.6f / %10.6f msecs\n",
            latency_min*1e3, 
            latency_avg*1e3, 
            latency_max*1e3 );
        printf (  "Bandwidth min / avg / max: %10.3f / %10.3f / %10.3f MByte/s\n",
            bandwidth_min/1e6, 
            bandwidth_avg/1e6, 
            bandwidth_max/1e6);
        printf(   "Ring:\n" );
        printf(   "On naturally ordered ring: latency= %13.6f msec, bandwidth= %13.6f MB/s\n", ring_lat*1e3, ring_bw/1e6);
        printf(   "On randomly  ordered ring: latency= %13.6f msec, bandwidth= %13.6f MB/s\n", rand_lat*1e3, rand_bw/1e6);
        printf (  "\n------------------------------------------------------------------\n");
        printf(   "\nBenchmark conditions:\n" );
        printf(   " The latency   measurements were done with %8d bytes\n", msg_length_for_lat);
        printf(   " The bandwidth measurements were done with %8d bytes\n", msg_length_for_bw);
        printf(   " The ring communication was done in both directions on %1d processes\n", size);

        printf(   " The Ping Pong measurements were done on \n");
        printf(   "  -  %10.0f pairs of processes for latency benchmarking, and \n", 1.0*number_of_pairs_for_lat);
        printf(   "  -  %10.0f pairs of processes for bandwidth benchmarking, \n", 1.0*number_of_pairs_for_bw);
        printf(   " out of %d*(%d-1) = %10.0f possible combinations on %1d processes.\n", size, size, 1.0*size*(size-1), size);
        printf(   " (1 MB/s = 10**6 byte/sec)\n" );
        printf(   "\n------------------------------------------------------------------\n");
        fflush( stdout );
    }

    MPI_Finalize();

    return 0;
}
