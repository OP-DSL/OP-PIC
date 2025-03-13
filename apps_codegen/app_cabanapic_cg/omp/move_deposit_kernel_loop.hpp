
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k2 {
enum CellAcc {
    jfx = 0 * 4,
    jfy = 1 * 4,
    jfz = 2 * 4,
};

enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void weight_current_to_accumulator_kernel(
        double* cell_acc,
        const double* q,
        const double& ux, const double& uy, const double& uz,
        const double& dx, const double& dy, const double& dz)
{
    double v0 = 0.0f, v1 = 0.0f, v2 = 0.0f, v3 = 0.0f, v4 = 0.0f, v5 = 0.0f;

    v5 = (*q) * ux * uy * uz * (1.0 / 3.0);

    v4 = (*q)*ux; v1 = v4*dy; v0 = v4-v1; v1 += v4; v4 = 1.0+dz; v2 = v0*v4; v3 = v1*v4; v4 = 1.0-dz; v0 *= v4; v1 *= v4; v0 += v5; v1 -= v5; v2 -= v5; v3 += v5;;
    cell_acc[CellAcc::jfx + 0] += v0;
    cell_acc[CellAcc::jfx + 1] += v1;
    cell_acc[CellAcc::jfx + 2] += v2;
    cell_acc[CellAcc::jfx + 3] += v3;

    v4 = (*q)*uy; v1 = v4*dz; v0 = v4-v1; v1 += v4; v4 = 1.0+dx; v2 = v0*v4; v3 = v1*v4; v4 = 1.0-dx; v0 *= v4; v1 *= v4; v0 += v5; v1 -= v5; v2 -= v5; v3 += v5;;
    cell_acc[CellAcc::jfy + 0] += v0;
    cell_acc[CellAcc::jfy + 1] += v1;
    cell_acc[CellAcc::jfy + 2] += v2;
    cell_acc[CellAcc::jfy + 3] += v3;

    v4 = (*q)*uz; v1 = v4*dx; v0 = v4-v1; v1 += v4; v4 = 1.0+dy; v2 = v0*v4; v3 = v1*v4; v4 = 1.0-dy; v0 *= v4; v1 *= v4; v0 += v5; v1 -= v5; v2 -= v5; v3 += v5;;
    cell_acc[CellAcc::jfz + 0] += v0;
    cell_acc[CellAcc::jfz + 1] += v1;
    cell_acc[CellAcc::jfz + 2] += v2;
    cell_acc[CellAcc::jfz + 3] += v3;
}

enum CellInterp {
    ex = 0,
    dexdy,
    dexdz,
    d2exdydz,
    ey,
    deydz,
    deydx,
    d2eydzdx,
    ez,
    dezdx,
    dezdy,
    d2ezdxdy,
    cbx,
    dcbxdx,
    cby,
    dcbydy,
    cbz,
    dcbzdz,
};

inline void move_deposit_kernel(
    char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    double* part_vel,
    double* part_pos,
    double* part_streak_mid,
    const double* part_weight,
    const double* cell_interp,
    double* cell_acc)
{
    if ((opp_move_hop_iter_one_flag))
    {
        const double& ex       = cell_interp[CellInterp::ex];
        const double& dexdy    = cell_interp[CellInterp::dexdy];
        const double& dexdz    = cell_interp[CellInterp::dexdz];
        const double& d2exdydz = cell_interp[CellInterp::d2exdydz];
        const double& ey       = cell_interp[CellInterp::ey];
        const double& deydz    = cell_interp[CellInterp::deydz];
        const double& deydx    = cell_interp[CellInterp::deydx];
        const double& d2eydzdx = cell_interp[CellInterp::d2eydzdx];
        const double& ez       = cell_interp[CellInterp::ez];
        const double& dezdx    = cell_interp[CellInterp::dezdx];
        const double& dezdy    = cell_interp[CellInterp::dezdy];
        const double& d2ezdxdy = cell_interp[CellInterp::d2ezdxdy];
        double cbx             = cell_interp[CellInterp::cbx];
        const double& dcbxdx   = cell_interp[CellInterp::dcbxdx];
        double cby             = cell_interp[CellInterp::cby];
        const double& dcbydy   = cell_interp[CellInterp::dcbydy];
        double cbz             = cell_interp[CellInterp::cbz];
        const double& dcbzdz   = cell_interp[CellInterp::dcbzdz];

        const double& dx = part_pos[Dim::x];             // Load position
        const double& dy = part_pos[Dim::y];             // Load position
        const double& dz = part_pos[Dim::z];             // Load position

        const double hax  = CONST_qdt_2mc[0] * ( ( ex + dy*dexdy ) + dz * ( dexdz + dy*d2exdydz ) );
        const double hay  = CONST_qdt_2mc[0] * ( ( ey + dz*deydz ) + dx * ( deydx + dz*d2eydzdx ) );
        const double haz  = CONST_qdt_2mc[0] * ( ( ez + dx*dezdx ) + dy * ( dezdy + dx*d2ezdxdy ) );

        cbx  = cbx + dx*dcbxdx;                     // Interpolate B
        cby  = cby + dy*dcbydy;
        cbz  = cbz + dz*dcbzdz;

        double ux = part_vel[Dim::x];             // Load velocity
        double uy = part_vel[Dim::y];             // Load velocity
        double uz = part_vel[Dim::z];             // Load velocity

        ux  += hax;                                 // Half advance E
        uy  += hay;
        uz  += haz;

        double v0   = CONST_qdt_2mc[0]/sqrt(1.0 + (ux*ux + (uy*uy + uz*uz)));
                                                    // Boris - scalars
        double v1   = cbx*cbx + (cby*cby + cbz*cbz);
        double v2   = (v0*v0)*v1;
        double v3   = v0*(1.0+v2*((1.0 / 3.0)+v2*(2.0 / 15.0)));
        double v4   = v3/(1.0+v1*(v3*v3));
        v4  += v4;

        v0   = ux + v3*( uy*cbz - uz*cby );         // Boris - uprime
        v1   = uy + v3*( uz*cbx - ux*cbz );
        v2   = uz + v3*( ux*cby - uy*cbx );
        ux  += v4*( v1*cbz - v2*cby );              // Boris - rotation
        uy  += v4*( v2*cbx - v0*cbz );
        uz  += v4*( v0*cby - v1*cbx );
        ux  += hax;                                 // Half advance E
        uy  += hay;
        uz  += haz;

        part_vel[Dim::x] = ux;                      // save new velocity
        part_vel[Dim::y] = uy;                      // save new velocity
        part_vel[Dim::z] = uz;                      // save new velocity

        /**/                                        // Get norm displacement
        v0   = 1.0/sqrt(1.0 + (ux*ux+ (uy*uy + uz*uz)));
        ux  *= CONST_cdt_d[Dim::x];
        uy  *= CONST_cdt_d[Dim::y];
        uz  *= CONST_cdt_d[Dim::z];

        ux  *= v0;
        uy  *= v0;
        uz  *= v0;

        v0   = dx + ux;                             // Streak midpoint (inbnds)
        v1   = dy + uy;
        v2   = dz + uz;

        v3   = v0 + ux;                             // New position
        v4   = v1 + uy;
        const double v5   = v2 + uz;

        // moving within the cell // Likely
        if (  v3<=1.0 &&  v4<=1.0 &&  v5<=1.0 && -v3<=1.0 && -v4<=1.0 && -v5<=1.0 )
        {
            part_pos[Dim::x] = v3;            // save new position
            part_pos[Dim::y] = v4;            // save new position
            part_pos[Dim::z] = v5;            // save new position

            const double q = part_weight[0] * CONST_qsp[0];

            weight_current_to_accumulator_kernel(
                cell_acc,
                &q,
                ux, uy, uz,
                v0, v1, v2);

            { opp_move_status_flag = OPP_MOVE_DONE; };
            return;
        }
        else
        {
            part_streak_mid[Dim::x] = ux;
            part_streak_mid[Dim::y] = uy;
            part_streak_mid[Dim::z] = uz;
        }
    }

    double s_dir[3];
    double v0, v1, v2, v3;
    int axis, face;

    double s_midx = part_pos[Dim::x]; // Old positions
    double s_midy = part_pos[Dim::y]; // Old positions
    double s_midz = part_pos[Dim::z]; // Old positions

    double s_dispx = part_streak_mid[Dim::x];   // distance moved during push
    double s_dispy = part_streak_mid[Dim::y];   // distance moved during push
    double s_dispz = part_streak_mid[Dim::z];   // distance moved during push

    s_dir[0] = (s_dispx>0) ? 1 : -1; // Try to get the movement direction
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;

    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.
    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx; // 3.4e38 is max of float
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and axis of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    //
    //   axis 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   axis 3        ... streak ends at end of the particle track
    /**/      v3=2,  axis=3;
    if(v0<v3) v3=v0, axis=0;
    if(v1<v3) v3=v1, axis=1;
    if(v2<v3) v3=v2, axis=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak
    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;
    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Accumulate the streak.  Note: accumulator values are 4 times
    // the total physical charge that passed through the appropriate
    // current quadrant in a time-step

    const double q = part_weight[0] * CONST_qsp[0];

    weight_current_to_accumulator_kernel(
        cell_acc,
        &q,
        s_dispx, s_dispy, s_dispz,
        s_midx, s_midy, s_midz);

    // Compute the remaining particle displacment
    part_streak_mid[Dim::x] -= s_dispx;
    part_streak_mid[Dim::y] -= s_dispy;
    part_streak_mid[Dim::z] -= s_dispz;

    // Compute the new particle offset
    part_pos[Dim::x] += (s_dispx + s_dispx);
    part_pos[Dim::y] += (s_dispy + s_dispy);
    part_pos[Dim::z] += (s_dispz + s_dispz);

    // If an end streak, return success (should be ~50% of the time)

    if( axis != 3 )
    {
        // Determine if the particle crossed into a local cell or if it
        // hit a boundary and convert the coordinate system accordingly.
        // Note: Crossing into a local cell should happen ~50% of the
        // time; hitting a boundary is usually a rare event.  Note: the
        // entry / exit coordinate for the particle is guaranteed to be
        // +/-1 _exactly_ for the particle.

        v0 = s_dir[axis];

        // TODO: this conditional could be better
        if (axis == 0) part_pos[Dim::x] = v0;
        if (axis == 1) part_pos[Dim::y] = v0;
        if (axis == 2) part_pos[Dim::z] = v0;

        // _exactly_ on the boundary.
        face = axis;
        if( v0>0 ) face += 3;

        (*opp_p2c) =  opp_c2c[face];

        // TODO: this conditional/branching could be better
        if (axis == 0) { part_pos[Dim::x] = -v0; /* printf("0\n"); */ }
        if (axis == 1) { part_pos[Dim::y] = -v0; /* printf("1\n"); */ }
        if (axis == 2) { part_pos[Dim::z] = -v0; /* printf("2\n"); */ }

        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
    else
    {
        { opp_move_status_flag = OPP_MOVE_DONE; };
    }
}
}

void opp_particle_move__move_deposit_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_vel | OPP_RW
    opp_arg arg1, // p_pos | OPP_RW
    opp_arg arg2, // p_streak_mid | OPP_RW
    opp_arg arg3, // p_weight | OPP_READ
    opp_arg arg4, // c_interp | OPP_READ
    opp_arg arg5 // c_acc | OPP_INC
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_deposit_kernel set_size %d", set->size);

    opp_profiler->start("move_deposit_kernel");

    const int nthreads = omp_get_max_threads();

    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    OPP_mesh_relation_data = (OPP_INT*)p2c_map->p2c_dat->data;

    opp_mpi_halo_exchanges(set, nargs, args);

    opp_create_thread_level_data<OPP_REAL>(args[5], 0);
        
    opp_mpi_halo_wait_all(nargs, args);

#ifdef LOG_HOPS
    std::vector<int> int_hops(nthreads, 0);
    std::vector<int> moreX_hops(nthreads, 0);
    OPP_move_moreX_hops = 0;
#endif

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int n, const int thread, OPP_REAL* arg5_thread_data) {

        OPP_INT *opp_p2c = OPP_mesh_relation_data + n;

        if (opp_p2c[0] == MAX_CELL_INDEX) {
            return;
        }

        OPP_INT *opp_c2c = nullptr;
        char move_flag = OPP_MOVE_DONE;
        bool iter_one_flag = true;

#ifdef LOG_HOPS
        int hops = 0;
#endif
        do {
            move_flag = OPP_MOVE_DONE;
            opp_c2c = c2c_map->map + (opp_p2c[0] * 6);

            opp_k2::move_deposit_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c, 
                (OPP_REAL *)args[0].data + (n * 3),
                (OPP_REAL *)args[1].data + (n * 3),
                (OPP_REAL *)args[2].data + (n * 3),
                (const OPP_REAL *)args[3].data + (n * 1),
                (const OPP_REAL *)args[4].data + (opp_p2c[0] * 18),
                arg5_thread_data + (opp_p2c[0] * 12)
            );
#ifdef LOG_HOPS
            hops++;
#endif
        } while (opp_check_part_move_status(move_flag, iter_one_flag, opp_p2c[0], n, thread));

#ifdef LOG_HOPS
        int_hops[thread] = (int_hops[thread] < hops) ? hops : int_hops[thread];
        if (hops > X_HOPS) moreX_hops[thread]++;
#endif  
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, nargs, args);
    const int total_count = OPP_iter_end - OPP_iter_start;


    opp_profiler->start("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("move_deposit_kernel_only");
    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = OPP_iter_start + ((size_t)total_count * thr) / nthreads;
        const size_t finish = OPP_iter_start + ((size_t)total_count * (thr+1)) / nthreads;
    
        OPP_REAL* arg5_thread_data = (OPP_REAL*)((*(args[5].dat->thread_data))[thr]);

        for (size_t i = start; i < finish; i++)
        {   
            multihop_mover(i, thr, arg5_thread_data);
        }
    }
    opp_profiler->end("move_deposit_kernel_only");

    opp_profiler->end("Mv_AllMv0");


    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {
        
        const std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
        opp_profiler->start(profName);
        
        opp_init_particle_move(set, nargs, args);

        // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
        opp_profiler->start("move_deposit_kernel_only");
        const int iter_count = (OPP_iter_end - OPP_iter_start);
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            const size_t start  = OPP_iter_start + ((size_t)iter_count * thr) / nthreads;
            const size_t finish = OPP_iter_start + ((size_t)iter_count * (thr+1)) / nthreads;
    
            OPP_REAL* arg5_thread_data = (OPP_REAL*)((*(args[5].dat->thread_data))[thr]);

            for (size_t i = start; i < finish; i++)
            {   
                multihop_mover(i, thr, arg5_thread_data);
            }
        }
        opp_profiler->end("move_deposit_kernel_only");

        opp_profiler->end(profName);
    }

    opp_reduce_thread_level_data<OPP_REAL>(args[5]);
#ifdef LOG_HOPS
    OPP_move_max_hops = *std::max_element(int_hops.begin(), int_hops.end());
    OPP_move_moreX_hops = std::accumulate(moreX_hops.begin(), moreX_hops.end(), 0);
#endif

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("move_deposit_kernel");
}
