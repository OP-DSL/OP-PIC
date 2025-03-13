
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_dat2_stride = -1;
OPP_INT opp_k2_dat3_stride = -1;
OPP_INT opp_k2_dat4_stride = -1;
OPP_INT opp_k2_dat5_stride = -1;
OPP_INT opp_k2_c2c_map_stride = -1;

OPP_INT* opp_k2_dat0_stride_s = nullptr;
OPP_INT* opp_k2_dat1_stride_s = nullptr;
OPP_INT* opp_k2_dat2_stride_s = nullptr;
OPP_INT* opp_k2_dat3_stride_s = nullptr;
OPP_INT* opp_k2_dat4_stride_s = nullptr;
OPP_INT* opp_k2_dat5_stride_s = nullptr;
OPP_INT* opp_k2_c2c_map_stride_s = nullptr;

namespace opp_k2 {

void opp_dev_move_kernel_sycl(opp_set set, const int nargs, opp_arg *args, opp_map p2c_map, opp_map c2c_map) {

    opp_set_stride(OPP_comm_iteration_d, OPP_comm_iteration_h, OPP_comm_iteration);
    opp_set_stride(opp_k2_dat0_stride_s, opp_k2_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat1_stride_s, opp_k2_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat2_stride_s, opp_k2_dat2_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat3_stride_s, opp_k2_dat3_stride, args[3].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat4_stride_s, opp_k2_dat4_stride, args[4].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat5_stride_s, opp_k2_dat5_stride, args[5].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif
    
    const int num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

    opp_profiler->start("move_kernel_only");

    opp_queue->submit([&](sycl::handler &cgh) {
        
        const OPP_INT* comm_iteration = OPP_comm_iteration_d;
        const OPP_INT* cell_set_size = OPP_cells_set_size_d;

        const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");

        OPP_INT *remove_count = (OPP_INT *)set->particle_remove_count_d;
        OPP_INT *remove_part_indices = (OPP_INT *)OPP_remove_particle_indices_d;
        OPP_INT *move_part_indices = (OPP_INT *)OPP_move_particle_indices_d;
        OPP_INT *move_cell_indices = (OPP_INT *)OPP_move_cell_indices_d;
        OPP_INT *move_count = (OPP_INT *)OPP_move_count_d;

        const OPP_INT* opp_k2_c2c_map_stride_sycl = opp_k2_c2c_map_stride_s;
        const OPP_INT* opp_k2_dat0_stride_sycl = opp_k2_dat0_stride_s;
        const OPP_INT* opp_k2_dat1_stride_sycl = opp_k2_dat1_stride_s;
        const OPP_INT* opp_k2_dat2_stride_sycl = opp_k2_dat2_stride_s;
        const OPP_INT* opp_k2_dat3_stride_sycl = opp_k2_dat3_stride_s;
        const OPP_INT* opp_k2_dat4_stride_sycl = opp_k2_dat4_stride_s;
        const OPP_INT* opp_k2_dat5_stride_sycl = opp_k2_dat5_stride_s;

        const OPP_REAL* CONST_cdt_d_sycl = CONST_cdt_d_s;
        const OPP_REAL* CONST_qdt_2mc_sycl = CONST_qdt_2mc_s;
        const OPP_REAL* CONST_qsp_sycl = CONST_qsp_s;

        OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_vel
        OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_pos
        OPP_REAL* dat2_sycl = (OPP_REAL*)args[2].data_d;     // p_streak_mid
        OPP_REAL* dat3_sycl = (OPP_REAL*)args[3].data_d;     // p_weight
        OPP_REAL* dat4_sycl = (OPP_REAL*)args[4].data_d;     // c_interp
        OPP_REAL** dat5_sycl = (OPP_REAL**)args[5].dat->thread_data_d;     // c_acc
        
        OPP_INT *p2c_map_sycl = (OPP_INT *)p2c_map->p2c_dat->data_d;
        const OPP_INT *c2c_map_sycl = (OPP_INT *)c2c_map->map_d; 

        const OPP_INT iter_start = OPP_iter_start;
        const OPP_INT iter_end = OPP_iter_end; 

        // user provided elemental kernel
        // -----------------------------------------------------------------------------------------
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

        auto  weight_current_to_accumulator_kernel_sycl = [=](
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
        };

        auto  move_deposit_kernel_sycl = [=](
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
                const double& ex       = cell_interp[(CellInterp::ex) * opp_k2_dat4_stride_sycl[0]];
                const double& dexdy    = cell_interp[(CellInterp::dexdy) * opp_k2_dat4_stride_sycl[0]];
                const double& dexdz    = cell_interp[(CellInterp::dexdz) * opp_k2_dat4_stride_sycl[0]];
                const double& d2exdydz = cell_interp[(CellInterp::d2exdydz) * opp_k2_dat4_stride_sycl[0]];
                const double& ey       = cell_interp[(CellInterp::ey) * opp_k2_dat4_stride_sycl[0]];
                const double& deydz    = cell_interp[(CellInterp::deydz) * opp_k2_dat4_stride_sycl[0]];
                const double& deydx    = cell_interp[(CellInterp::deydx) * opp_k2_dat4_stride_sycl[0]];
                const double& d2eydzdx = cell_interp[(CellInterp::d2eydzdx) * opp_k2_dat4_stride_sycl[0]];
                const double& ez       = cell_interp[(CellInterp::ez) * opp_k2_dat4_stride_sycl[0]];
                const double& dezdx    = cell_interp[(CellInterp::dezdx) * opp_k2_dat4_stride_sycl[0]];
                const double& dezdy    = cell_interp[(CellInterp::dezdy) * opp_k2_dat4_stride_sycl[0]];
                const double& d2ezdxdy = cell_interp[(CellInterp::d2ezdxdy) * opp_k2_dat4_stride_sycl[0]];
                double cbx             = cell_interp[(CellInterp::cbx) * opp_k2_dat4_stride_sycl[0]];
                const double& dcbxdx   = cell_interp[(CellInterp::dcbxdx) * opp_k2_dat4_stride_sycl[0]];
                double cby             = cell_interp[(CellInterp::cby) * opp_k2_dat4_stride_sycl[0]];
                const double& dcbydy   = cell_interp[(CellInterp::dcbydy) * opp_k2_dat4_stride_sycl[0]];
                double cbz             = cell_interp[(CellInterp::cbz) * opp_k2_dat4_stride_sycl[0]];
                const double& dcbzdz   = cell_interp[(CellInterp::dcbzdz) * opp_k2_dat4_stride_sycl[0]];

                const double& dx = part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]];             // Load position
                const double& dy = part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]];             // Load position
                const double& dz = part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]];             // Load position

                const double hax  = CONST_qdt_2mc_sycl[0] * ( ( ex + dy*dexdy ) + dz * ( dexdz + dy*d2exdydz ) );
                const double hay  = CONST_qdt_2mc_sycl[0] * ( ( ey + dz*deydz ) + dx * ( deydx + dz*d2eydzdx ) );
                const double haz  = CONST_qdt_2mc_sycl[0] * ( ( ez + dx*dezdx ) + dy * ( dezdy + dx*d2ezdxdy ) );

                cbx  = cbx + dx*dcbxdx;                     // Interpolate B
                cby  = cby + dy*dcbydy;
                cbz  = cbz + dz*dcbzdz;

                double ux = part_vel[(Dim::x) * opp_k2_dat0_stride_sycl[0]];             // Load velocity
                double uy = part_vel[(Dim::y) * opp_k2_dat0_stride_sycl[0]];             // Load velocity
                double uz = part_vel[(Dim::z) * opp_k2_dat0_stride_sycl[0]];             // Load velocity

                ux  += hax;                                 // Half advance E
                uy  += hay;
                uz  += haz;

                double v0   = CONST_qdt_2mc_sycl[0]/sqrt(1.0 + (ux*ux + (uy*uy + uz*uz)));
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

                part_vel[(Dim::x) * opp_k2_dat0_stride_sycl[0]] = ux;                      // save new velocity
                part_vel[(Dim::y) * opp_k2_dat0_stride_sycl[0]] = uy;                      // save new velocity
                part_vel[(Dim::z) * opp_k2_dat0_stride_sycl[0]] = uz;                      // save new velocity

                /**/                                        // Get norm displacement
                v0   = 1.0/sqrt(1.0 + (ux*ux+ (uy*uy + uz*uz)));
                ux  *= CONST_cdt_d_sycl[Dim::x];
                uy  *= CONST_cdt_d_sycl[Dim::y];
                uz  *= CONST_cdt_d_sycl[Dim::z];

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
                    part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]] = v3;            // save new position
                    part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]] = v4;            // save new position
                    part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]] = v5;            // save new position

                    const double q = part_weight[(0) * opp_k2_dat3_stride_sycl[0]] * CONST_qsp_sycl[0];

                    weight_current_to_accumulator_kernel_sycl(
                        cell_acc,
                        &q,
                        ux, uy, uz,
                        v0, v1, v2);

                    { opp_move_status_flag = OPP_MOVE_DONE; };
                    return;
                }
                else
                {
                    part_streak_mid[(Dim::x) * opp_k2_dat2_stride_sycl[0]] = ux;
                    part_streak_mid[(Dim::y) * opp_k2_dat2_stride_sycl[0]] = uy;
                    part_streak_mid[(Dim::z) * opp_k2_dat2_stride_sycl[0]] = uz;
                }
            }

            double s_dir[3];
            double v0, v1, v2, v3;
            int axis, face;

            double s_midx = part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]]; // Old positions
            double s_midy = part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]]; // Old positions
            double s_midz = part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]]; // Old positions

            double s_dispx = part_streak_mid[(Dim::x) * opp_k2_dat2_stride_sycl[0]];   // distance moved during push
            double s_dispy = part_streak_mid[(Dim::y) * opp_k2_dat2_stride_sycl[0]];   // distance moved during push
            double s_dispz = part_streak_mid[(Dim::z) * opp_k2_dat2_stride_sycl[0]];   // distance moved during push

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

            const double q = part_weight[(0) * opp_k2_dat3_stride_sycl[0]] * CONST_qsp_sycl[0];

            weight_current_to_accumulator_kernel_sycl(
                cell_acc,
                &q,
                s_dispx, s_dispy, s_dispz,
                s_midx, s_midy, s_midz);

            // Compute the remaining particle displacment
            part_streak_mid[(Dim::x) * opp_k2_dat2_stride_sycl[0]] -= s_dispx;
            part_streak_mid[(Dim::y) * opp_k2_dat2_stride_sycl[0]] -= s_dispy;
            part_streak_mid[(Dim::z) * opp_k2_dat2_stride_sycl[0]] -= s_dispz;

            // Compute the new particle offset
            part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]] += (s_dispx + s_dispx);
            part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]] += (s_dispy + s_dispy);
            part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]] += (s_dispz + s_dispz);

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
                if (axis == 0) part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]] = v0;
                if (axis == 1) part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]] = v0;
                if (axis == 2) part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]] = v0;

                // _exactly_ on the boundary.
                face = axis;
                if( v0>0 ) face += 3;

                (*opp_p2c) =  opp_c2c[(face) * opp_k2_c2c_map_stride_sycl[0]];

                // TODO: this conditional/branching could be better
                if (axis == 0) { part_pos[(Dim::x) * opp_k2_dat1_stride_sycl[0]] = -v0; /* printf("0\n"); */ }
                if (axis == 1) { part_pos[(Dim::y) * opp_k2_dat1_stride_sycl[0]] = -v0; /* printf("1\n"); */ }
                if (axis == 2) { part_pos[(Dim::z) * opp_k2_dat1_stride_sycl[0]] = -v0; /* printf("2\n"); */ }

                { opp_move_status_flag = OPP_NEED_MOVE; };
            }
            else
            {
                { opp_move_status_flag = OPP_MOVE_DONE; };
            }
        };

        // -----------------------------------------------------------------------------------------
        auto opp_move_kernel = [=](sycl::nd_item<1> item) {
            
            const int tid = item.get_global_linear_id();
            const int n = tid + iter_start;

            if (n < iter_end) {
                OPP_INT *opp_p2c = (p2c_map_sycl + n);
                if (opp_p2c[0] == MAX_CELL_INDEX) {
                    return;
                }
                
                char move_flag = OPP_NEED_MOVE;
                bool iter_one_flag = (comm_iteration[0] > 0) ? false : true;

                OPP_REAL arg5_p2c_local[12];
                OPP_REAL* tmp5_sycl = dat5_sycl[item.get_local_id(0) % array_count];

                do {
                    const OPP_INT p2c = opp_p2c[0];
                    const OPP_INT* opp_c2c = c2c_map_sycl + p2c;

                    for (int d = 0; d < 12; ++d)
                        arg5_p2c_local[d] = OPP_REAL_ZERO;

                    move_deposit_kernel_sycl(
                        move_flag, iter_one_flag, opp_c2c, opp_p2c,
                        dat0_sycl + n, // p_vel 
                        dat1_sycl + n, // p_pos 
                        dat2_sycl + n, // p_streak_mid 
                        dat3_sycl + n, // p_weight 
                        dat4_sycl + p2c, // c_interp 
                        arg5_p2c_local // c_acc 
         
                    ); 
                    
                    for (int d = 0; d < 12; ++d)
                        opp_atomic_fetch_add(tmp5_sycl + p2c + (d * opp_k2_dat5_stride_sycl[0]), arg5_p2c_local[d]);
                
            
                } while (opp_part_check_status_device(move_flag, iter_one_flag, opp_p2c, n, 
                            remove_count, remove_part_indices, 
                            move_part_indices, move_cell_indices, move_count, cell_set_size));
            }        
        };

        // -----------------------------------------------------------------------------------------
        cgh.parallel_for<class opp_particle_move>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), opp_move_kernel);
    });

    OPP_DEVICE_SYNCHRONIZE();

    opp_profiler->end("move_kernel_only");
}
} // end of namespace

//--------------------------------------------------------------
void opp_particle_move__move_deposit_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_vel | OPP_RW
    opp_arg arg1,   // p_pos | OPP_RW
    opp_arg arg2,   // p_streak_mid | OPP_RW
    opp_arg arg3,   // p_weight | OPP_READ
    opp_arg arg4,   // c_interp | OPP_READ
    opp_arg arg5   // c_acc | OPP_INC
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_deposit_kernel set_size %d", set->size);

    opp_profiler->start("move_deposit_kernel");

    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;

    opp_set_stride(OPP_cells_set_size_d, OPP_cells_set_size, set->cells_set->size);
    opp_set_stride(opp_k2_c2c_map_stride_s, opp_k2_c2c_map_stride, c2c_stride);

    opp_mpi_halo_wait_all(nargs, args);

    opp_init_particle_move(set, nargs, args);

    opp_create_thread_level_data<OPP_REAL>(args[5]);

    opp_k2::opp_dev_move_kernel_sycl(set, nargs, args, p2c_map, c2c_map);

    while (opp_finalize_particle_move(set)) {
        
        opp_init_particle_move(set, nargs, args);

        opp_k2::opp_dev_move_kernel_sycl(set, nargs, args, p2c_map, c2c_map);
    }
    opp_reduce_thread_level_data<OPP_REAL>(args[5]);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();
 
    opp_profiler->end("move_deposit_kernel");
}
