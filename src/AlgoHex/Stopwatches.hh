/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/Util/ScopedStopWatch.hh>

#include <AlgoHex/Config/Export.hh>

namespace AlgoHex::sw
{

using HSW = ::AlgoHex::HierarchicalStopWatch;

extern ALGOHEX_EXPORT HSW total;

extern ALGOHEX_EXPORT HSW field_opt;
extern ALGOHEX_EXPORT HSW projection;
extern ALGOHEX_EXPORT HSW linsys_setup;
extern ALGOHEX_EXPORT HSW linsys_setup_smooth;
extern ALGOHEX_EXPORT HSW linsys_setup_normal_constr;
extern ALGOHEX_EXPORT HSW linsys_setup_full_constr;
extern ALGOHEX_EXPORT HSW linsys_setup_lin;
extern ALGOHEX_EXPORT HSW linsys_solve;

extern ALGOHEX_EXPORT HSW sg_extraction;
extern ALGOHEX_EXPORT HSW refinement;
extern ALGOHEX_EXPORT HSW interpolate_cq;
extern ALGOHEX_EXPORT HSW extract_sg;

extern ALGOHEX_EXPORT HSW lhfg_total;
extern ALGOHEX_EXPORT HSW singular_graph_repair;
extern ALGOHEX_EXPORT HSW sg_relocation;
extern ALGOHEX_EXPORT HSW edge_split;
extern ALGOHEX_EXPORT HSW edge_collapse;
extern ALGOHEX_EXPORT HSW edge_swap;
//    extern ALGOHEX_EXPORT HSW garbage_collect;
extern ALGOHEX_EXPORT HSW singularity_smooth;
extern ALGOHEX_EXPORT HSW split_for_dof;
extern ALGOHEX_EXPORT HSW sgaf_field_opt_local;
extern ALGOHEX_EXPORT HSW sgaf_field_opt_global;
extern ALGOHEX_EXPORT HSW sgaf_visualization;
extern ALGOHEX_EXPORT HSW check_for_termination;

extern ALGOHEX_EXPORT HSW post_remesh;

extern ALGOHEX_EXPORT HSW frame_field_int_opt;
extern ALGOHEX_EXPORT HSW frame_field_int_dual_opt;

extern ALGOHEX_EXPORT HSW parameterization;

extern ALGOHEX_EXPORT HSW hex_extraction;

} // namespace AlgoHex
