/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "Stopwatches.hh"


namespace AlgoHex::sw
{

HSW total{"Total time"};
HSW field_opt{"Field generation", total};
HSW projection{"SpH Projection", field_opt};
HSW linsys_setup{"Linear system setup", field_opt};
HSW linsys_setup_smooth{"Smoothness term", linsys_setup};
HSW linsys_setup_normal_constr{"Normal constraints", linsys_setup};
HSW linsys_setup_full_constr{"Full constraints", linsys_setup};
HSW linsys_setup_lin{"Local linearisation", linsys_setup};
HSW linsys_solve{"Linear system solve", field_opt};

HSW sg_extraction{"Singular graph extraction", total};
HSW refinement{"Refine singular region", sg_extraction};
HSW interpolate_cq{"Interpolate for cell quaternions", sg_extraction};
HSW extract_sg{"Extract singular edges", sg_extraction};

HSW lhfg_total{"Locally hexmeshable field generation", total};
HSW singular_graph_repair{"Singular graph repair", lhfg_total};
HSW sg_relocation{"Singularity relocation", lhfg_total};
HSW edge_split{"Split edges", sg_relocation};
HSW edge_collapse{"Collapse edges", sg_relocation};
HSW edge_swap{"Swap edges", sg_relocation};
//    HSW garbage_collect{"Collect garbage", sg_relocation};
HSW singularity_smooth{"Smooth singular graph", sg_relocation};
HSW split_for_dof{"Split for field alignment DOF", sg_relocation};

HSW check_for_termination{"Check for termination", lhfg_total};


HSW sgaf_field_opt_local{"Optimize field locally", lhfg_total};
HSW sgaf_field_opt_global{"Optimize field with alignment", lhfg_total};
HSW sgaf_visualization{"Save mesh for visualization", lhfg_total};

HSW post_remesh{"Post remeshing", total};

HSW frame_field_int_opt{"FrameField Integrability Optimization", total};
HSW frame_field_int_dual_opt{"FrameField Dual Integrability Optimization", total};

HSW parameterization{"Parameterization", total};

HSW hex_extraction{"Hexmesh extraction", total};
} // namespace AlgoHex::sw
